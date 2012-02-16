# Generator base classes: Generator, ctsGen, discGen
from __future__ import division

from allimports import *
from PyDSTool.utils import *
from PyDSTool.common import *
from PyDSTool.parseUtils import symbolMapClass
import PyDSTool.Events as Events

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, transpose, shape
import math, random
from copy import copy, deepcopy
try:
    # use pscyo JIT byte-compiler optimization, if available
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

# -----------------------------------------------------------------------------
# Version History:
#
# 21 June 2006
#   Added _addEvents and _makeBoundsEvents methods.
#
# 28 January 2006
#   Added __deepcopy__ method.
#
# 22 January 2006
#   Added trajevents structure's place holder in base class.
#   Altered behaviour of getEvents to return point set of event time + variable
#     data, and renamed previous version to getEventTimes.
#
# 30 September 2005
#   Added info method.
#   Added GenSpecHelper class.
#   Added methods to return spec, auxspec, auxfn and event definitions
#     (for python-based FuncSpecs)
#
# 12 July 2005
#   Added error E_NONUNIQUETERM and associated error message.
#
# -----------------------------------------------------------------------------

__all__ = ['ctsGen', 'discGen', 'theGenSpecHelper', 'Generator',
           'genDB']

# -----------------------------------------------------------------------------

class genDBClass(object):
    """This class keeps a tab of which non-Python Generators have been
    created in a session. A single global instance of this class is created,
    and prevents the user from re-using an instance of a low-level
    DLL (such as that for a Dopri_ODEsystem vector field) in more than
    one Generator."""

    def __init__(self):
        self.database = {}

    def check(self, gen):
        """Look up generator instance and verify its type match"""
        if gen.funcspec.targetlang == 'python':
            # these do not need to be checked
            return
        try:
            if gen._solver.rhs in self.database:
                return className(gen) == self.database[gen._solver.rhs][1]
            else:
                return True
        except AttributeError:
            # these do not need to be checked
            return

    def unregister(self, gen):
        if gen.funcspec.targetlang == 'python':
            # these do not need to be checked
            return
        try:
            del self.database[gen._solver.rhs]
        except (AttributeError, KeyError):
            # invalid type of generator or not present
            # so do nothing
            return

    def register(self, gen):
        if gen.funcspec.targetlang == 'python':
            # these do not need to be checked
            return
        # verify name of vector field in database
        try:
            if gen._solver.rhs not in self.database:
                self.database[gen._solver.rhs] = (gen.name, className(gen))
            else:
                raise PyDSTool_KeyError, "Generator name %s already exists"%gen.name
        except AttributeError:
            # these do not need to be checked
            return

    def __repr__(self):
        return "Generator internal database class"

    __str__ = __repr__


    def clearall(self):
        self.database = {}


# use single instance of nameResolver per session
global genDB
genDB = genDBClass()


# ----------------------------------------------------------------------


class Generator(object):
    """
    Trajectory Generator abstract class.
    """

    def __init__(self, kw):
        # _funcreg stores function instances (and base classes) that are
        # passed from self.funcspec, as they are not defined in __main__.
        # A Generator object cannot be deep copied or
        # pickled without the additional __getstate__ and __setstate__
        # methods which know how to reconstruct these functions.
        self._funcreg = {}
        # initialization keyword keys
        self.needKeys = ['name']
        self.optionalKeys = ['globalt0','checklevel','model','abseps']
        try:
            # sometimes certain keys are checked prior to calling this base
            # class so we don't want to tread on the feet of those __init__
            # methods
            dummy = self.foundKeys
            # if this works then we won't reset to zero
        except:
            self.foundKeys = 0
        try:
            self.name = kw['name']
            self.foundKeys += 1
        except KeyError:
            raise PyDSTool_KeyError, "name must be supplied as a keyword in arguments"
        # dimension value is set later
        self.dimension = None
        # Regular Variable objects, and auxiliary variables.
        # These are callable and can even appear in
        # inputs. They are defined internally by the FuncSpec object.
        # These self.variables are generally placeholders containing domain
        # information only, except for the special Generator types LookupTable
        # and InterpTable. The actual variable contents from a computed
        # trajectory are created locally during computeTraj and immediately
        # exported to a Trajectory object, not affecting these variable objects.
        self.variables = {}
        # initial conditions for each regular and auxiliary variable
        self.initialconditions = {}
        # Generator's functional specification for variables, e.g.
        # right-hand side of ODE, or data lists for look-up table
        self.funcspec = None
        # Generator's internal pars - for the functional specification
        self.pars = {}
        # Place-holder for event structure, if used
        self.eventstruct = None
        # Place-holder for trajectory events, if used
        self.trajevents = None
        # Autonomous external inputs (dictionary of callable Variables)
        self.inputs = {}
        # Local independent variable interval and validity range (relative to global t0)
        self.indepvariable = None
        # Sorted list of callable names
        self._callnames = []
        # Registry of object names and types
        self._registry = {}
        # Absolute tolerance for Interval endpoints
        if 'abseps' in kw:
            self._abseps = kw['abseps']
            self.foundKeys += 1
        else:
            self._abseps = 1e-14
        # `Checklevel` code determining response to uncertain float comparisons
        if 'checklevel' in kw:
            if kw['checklevel'] not in range(4):
                    raise ValueError('Invalid interval endpoint checking option')
            self.checklevel = kw['checklevel']
            self.foundKeys += 1
        else:
            # default to level 2 (warn on uncertain, error on not contained)
            self.checklevel = 2
        # Warnings, e.g. for bounds violations not causing exceptions, event
        # flags
        self.warnings = []
        self._warnmessages = warnmessages
        self._warnfields = warnfields

        # Errors, e.g. failed numerical integration
        self.errors = []
        self._errmessages = {E_NONUNIQUETERM: 'More than one terminal event found',
                             E_COMPUTFAIL: 'Computation of trajectory failed'}
        self._errorfields = {E_NONUNIQUETERM: ['t', 'event list'],
                             E_COMPUTFAIL: ['t', 'error info']}
        # traceback may store information about variable state, pars, etc.
        # at time of an error that breaks the solver
        self.traceback = {}
        # global independent variable reference (for non-autonomous systems,
        # especially when a Generator is embedded in a hybrid Model object)
        if 'globalt0' in kw:
            assert isinstance(kw['globalt0'], float) or \
                isinstance(kw['globalt0'], int), 'Incorrect type of globalt0'
            self.globalt0 = kw['globalt0']
            self.foundKeys += 1
        else:
            self.globalt0 = 0.
        # If part of a Model, keep a reference to which one this
        # Generator belongs to
        if 'model' in kw:
            assert isinstance(kw['model'], str), "model tag must be a string"
            self._modeltag = kw['model']
            self.foundKeys += 1
        else:
            self._modeltag = None
        # Indicator of whether a trajectory has been successfully computed
        self.defined = False


    def getEventTimes(self, globalindepvar=True):
        """Produce dictionary of lists of all flagged events' independent
        variable values, for each event (whether terminal or not).

        The events are not guaranteed to be ordered by the value of the
        independent variable.
        """

        result = {}.fromkeys(self.eventstruct.events.keys(), None)
        for evname in self.eventstruct.events.keys():
            result[evname] = []
        if len(self.warnings) > 0:
            if globalindepvar:
                for (w, d) in self.warnings:
                    if w == W_TERMEVENT or w == W_NONTERMEVENT:
                        if isinstance(d[1], list):
                            for evname in d[1]:
                                result[evname].append(d[0]+self.globalt0)
                        else:
                            result[d[1]].append(d[0]+self.globalt0)
            else:
                for (w, d) in self.warnings:
                    if w == W_TERMEVENT or w == W_NONTERMEVENT:
                        if isinstance(d[1], list):
                            for evname in d[1]:
                                result[evname].append(d[0])
                        else:
                            result[d[1]].append(d[0])
        for evname in result:
            result[evname].sort()
        return result


    def getEvents(self, globalindepvar=True):
        """Produce dictionary of pointlists of all flagged events' independent
        and dependent variable values, for each event (whether terminal or not).

        The events are not guaranteed to be ordered by the value of the
        independent variable.
        """

        if globalindepvar and self.globalt0 != 0:
            result = {}
            for (evname, evptset) in self.trajevents.iteritems():
                result[evname] = copy(evptset)
                try:
                    result[evname].indepvararray += self.globalt0
                except AttributeError:
                    # empty pointset
                    pass
            return result
        else:
            return copy(self.trajevents)


    def haveJacobian(self):
        """Default method. Can be overridden by subclasses."""
        return False


    def haveJacobian_pars(self):
        """Default method. Can be overridden by subclasses."""
        return False


    def info(self, verboselevel=1):
        print self._infostr(verboselevel)


    def _infostr(self, verbose=1):
        """Return detailed information about the Generator
        specification."""

        if verbose == 0:
            outputStr = "Generator "+self.name
        else:
            outputStr = '**************************************************'
            outputStr += '\n           Generator  '+self.name
            outputStr +='\n**************************************************'
            outputStr +='\nType : ' + className(self)
            outputStr +='\nIndependent variable interval: ' + str(self.indepvariable.depdomain)
            outputStr +='\nGlobal t0 = ' + str(self.globalt0)
            outputStr +='\nInterval endpoint check level = ' + str(self.checklevel)
            outputStr +='\nDimension = ' + str(self.dimension)
            outputStr +='\n'
            if isinstance(self.funcspec, FuncSpec):
                outputStr += self.funcspec._infostr(verbose)
        if verbose == 2:
            outputStr += '\nVariables` validity intervals:'
            for v in self.variables.values():
                outputStr += '\n  ' + str(v.depdomain)
            if self.eventstruct is not None:
                outputStr += '\nEvents defined:'
                outputStr += '\n  ' + str(self.eventstruct.events.keys())
        if self._modeltag is not None:
            outputStr += '\nAssociated Model: ' + self._modeltag.name
        if verbose > 0:
            outputStr += '\n'
        return outputStr


    def showEventSpec(self):
        for evname, ev in self.eventstruct.events.iteritems():
            print evname + ":\n" + ev._funcstr
            print "\n"


    def showSpec(self):
        print self.funcspec.spec[0]


    def showAuxSpec(self):
        print self.funcspec.auxspec[0]


    def showAuxFnSpec(self, auxfnname=None):
        if auxfnname is None:
            retdict = {}
            for aname, aspec in self.funcspec.auxfns.iteritems():
                retdict[aname] = aspec[0]
            info(retdict)
        else:
            try:
                print self.funcspec.auxfns[auxfnname][0]
            except KeyError:
                raise NameError("Aux function %s not found"%auxfnname)


    def __repr__(self):
        return self._infostr(verbose=0)


    __str__ = __repr__


    def validateSpec(self):
        try:
            assert self.dimension > 0
            assert len(self.variables) == self.dimension
            # don't assert self.pars because not all systems need them
            assert self.indepvariable.name == 't'
            assert self.funcspec
            assert self.checklevel in range(4)
            #  check that all names in individual dicts are all in _registry
            if self.pars:
                for name in self.pars:
                    t = type(self.pars[name])
                    assert t in [Int, Float]
                    assert t == self._registry[name]
            if self.inputs:
                for subjectname, obj in self.inputs.iteritems():
                    # test for containment of input's interval in independent
                    # variable interval
                    # (use checklevel = 1 for this o/w could get errors)
                    assert self.contains(obj.indepdomain,
                                         self.indepvariable.indepdomain, 0)
                    # test that types entered in registry are still correct
                    assert type(obj) == self._registry[subjectname]
            for name in self.variables:
                assert self.variables[name].__class__ == self._registry[name]
            dummy = self.indepvariable(self._tdata[0]) # exception if this call is ill-defined
            # check consistency with FuncSpec type of self.funcspec
            # (unnecessary for dictionary version of FuncSpec)
            if isinstance(self.funcspec, FuncSpec):
                varnames = self.variables.keys()
                fsvars = self.funcspec.vars
                if len(varnames) > 1:
                    varnames.sort()
                    fsvars.sort()
                    assert varnames == fsvars, ('Inconsistency with funcspec '
                                                'variable names')
                else:
                    assert varnames == fsvars
                parnames = self.pars.keys()
                fspars = self.funcspec.pars
                if len(parnames) > 1:
                    parnames.sort()
                    fspars.sort()
                    assert parnames == fspars, ('Inconsistency with funcspec '
                                                'parameter names')
                else:
                    assert parnames == fspars
                if self.inputs:
                    inputnames = self.inputs.keys()
                    fsinputs = self.funcspec.inputs
                    if len(inputnames) > 1:
                        inputnames.sort()
                        fsinputs.sort()
                        assert inputnames == fsinputs, ('Inconsistency with funcspec'
                                                        ' input names')
                    else:
                        assert inputnames == fsinputs
            else:
                assert len(self.funcspec) == self.dimension
        except:
            print 'Invalid system specification'
            raise


    # call this after all expected keywords have been processed
    def checkArgs(self, kw):
        if len(kw) == self.foundKeys:
            for name in self.needKeys:
                if name not in kw:
                    raise PyDSTool_KeyError('Necessary key missing: ' + name)
            for name in kw:
                if name not in self.needKeys + self.optionalKeys:
                    raise PyDSTool_KeyError('Key name ' + name + ' is invalid')
        else:
            print 'Keywords supplied:\n\t' + str(kw.keys())
            print '# keywords found: ' + str(self.foundKeys)
            print 'Needed:\n\t' + str(self.needKeys)
            print 'Optional:\n\t' + str(self.optionalKeys)
            raise PyDSTool_KeyError, 'Invalid keyword arguments for this class'


    def _register(self, items):
        """_register names and types of sub-system variables (including
        Generator variables), pars and external inputs.

        Names must be unique for the Generator.
        """

        if isinstance(items, dict):
            # for parameter and variable dictionaries
            for name, v in items.iteritems():
                if isinstance(self.funcspec, FuncSpec):
                    assert name in self.funcspec.vars \
                       or name in self.funcspec.auxvars \
                       or name in self.funcspec.pars \
                       or name in self.funcspec.inputs, \
                       ("Generator = '"
                       +self.name+"': name "+name+" not found in "
                       "functional specification declaration")
                if name not in self._registry:
                    self._registry[name] = type(v)
                else:
                    raise ValueError, ('The name `' + name + '` of type `'
                                       + type(v).__name__ +
                                       '` already exists in the registry')
                if isinstance(v, Variable) and \
                   name in self.variables or name in self.inputs:
                    self._callnames.append(name)
            self._callnames.sort()
        elif isinstance(items, Variable) and items.name == 't':
            # for self.indepvariable
            if items.name not in self._registry:
                self._registry[items.name] = type(items)
            else:
                raise ValueError, ('The reserved name `t` has already'
                                   ' been declared to the registry')
        else:
            raise TypeError, ('Expected dictionary or independent variable in '
                              'argument to _register()')


    def _addEvents(self, evs):
        if isinstance(evs, list):
            map(self.eventstruct.add, copy(evs))
        elif isinstance(evs, Event):
            self.eventstruct.add(evs)  # singleton
        else:
            raise TypeError, ('Unsupported type of argument for event '
                              'structure')


    def _makeBoundsEvents(self, precise=True, eventtol=1e-6,
                          activatedbounds=None):
        events = []
        # pars + vars (pars used only by PyCont during continuation)
        alldoms = copy(self._pdomain)
        alldoms.update(self._xdomain)
        allnames = self.funcspec.vars + self.funcspec.pars
        if activatedbounds is None:
            activatedbounds = {}.fromkeys(allnames, (True,True))
        for xname, xdom in alldoms.iteritems():
            if xname not in allnames:
                # don't make bound constraints for non-state variables
                continue
            xdlo = xdom[0]
            xdhi = xdom[1]
            evname = xname+"_domlo"
            evargs = {'term': True,
                      'precise': precise,
                      'name': evname,
                      'eventtol': eventtol,
                      'eventdelay': eventtol*10,
                      'eventinterval': eventtol*20,
                      'xdomain': self._xdomain,
                      'pdomain': self._pdomain}
            try:
                evargs['active'] = activatedbounds[xname][0] and isfinite(xdlo)
            except KeyError:
                evargs['active'] = False
            evstr = xname + '-' + 'getbound("%s", 0)'%xname
            ev = Events.makeZeroCrossEvent(evstr, -1, evargs, [xname],
                                       targetlang=self.funcspec.targetlang,
                                       reuseterms=self.funcspec.reuseterms)
            events.append(ev)
            evname = xname+"_domhi"
            evargs = {'term': True,
                      'precise': precise,
                      'name': evname,
                      'eventtol': eventtol,
                      'eventdelay': eventtol*10,
                      'eventinterval': eventtol*20,
                      'xdomain': self._xdomain,
                      'pdomain': self._pdomain}
            try:
                evargs['active'] = activatedbounds[xname][1] and isfinite(xdhi)
            except KeyError:
                evargs['active'] = False
            evstr = xname + '-' + 'getbound("%s", 1)'%xname
            ev = Events.makeZeroCrossEvent(evstr, 1, evargs, [xname],
                                       targetlang=self.funcspec.targetlang,
                                       reuseterms=self.funcspec.reuseterms)
            events.append(ev)
        if events != []:
            self._addEvents(events)


    def set(self, **kw):
        """Set generic parameters."""
        # Note: globalt0 may be unused by many classes of Generator
        if len(kw) > 0:
            if 'globalt0' in kw:
                self.globalt0 = kw['globalt0']
            if 'checklevel' in kw:
                self.checklevel = kw['checklevel']
            if remain(kw.keys(), ['globalt0', 'checklevel']) != []:
                raise PyDSTool_KeyError, 'Invalid keywords passed'


    def setEventICs(self, ics, gt0=0):
        """Set initialconditions attribute of all generator's events, in
        case event uses auxiliary functions that access this information."""
        for ev in self.eventstruct.events.values():
            ev.initialconditions = ics
            ev.globalt0 = gt0


    # Auxiliary functions for user-defined code to call

    def _auxfn_globalindepvar(self, parsinps, t):
        return self.globalt0 + t


    def _auxfn_initcond(self, parsinps, varname):
        return self.initialconditions[varname]


    def _auxfn_heav(self, parsinps, x):
        if x>0:
            return 1
        else:
            return 0


    def _auxfn_if(self, parsinps, c, e1, e2):
        if c:
            return e1
        else:
            return e2


    def _auxfn_getindex(self, parsinps, varname):
        return self._var_namemap[varname]


    def clearWarnings(self):
        self.warnings = []
        if self.inputs:
            for iname in self.inputs:
                self.inputs[iname].clearWarnings()


    def showWarnings(self):
        if len(self.warnings)>0:
            print 'Warnings:'
            for (w, d) in self.warnings:
                dstr = ''
                for i in range(len(d)):
                    dentry = d[i]
                    dstr += self._warnfields[w][i] + ' = ' + str(dentry) + ", "
                dstr = dstr[:-2]  # drop trailing comma
                print 'Warning code %s:  %s\n Info:  %s' %(w, \
                                                self._warnmessages[w], dstr)


    def clearErrors(self):
        self.errors = []


    def showErrors(self):
        if len(self.errors)>0:
            print 'Errors:'
            for (e, d) in self.errors:
                dstr = ''
                for i in range(len(d)):
                    dentry = d[i]
                    dstr += self._errorfields[e][i] + ' = ' + str(dentry) + ", "
                dstr = dstr[:-2]  # drop trailing comma
                print 'Error code %s:  %s\n Info:\n  %s' %(e, \
                                                    self._errmessages[e], dstr)


    def _generate_ixmaps(self, gentypes=None):
        """Generate indices mapping.

        This creates a mapping from the names of variables,
        pars and inputs, to indices in the arrays used for
        refering to the internal (dynamic) call methods."""

        self.validateSpec()
        if gentypes is not None:
            if isinstance(gentypes, str):
                gentypes = [gentypes]
            for s in gentypes:
                assert s in ['variables', 'inputs', 'pars'], \
                       ('Incorrect type string for _generate_ixmaps')
        else:
            # default to all
            gentypes = ['variables', 'inputs', 'pars']
        # ixmap (list) : int -> str
        # namemap (dict) : str -> int
        if 'variables' in gentypes:
            self._var_ixmap = sortedDictKeys(self.variables,
                                             self.funcspec.vars)
            self._var_namemap = invertMap(self._var_ixmap)
        if 'pars' in gentypes:
            if self.pars:
                self._parameter_ixmap = sortedDictKeys(self.pars)
                self._parameter_namemap = invertMap(self._parameter_ixmap)
            else:
                self._parameter_ixmap = []
                self._parameter_namemap = {}
        if 'inputs' in gentypes:
            if self.inputs:
                self._inputs_ixmap = \
                    sortedDictKeys(self.inputs)
                self._inputs_namemap = invertMap(self._inputs_ixmap)
            else:
                self._inputs_ixmap = []
                self._inputs_namemap = {}


    def contains(self, interval, val, checklevel=2):
        # NB. val may be another interval
        if checklevel == 0:
            # level 0 -- no bounds checking at all
            # code should avoid calling this function with checklevel = 0
            # if possible, but this case is left here for completeness and
            # consistency
            return True
        elif checklevel == 2:
            # level 2 -- warn on uncertain and continue
            testresult = interval.contains(val)
            if testresult is contained:
                return True
            elif testresult is uncertain:
                self.warnings.append((W_UNCERTVAL, (val,interval)))
                return True
            else:
                return False
        elif checklevel == 1:
            # level 1 -- ignore uncertain cases (treat as contained)
            if interval.contains(val) is not notcontained:
                return True
            else:
                return False
        else:
            # level 3 -- exception will be raised for uncertain case
            if val in interval:
                return True
            else:
                return False


    # Methods for pickling protocol
    def __getstate__(self):
        d = copy(self.__dict__)
        for fname, finfo in self._funcreg.iteritems():
            try:
                del d[fname]
            except KeyError:
                pass
        return d


    def __setstate__(self, state):
        self.__dict__.update(state)
        if self._funcreg != {}:
            self.addMethods()


    def __del__(self):
        # delete object-specific class methods etc. before deleting
        # to avoid crowding namespace
        try:
            for fname, finfo in self._funcreg.iteritems():
                try:
                    delattr(eval(finfo[0]), fname)
                except AttributeError:
                    pass
                except NameError:
                    # not sure what happens here, but some other names
                    # may be deleted before all references to them have
                    # been deleted, but it's very non-fatal to ignore.
                    pass
            if hasattr(self, 'eventstruct'):
                if self.eventstruct is not None:
                    self.eventstruct.__del__()
            if self.indepvariable is not None:
                del self.indepvariable
            for v in self.variables.values():
                v.__del__()
            if hasattr(self, 'inputs'):
                for v in self.inputs.values():
                    v.__del__()
        except AttributeError:
            # self does not have _funcreg
            pass
        except NameError:
            # see above notes for NameError catch
            pass


    def __copy__(self):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


    def __deepcopy__(self, memo=None, _nil=[]):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


#--------------------------------------------------------------------------


class ctsGen(Generator):
    "Abstract class for continuously-parameterized trajectory generators."

    def validateSpec(self):
        # only check that domain is cts, range may be a finite subset of points
        assert isinputcts(self.indepvariable), ("self.indepvariable must be continuously-"
                                         "defined for this class")

    def __del__(self):
        Generator.__del__(self)



class discGen(Generator):
    "Abstract class for discretely-parameterized trajectory generators."

    def validateSpec(self):
        assert isdiscrete(self.indepvariable), ("self.indepvariable must be discretely-"
                                         "defined for this class")


    def __del__(self):
        Generator.__del__(self)


#--------------------------------------------------------------------------


class GenSpecInfoObj(object):
    # empty class struct for GenSpecHelper
    pass
    

class GenSpecHelper(object):
    """Generator specification helper - abstract class.

    Used to help ModelConstructor translate abstract model specifications
    into concrete specifications specific to individual Generators."""

    def __init__(self):
        self.gshDB = {}

    def add(self, genClass, symbolMapDict, lang, specType='RHSfuncSpec'):
        genName = className(genClass)
        if genName in self.gshDB:
            raise ValueError("Generator %s has already been declared"%genName)
        else:
            infoObj = GenSpecInfoObj()
            infoObj.genClass = genClass
            infoObj.symbolMap = symbolMapClass(symbolMapDict)
            infoObj.lang = lang
            infoObj.specType = specType
            if compareBaseClass(genClass, ctsGen):
                infoObj.domain = Continuous
            else:
                infoObj.domain = Discrete
            self.gshDB[genName] = infoObj

    def __call__(self, subject):
        try:
            if isinstance(subject, str):
                return self.gshDB[subject]
            else:
                return self.gshDB[className(subject)]
        except KeyError:
            raise KeyError("Generator %s was not found in database"%str(subject))

    def __contains__(self, subject):
        return subject in self.gshDB or className(subject) in self.gshDB



global theGenSpecHelper
theGenSpecHelper = GenSpecHelper()
