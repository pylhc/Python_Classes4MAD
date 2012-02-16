# ODEsystem base class

from allimports import *
from baseclasses import ctsGen, theGenSpecHelper
from PyDSTool.utils import *
from PyDSTool.common import *
from PyDSTool.Symbolic import ensureStrArgDict
from PyDSTool.parseUtils import NAMESEP
from PyDSTool.Variable import Variable, iscontinuous
from PyDSTool.Trajectory import Trajectory
from PyDSTool.Points import Pointset

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, transpose, shape
import math, random
import os, types
from copy import copy, deepcopy
try:
    # use pscyo JIT byte-compiler optimization, if available
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

# -----------------------------------------------------------------------------
# Version History
#
# 11 September 2006
#   Fixed creation of external input Variable object to properly use Pointsets.
#
# 30 December 2005
#   Vfield renamed to Rhs, because for Radau (and possibly other integrators
#     in the future) that use implicitly defined vector fields, this function
#     will only return the value of the RHS of the definition.
#
# 15 November 2005
#   Added haveMass() method.
#
# 10 September 2005
#   Added haveJacobian_pars() method.
#
# 28 August 2005
#   Added support for reusable terms in function specifications, provided by
#     user at initialization. 2 passes are possible (2nd order substitutions)
#   Added support for arbitrary code to be inserted at beginning of vector
#     field function definition. Code must be in appropriate target language.
#   Added C compiler defaults.
#
# 17 July 2005
#   Added support for singleton xdomains
#
# 14 July 2005
#   Added initial conditions checks in validateICs().
#
# -----------------------------------------------------------------------------

class ODEsystem(ctsGen):
    """Abstract class for ODE system solvers."""

    def __init__(self, kw):
        ctsGen.__init__(self, kw)
        self._errmessages[E_COMPUTFAIL] = 'Integration failed'
        self.needKeys.extend(['varspecs'])
        self.optionalKeys.extend(['tdomain', 'xdomain', 'inputs', 'tdata',
                          'ics', 'events', 'compiler', 'enforcebounds',
                          'activatedbounds', 'checklevel', 'algparams',
                          'auxvars', 'vars', 'pars', 'fnspecs', 'pdomain',
                          'reuseterms', 'vfcodeinsert_start', 'vfcodeinsert_end'])
        fs_args = {}
        if 'varspecs' in kw:
            self.foundKeys += 1
            varspecs = ensureStrArgDict(kw['varspecs'])
        else:
            raise PyDSTool_KeyError("Keyword 'varspecs' missing from "
                                    "argument")
        if 'tdomain' in kw:
            self._tdomain = kw['tdomain']
            if self._tdomain[0] >= self._tdomain[1]:
                print "Time domain specified: [%s, %s]"%(self._tdomain[0], self._tdomain[1])
                raise PyDSTool_ValueError('tdomain values must be in order of increasing size')
            self.foundKeys += 1
        else:
            self._tdomain = [-Inf, Inf]
        if 'tdata' in kw:
            self._tdata = kw['tdata']
            if self._tdata[0] >= self._tdata[1]:
                raise PyDSTool_ValueError('tdata values must be in order of increasing size')
            # _tdata is made into list to be consistent with
            # other uses of it in other Generators...
            if self._tdomain[0] > self._tdata[0]:
                print 'tdata cannot be specified below smallest '\
                      'value in tdomain\n (possibly due to uncertain'\
                      'bounding). It has been automatically adjusted\n'\
                      ' from\n ', self._tdata[0], 'to', self._tdomain[0],\
                      '(difference of', self._tdomain[0]-self._tdata[0], ')'
                if self._modeltag:
                      print 'Try reducing step size in model.'
                self._tdata[0] = self._tdomain[0]
            if self._tdomain[1] < self._tdata[1]:
                print 'tdata cannot be specified above largest '\
                      'value in tdomain\n (possibly due to uncertain '\
                      'bounding). It has been automatically adjusted\n'\
                      ' from\n ', self._tdata[1], 'to', self._tdomain[1],\
                      '(difference of', self._tdata[1]-self._tdomain[1], ')'
                if self._modeltag:
                      print 'Try reducing step size in model.'
                self._tdata[1] = self._tdomain[1]
            self.foundKeys += 1
        else:
            self._tdata = self._tdomain  # default needed
        if 'inputs' in kw:
            inputs = copy(kw['inputs'])
            self.inputs = {}
            if isinstance(inputs, Trajectory):
                # extract the variables
                self.inputs = inputs.variables
            elif isinstance(inputs, Variable):
                self.inputs = {inputs.name: inputs}
            elif isinstance(inputs, Pointset):
                # turn into Variables with linear interpoolation between
                # independent variable values
                for n in inputs.coordnames:
                    x_array = inputs[n]
                    self.inputs[n] = Variable(interp1d(inputs.indepvararray,
                                                       x_array), 't',
                                         Interval(n, 'float', extent(x_array)),
                                         name=n)
            elif isinstance(inputs, dict):
                self.inputs = inputs
                # ensure values are Variables or Pointsets
                for k, v in self.inputs.iteritems():
                    if not isinstance(v, Variable):
                        try:
                            self.inputs[k]=Variable(v)
                        except:
                            raise TypeError, "Invalid specification of inputs"
            else:
                raise TypeError, "Invalid specification of inputs"
            for v in self.inputs.values():
                if not iscontinuous(v):
                    raise ValueError, \
                          ("External inputs for ODE system must be "
                           "continuously defined")
            fs_args['inputs'] = self.inputs.keys()
            self._register(self.inputs)
            self.foundKeys += 1
            self._extInputsChanged = True
        else:
            self.inputs = {}
            self._extInputsChanged = False
        if 'ics' in kw:
            self._xdatadict = {}
            for k, v in dict(kw['ics']).iteritems():
                self._xdatadict[str(k)] = v
            self.initialconditions = self._xdatadict
            for name in remain(varspecs.keys(),
                               self._xdatadict.keys()):
                self.initialconditions[name] = NaN
            self.foundKeys += 1
        else:
            self._xdatadict = {}
            for name in varspecs:
                self.initialconditions[name] = NaN
        if 'auxvars' in kw:
            assert 'vars' not in kw, ("Cannot use both 'auxvars' and 'vars' "
                                      "keywords")
            if isinstance(kw['auxvars'], list):
                auxvars = [str(v) for v in kw['auxvars']]
            else:
                auxvars = [str(kw['auxvars'])]
            vars = remain(varspecs.keys(), auxvars)
            self.foundKeys += 1
        elif 'vars' in kw:
            assert 'auxvars' not in kw, \
                   "Cannot use both 'auxvars' and 'vars' keywords"
            if isinstance(kw['vars'], list):
                vars = [str(v) for v in kw['vars']]
            else:
                vars = [str(kw['vars'])]
            auxvars = remain(varspecs.keys(), vars)
            self.foundKeys += 1
        else:
            auxvars = []
            vars = varspecs.keys()
        if auxvars != []:
            fs_args['auxvars'] = auxvars
        if 'xdomain' in kw:
            self._xdomain = {}
            for k, v in dict(kw['xdomain']).iteritems():
                name = str(k)
                if type(v) in [list, Array, NArray]:
                    assert len(v) == 2, \
                           "Invalid size of domain specification for "+name
                    if v[0] >= v[1]:
                        raise PyDSTool_ValueError('xdomain values must be in'
                                                  'order of increasing size')
                    else:
                        self._xdomain[name] = copy(v)
                elif type(v) in [float, int]:
                    self._xdomain[name] = [v, v]
                else:
                    raise PyDSTool_TypeError('Invalid type for xdomain spec'
                                             ' '+name)
            for name in remain(varspecs.keys(), self._xdomain.keys()):
                self._xdomain[name] = [-Inf, Inf]
            self.foundKeys += 1
        else:
            self._xdomain = {}
            for name in varspecs:
                self._xdomain[name] = [-Inf, Inf]
        if 'reuseterms' in kw:
            fs_args['reuseterms'] = kw['reuseterms']
            self.foundKeys += 1
        fs_args.update({'vars': vars,
                       'varspecs': varspecs,
                       'name': self.name
                       })
        if 'algparams' in kw:
            self._algparams = copy(kw['algparams'])
            self.foundKeys += 1
        else:
            self._algparams = {}
        if 'pars' in kw:
            self.pars = {}
            if type(kw['pars']) == list:
                # may be a list of symbolic definitions
                for p in kw['pars']:
                    try:
                        self.pars[p.name] = p.tonumeric()
                    except (AttributeError, TypeError):
                        raise TypeError, "Invalid parameter symbolic definition"
            else:
                for k, v in dict(kw['pars']).iteritems():
                    self.pars[str(k)] = v
            fs_args['pars'] = self.pars.keys()
            self._register(self.pars)
            self.foundKeys += 1
        if 'pdomain' in kw:
            if self.pars:
                self._pdomain = {}
                for k, v in dict(kw['pdomain']).iteritems():
                    assert len(v) == 2, \
                               "Invalid size of domain specification for "+k
                    self._pdomain[str(k)] = v
                for name in self._pdomain:
                    if self._pdomain[name][0] >= self._pdomain[name][1]:
                        raise PyDSTool_ValueError('pdomain values must be in order of increasing size')
                for name in remain(self.pars.keys(), self._pdomain.keys()):
                    self._pdomain[name] = [-Inf, Inf]
                self.foundKeys += 1
            else:
                raise ValueError('Cannot specify pdomain because no pars declared')
        else:
            if self.pars:
                self._pdomain = {}
                for pname in self.pars:
                    self._pdomain[pname] = [-Inf, Inf]
        if self.pars:
            self.parameterDomains = {}
            for pname in self._pdomain:
                self.parameterDomains[pname] = Interval(pname, 'float',
                                                        self._pdomain[pname],
                                                        self._abseps)
                try:
                    cval = self.parameterDomains[pname].contains(self.pars[pname])
                except KeyError:
                    raise ValueError("Parameter %s is missing a value"%pname)
                if self.checklevel < 3:
                    if cval is not notcontained:
                        if cval is uncertain and self.checklevel == 2:
                            print 'Warning: Parameter value at bound'
                    else:
                        raise PyDSTool_ValueError, 'Parameter value out of bounds'
                else:
                    if cval is uncertain:
                        raise PyDSTool_UncertainValueError, 'Parameter value at bound'
                    elif cval is notcontained:
                        raise PyDSTool_ValueError, 'Parameter value out of bounds'  
        if 'fnspecs' in kw:
            fs_args['fnspecs'] = ensureStrArgDict(kw['fnspecs'])
            self.foundKeys += 1
        fs_args['targetlang'] = theGenSpecHelper(self).lang
        if 'compiler' in kw:
            if fs_args['targetlang'] == 'python':
                print "Warning: redundant option 'compiler' for python target"
            self._compiler = kw['compiler']
            self.foundKeys += 1
        else:
            osname = os.name
            # os-specific defaults for C compiler
            if osname == 'nt':
                self._compiler = 'mingw32'
            elif osname == 'mac':
                self._compiler = 'mwerks'
            elif osname == 'posix' or osname == 'unix':
                self._compiler = 'unix'
            elif osname == 'os2emx':
                self._compiler = 'emx'
            else:
                self._compiler = ''
        if 'vfcodeinsert_start' in kw:
            fs_args['codeinsert_start'] = kw['vfcodeinsert_start']
            self.foundKeys += 1
        if 'vfcodeinsert_end' in kw:
            fs_args['codeinsert_end'] = kw['vfcodeinsert_end']
            self.foundKeys += 1
        self.funcspec = RHSfuncSpec(fs_args)
        # Holder and interface for events
        self.eventstruct = EventStruct()
        if 'enforcebounds' in kw:
            if 'activatedbounds' in kw:
                ab = kw['activatedbounds']
                self.foundKeys += 1
            else:
                ab = None
            if kw['enforcebounds']:
                self._makeBoundsEvents(precise=True, activatedbounds=ab)
            self.foundKeys += 1
        if 'events' in kw:
            self._addEvents(kw['events'])
            self.foundKeys += 1
        self.checkArgs(kw)
        self.numpars = len(self.pars)
        tindepdomain = Interval('t_domain', 'float', self._tdomain,
                                self._abseps)
        tdepdomain = Interval('t', 'float', self._tdata, self._abseps)
        self.indepvariable = Variable(listid, tindepdomain, tdepdomain, 't')
        self._register(self.indepvariable)
        self.dimension = len(self.funcspec.vars)
        for xname in self.funcspec.vars + self.funcspec.auxvars:
            # Add a temporary dependent variable domain, for validation testing
            # during integration
            self.variables[xname] = Variable(indepdomain=tdepdomain,
                                             depdomain=Interval(xname, 'float',
                                                          self._xdomain[xname],
                                                          self._abseps))
        self._register(self.variables)
        self._generate_ixmaps()
        # Introduce any python-specified code to the local namespace
        if fs_args['targetlang'] == 'python':
            self.addMethods()
        # all registration completed
        self.validateSpec()


    def prepDirection(self, dirn):
        """Common pre-integration tasks go here"""
        if dirn in ['f', 'forward', 1]:
            self._dircode = 1
            continue_integ = False
        elif dirn in ['b', 'backward', -1]:
            self._dircode = -1
            continue_integ = False
        elif dirn in ['c', 'continue', 0]:
            if self.defined and self._solver is not None:
                continue_integ = True
            else:
                # just treat as initial forwards integration
                continue_integ = False
                self._dircode = 1
        else:
            raise ValueError, 'Invalid direction code in argument'
        return continue_integ



    def addMethods(self, usePsyco=False):
        """Add Python-specific functions to this object's methods,
        accelerating them with psyco, if it is available."""

        # Add the auxiliary function specs to this Generator's namespace
        for auxfnname in self.funcspec.auxfns:
            fninfo = self.funcspec.auxfns[auxfnname]
            if not hasattr(self, fninfo[1]):
                # user-defined auxiliary functions
                # (built-ins are provided explicitly)
                try:
                    exec fninfo[0]
                except:
                    print 'Error in supplied auxiliary function code'
#                print "\n", fninfo[0], "\n"
                self._funcreg[fninfo[1]] = ('self', fninfo[0])
                #setattr(self.__class__, fninfo[1], eval(fninfo[1]))
                setattr(self, fninfo[1], types.MethodType(locals()[fninfo[1]],
                                                           self,
                                                           self.__class__))
            # bind all the auxfns here
            if HAVE_PSYCO and usePsyco:
                psyco.bind(getattr(self, fninfo[1]))
        # Add the spec function to this Generator's namespace
        fninfo = self.funcspec.spec
        try:
            exec fninfo[0]
        except:
            print 'Error in supplied functional specification code'
            raise
#        print "\n", fninfo[0], "\n"
        self._funcreg[fninfo[1]] = ('self', fninfo[0])
        setattr(self, fninfo[1], types.MethodType(locals()[fninfo[1]],
                                                           self,
                                                           self.__class__))
        if HAVE_PSYCO and usePsyco:
            psyco.bind(getattr(self, fninfo[1]))
        # Add the auxiliary spec function (if present) to this
        # Generator's namespace
        if self.funcspec.auxspec != '':
            fninfo = self.funcspec.auxspec
            try:
                exec fninfo[0]
            except:
                print 'Error in supplied auxiliary variable code'
                raise
#            print "\n", fninfo[0], "\n"
            self._funcreg[fninfo[1]] = ('self', fninfo[0])
            setattr(self, fninfo[1], types.MethodType(locals()[fninfo[1]],
                                                           self,
                                                           self.__class__))
            if HAVE_PSYCO and usePsyco:
                psyco.bind(getattr(self, fninfo[1]))


    def haveJacobian(self):
        """Report whether ODE system has an explicit user-specified Jacobian
        associated with it."""
        return 'Jacobian' in self.funcspec.auxfns


    def haveJacobian_pars(self):
        """Report whether ODE system has an explicit user-specified Jacobian
        with respect to pars associated with it."""
        return 'Jacobian_pars' in self.funcspec.auxfns


    def haveMass(self):
        """Report whether ODE system has an explicit user-specified mass matrix
        associated with it."""
        return 'massMatrix' in self.funcspec.auxfns


    def checkInitialConditions(self, checkauxvars=False):
        for xname, val in self.initialconditions.iteritems():
            if xname not in self.funcspec.vars and not checkauxvars:
                # auxvars do not need initial conditions unless
                # explicitly requested (e.g. for user call to RHS
                # function of C code vector field before trajectory
                # has been computed)
                continue
            try:
                if not isfinite(val):
                    raise ValueError, ("Initial condition for "+xname+" has been "
                                    "incorrectly initialized")
            except TypeError:
                print "Found: ", val
                print "of type: ", type(val)
                raise TypeError("Invalid type for %s`s initial"%xname \
                                + "condition value")
            if not self.contains(self.variables[xname].depdomain,
                                 val, self.checklevel):
                print "Bounds: ", self.variables[xname].depdomain.get()
                print "Variable value: ", val
                raise ValueError, ("Initial condition for "+xname+" has been "
                                   "set outside of prescribed bounds")


    def set(self, **kw):
        """Set ODE system parameters"""

        validKeys = ['globalt0', 'xdomain', 'tdata', 'tdomain', 'checklevel',
                     'ics', 'pars', 'algparams', 'inputs', 'pdomain']
        assert remain(kw.keys(), validKeys) == [], "Invalid keys in argument"
        if 'globalt0' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, globalt0=kw['globalt0'])
            # set this True so that can adjust the input time arrays
            # to new globalt0
            self._extInputsChanged = True
        if 'checklevel' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, checklevel=kw['checklevel'])
        # optional keys for this call are ['pars', 'tdomain', 'ics',
        #   'algparams', 'tdata', 'xdomain', 'pdomain', 'inputs']
        if 'ics' in kw:
            for k_temp, v in kw['ics'].iteritems():
                k = k_temp.replace(NAMESEP, "_")
                if k in self.funcspec.vars+self.funcspec.auxvars:
                    self._xdatadict[k] = v
                else:
                    raise ValueError, 'Illegal variable name %s'%k
            self.initialconditions.update(self._xdatadict)
        if 'tdata' in kw:
            self._tdata = kw['tdata']
        if 'tdomain' in kw:
            self._tdomain = kw['tdomain']
            self.indepvariable.indepdomain.set(self._tdomain)
        if self._tdomain[0] > self._tdata[0]:
            print 'tdata cannot be specified below smallest '\
                  'value in tdomain\n (possibly due to uncertain bounding).'\
                  ' It has been automatically adjusted from\n ', \
                  self._tdata[0], 'to', self._tdomain[0], '(difference of', \
                  self._tdomain[0]-self._tdata[0], ')'
            if self._modeltag:
                  print 'Try reducing step size in model.'
            self._tdata[0] = self._tdomain[0]
        if self._tdomain[1] < self._tdata[1]:
            print 'tdata cannot be specified above largest '\
                  'value in tdomain\n (possibly due to uncertain bounding).'\
                  ' It has been automatically adjusted from\n ', \
                  self._tdomain[1], 'to', \
                  self._tdomain[1], '(difference of', \
                  self._tdata[1]-self._tdomain[1], ')'
            if self._modeltag:
                  print 'Try reducing step size in model.'
            self._tdata[1] = self._tdomain[1]
        self.indepvariable.depdomain.set(self._tdata)
        if 'xdomain' in kw:
            for k_temp, v in kw['xdomain'].iteritems():
                k = k_temp.replace(NAMESEP, "_")   # use internal names
                if k in self.funcspec.vars+self.funcspec.auxvars:
                    if type(v) in [list, Array, NArray]:
                        assert len(v) == 2, \
                               "Invalid size of domain specification for "+k
                        if v[0] >= v[1]:
                            raise PyDSTool_ValueError('xdomain values must be'
                                                      'in order of increasing '
                                                      'size')
                        else:
                            self._xdomain[k] = copy(v)
                    elif type(v) in [float, int]:
                        self._xdomain[k] = [v, v]
                    else:
                        raise PyDSTool_TypeError('Invalid type for xdomain spec'
                                                 ' '+k)
                else:
                    raise ValueError, 'Illegal variable name'
                try:
                    self.variables[k].depdomain.set(v)
                except TypeError:
                    raise TypeError, ('xdomain must be a dictionary of variable'
                                      ' names -> valid interval 2-tuples or '
                                      'singletons')
                for ev in self.eventstruct.events.values():
                    ev._xdomain[k] = self._xdomain[k]
        if 'pdomain' in kw:
            for k_temp, v in kw['pdomain'].iteritems():
                k = k_temp.replace(NAMESEP, "_")   # use internal names
                if k in self.funcspec.pars:
                    if type(v) in [list, Array, NArray]:
                        assert len(v) == 2, \
                               "Invalid size of domain specification for "+k
                        if v[0] >= v[1]:
                            raise PyDSTool_ValueError('pdomain values must be'
                                                      'in order of increasing '
                                                      'size')
                        else:
                            self._pdomain[k] = copy(v)
                    elif type(v) in [float, int]:
                        self._pdomain[k] = [v, v]
                    else:
                        raise PyDSTool_TypeError('Invalid type for pdomain spec'
                                                 ' '+k)
                else:
                    raise ValueError, 'Illegal parameter name'
                try:
                    self.parameterDomains[k].depdomain.set(v)
                except TypeError:
                    raise TypeError, ('xdomain must be a dictionary of parameter'
                                      ' names -> valid interval 2-tuples or '
                                      'singletons')
                for ev in self.eventstruct.events.values():
                    ev._pdomain[k] = self._pdomain[k]
        if 'pars' in kw:
            assert self.numpars > 0, ('No pars were declared for this '
                                      'model')
            for k_temp, v in kw['pars'].iteritems():
                k = k_temp.replace(NAMESEP, "_")
                if k in self.pars:
                    cval = self.parameterDomains[k].contains(v)
                    if self.checklevel < 3:
                        if cval is not notcontained:
                            self.pars[k] = v
                            if cval is uncertain and self.checklevel == 2:
                                print 'Warning: Parameter value at bound'
                        else:
                            raise PyDSTool_ValueError, 'Parameter value out of bounds'
                    else:
                        if cval is contained:
                            self.pars[k] = v
                        elif cval is uncertain:
                            raise PyDSTool_UncertainValueError, 'Parameter value at bound'
                        else:
                            raise PyDSTool_ValueError, 'Parameter value out of bounds'
                else:
                    raise PyDSTool_ValueError, 'Illegal parameter name'
        if 'algparams' in kw:
            for k, v in kw['algparams'].iteritems():
                self._algparams[k] = v
        if 'inputs' in kw:
            assert self.inputs, ('Cannot provide inputs after '
                                             'initialization without them')
            inputs = copy(kw['inputs'])
            _inputs = {}
            if isinstance(inputs, Trajectory):
                # extract the variables
                _inputs = inputs.variables
            elif isinstance(inputs, Variable):
                _inputs = {inputs.name: inputs}
            elif isinstance(inputs, Pointset):
                # turn into Variables
                _inputs = dict(zip(inputs.coordnames,
                            [Variable(inputs[n],name=n) for n in inputs.coordnames]))
            elif isinstance(inputs, dict):
                _inputs = inputs
                # ensure values are Variables or Pointsets
                for k, v in _inputs.iteritems():
                    if not isinstance(v, Variable):
                        try:
                            _inputs[k]=Variable(v)
                        except:
                            raise TypeError, "Invalid specification of inputs"
            else:
                raise TypeError, "Invalid specification of inputs"
            for v in _inputs.values():
                if not iscontinuous(v):
                    raise ValueError, "External inputs for ODE system must be continously defined"
            if _inputs:
                for i in _inputs:
                    assert i in self.inputs, 'Incorrect input name provided'
                    self.inputs[i] = _inputs[i]
                # re-calc inputs ixmaps
                self._generate_ixmaps('inputs')
            self._extInputsChanged = True                    


    def compute(self, dirn):
        """This is an abstract class."""

        # Don't call this from a concrete sub-class
        if self.__class__ is ODEsystem:
            raise NotImplementedError, ('ODEsystem is an abstract class. '
                              'Use a concrete sub-class of ODEsystem.')
        else:
            raise NotImplementedError, ("This Generator does not support "
                              "calls to compute() or computeTraj()")


    # for backwards compatibility
    setPars = set
    computeTraj = compute


    def validateSpec(self):
        ctsGen.validateSpec(self)
        if remain(self.initialconditions.keys(),
                      self.variables.keys()) != []:
            print "IC names defined:", self.initialconditions.keys()
            print "Varnames defined:", self.variables.keys()
            raise ValueError, \
                      "Mismatch between initial condition and variable names"


    def validateICs(self):
        assert hasattr(self, 'initialconditions')
        for v in self.funcspec.vars:
            try:
                assert self.initialconditions[v] is not NaN, ('Must specify'
                                    ' initial conditions for all variables')
            except KeyError:
                raise KeyError, ("Variable name "+v+" not found in initial"
                                 " conditions")


    def Rhs(self, t, xdict, pdict, idict):
        # Don't call this from a concrete sub-class
        if self.__class__ is ODEsystem:
            raise NotImplementedError, ('ODEsystem is an abstract class. '
                              'Use a concrete sub-class of ODEsystem.')
        else:
            raise NotImplementedError, ("This Generator does not support "
                              "calls to Rhs()")


    def Jacobian(self, t, xdict, pdict, idict):
        # Don't call this from a concrete sub-class
        if self.__class__ is ODEsystem:
            raise NotImplementedError, ('ODEsystem is an abstract class. '
                              'Use a concrete sub-class of ODEsystem.')
        else:
            raise NotImplementedError, ("This Generator does not support "
                              "calls to Jacobian()")


    def JacobianP(self, t, xdict, pdict, idict):
        # Don't call this from a concrete sub-class
        if self.__class__ is ODEsystem:
            raise NotImplementedError, ('ODEsystem is an abstract class. '
                              'Use a concrete sub-class of ODEsystem.')
        else:
            raise NotImplementedError, ("This Generator does not support "
                              "calls to JacobianP()")


    def AuxFunc(self, t, xdict, pdict, idict):
        # Don't call this from a concrete sub-class
        if self.__class__ is ODEsystem:
            raise NotImplementedError, ('ODEsystem is an abstract class. '
                              'Use a concrete sub-class of ODEsystem.')
        else:
            raise NotImplementedError, ("This Generator does not support "
                              "calls to AuxFunc()")


    # Methods for pickling protocol
    def __getstate__(self):
        d = copy(self.__dict__)
        for fname, finfo in self._funcreg.iteritems():
            try:
                del d[fname]
            except KeyError:
                pass
        # can't pickle module objects, so ensure that _solver deleted:
        # Dopri and Radau don't have _solver in funcreg
        try:
            del d['_solver']
        except KeyError:
            pass
        return d


    def __setstate__(self, state):
        self.__dict__.update(state)
        self._solver = None
        if self._funcreg != {}:
            self.addMethods()


    def __del__(self):
        ctsGen.__del__(self)


