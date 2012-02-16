# MapSystem
# For the MapSystem class, maps have an explicit "time" step, generating
#  values for x at a specifc "time". Purely abstract time (i.e., iteration
#  steps) is represented using integers. (We could make LookupTable a
#  0-param sequence of maps with explicit time range?)
from __future__ import division

from allimports import *
from baseclasses import Generator, discGen, theGenSpecHelper
from PyDSTool.utils import *
from PyDSTool.common import *
from PyDSTool.Symbolic import ensureStrArgDict
from PyDSTool.parseUtils import NAMESEP
from PyDSTool.Variable import Variable
from PyDSTool.Trajectory import Trajectory
from PyDSTool.Points import Pointset

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, transpose, shape
import math, random, types
from copy import copy, deepcopy
try:
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

# -----------------------------------------------------------------------------
# Version History
#
# 21 June 2006
#   Added reuseterms support
#
# 14 Mar 2006
#   Added support for verbose code insertion using vfcodeinsert_start and _end.
#   Added support for numerical definition of the map, by supplying 'system'.
#
# 22 Jan 2006
#   Added trajevents structure and cleaned up terminal event processing.
#
# 15 Oct 2005
#   Added checkInitialConditions method, identical to that for ODE solvers.
#
# 15 Sept 2005
#   Added support for Psyco JIT byte-compiler optimization, if available.
#     (x86 platforms only).
#
# 17 July 2005
#   Added support for singleton xdomains.
#
# -----------------------------------------------------------------------------

class MapSystem(discGen):
    """Discrete dynamical systems, as maps (difference equations).

    """

    def __init__(self, kw):
        discGen.__init__(self, kw)
        self._errmessages[E_COMPUTFAIL] = 'Computation of trajectory failed'
        self.needKeys.extend(['varspecs'])
        self.optionalKeys.extend(['tdomain', 'xdomain', 'inputs', 'tdata',
                          'ics', 'events', 'system', 'ignorespecial', #'compiler',
                          'auxvars', 'vars', 'fnspecs', 'ttype', 'xtype',
                          'tstep', 'checklevel', 'pars', 'pdomain',
                          'vfcodeinsert_start', 'vfcodeinsert_end',
                          'enforcebounds', 'activatedbounds', 'reuseterms'])
        fs_args = {}
        if 'varspecs' in kw:
            self.foundKeys += 1
            varspecs = ensureStrArgDict(kw['varspecs'])
        else:
            raise PyDSTool_KeyError("Keyword 'varspecs' missing from "
                                    "argument")
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
        if 'vfcodeinsert_start' in kw:
            fs_args['codeinsert_start'] = kw['vfcodeinsert_start']
            self.foundKeys += 1
        if 'vfcodeinsert_end' in kw:
            fs_args['codeinsert_end'] = kw['vfcodeinsert_end']
            self.foundKeys += 1
        if 'system' in kw:
            self._solver = kw['system']
            try:
                fs_args['ignorespecial'] = [self._solver.name]
            except:
                raise TypeError, "Invalid solver system provided"
            self.foundKeys += 1
            if self.pars:
                # automatically pass par values on to embedded system
                # when Rhs called
                parlist = self.pars.keys()
                parstr = "".join(["'%s': %s, "%(parname,parname) \
                                  for parname in parlist])
            if 'codeinsert_start' in fs_args:
                fs_args['codeinsert_start'] = \
                    '    %s.set(pars={%s})\n'%(self._solver.name, parstr) \
                    + fs_args['codeinsert_start']
            else:
                fs_args['codeinsert_start'] = \
                    '    %s.set(pars={%s})\n'%(self._solver.name, parstr)
        else:
            fs_args['ignorespecial'] = []
            self._solver = None
        if 'ignorespecial' in kw:
            fs_args['ignorespecial'].extend(kw['ignorespecial'])
            self.foundKeys += 1
        if 'tdomain' in kw:
            self._tdomain = kw['tdomain']
            if self._tdomain[0] >= self._tdomain[1]:
                raise PyDSTool_ValueError('tdomain values must be in order of increasing size')
            self.foundKeys += 1
        else:
            self._tdomain = [-Inf, Inf]
        if 'ttype' in kw:
            indepvartype = kw['ttype']
            assert indepvartype in _pynametype.keys(), 'Invalid ttype'
            self.foundKeys += 1
        else:
            indepvartype = 'float'
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
                self._tdata[0] = self._tdomain[0]
            if self._tdomain[1] < self._tdata[1]:
                print 'tdata cannot be specified above largest '\
                      'value in tdomain\n (possibly due to uncertain '\
                      'bounding). It has been automatically adjusted\n'\
                      ' from\n ', self._tdata[1], 'to', self._tdomain[1],\
                      '(difference of', self._tdata[1]-self._tdomain[1], ')'
                self._tdata[1] = self._tdomain[1]
            self.foundKeys += 1
        else:
            self._tdata = self._tdomain  # default needed
        if 'tstep' in kw:
            self.tstep = kw['tstep']
            if self.tstep > self._tdata[1]-self._tdata[0]:
                raise PyDSTool_ValueError('tstep too large')
            if indepvartype == 'int' and round(self.tstep) != self.tstep:
                raise PyDSTool_ValueError('tstep must be an integer for integer ttype')
            self.foundKeys += 1
        else:
            if indepvartype == 'int':
                # default to 1 for integer types
                self.tstep = 1
            else:
                # no reasonable default - so raise error
                raise PyDSTool_KeyError('tstep key needed for float ttype')
        if 'inputs' in kw:
            inputs = copy(kw['inputs'])
            self.inputs = {}
            if isinstance(inputs, Trajectory):
                # extract the variables
                self.inputs = inputs.variables
            elif isinstance(inputs, Variable):
                self.inputs = {inputs.name: inputs}
            elif isinstance(inputs, Pointset):
                # turn into Variables but only define at Pointset's
                # own independent variable values
                self.inputs = dict(zip(inputs.coordnames,
                            [Variable(inputs[[n]],name=n) for n in inputs.coordnames]))
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
            # don't check that variables have discrete domains in case wish
            # to use cts-valued ones sampled at discrete times
            fs_args['inputs'] = self.inputs.keys()
            self._register(self.inputs)
            self.foundKeys += 1
        else:
            self.inputs = {}
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
        if 'xtype' in kw:
            depvartype = kw['xtype']
            assert depvartype in _pynametype.keys(), 'Invalid xtype'
            self.foundKeys += 1
        else:
            depvartype = 'float'
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
        if 'pdomain' in kw:
            if self.pars:
                self._pdomain = {}
                for k, v in dict(kw['pdomain']).iteritems():
                    assert len(v) == 2, "Invalid interval for parameter %s"%k
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
##        fs_args['targetlang'] = 'python'
##        if 'compiler' in kw:
##            if fs_args['targetlang'] == 'python':
##                print "Warning: redundant option 'compiler' for python target"
##            self._compiler = kw['compiler']
##            self.foundKeys += 1
##        else:
##            self._compiler = ''
##        if fs_args['targetlang'] != 'python' and inputs:
##            raise PyDSTool_ValueError('Cannot use external inputs with non-Python target language')
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
                self._makeBoundsEvents(precise=False, activatedbounds=ab)
            self.foundKeys += 1
        if 'events' in kw:
            self._addEvents(kw['events'])
            self.foundKeys += 1
        self.checkArgs(kw)
        self.numpars = len(self.pars)
        tindepdomain = Interval('t_domain', indepvartype, self._tdomain,
                                self._abseps)
        tdepdomain = Interval('t', indepvartype, self._tdata, self._abseps)
        self.indepvariable = Variable(listid, tindepdomain, tdepdomain, 't')
        self._register(self.indepvariable)
        self.dimension = len(self.funcspec.vars)
        for xname in self.funcspec.vars + self.funcspec.auxvars:
            # Add a temporary dependent variable domain, for validation testing
            # during integration
            self.variables[xname] = Variable(indepdomain=tdepdomain,
                                         depdomain=Interval(xname, depvartype,
                                                          self._xdomain[xname],
                                                          self._abseps))
        self._register(self.variables)
        self._generate_ixmaps()
        # Introduce any python-specified code to the local namespace
##        if fs_args['targetlang'] == 'python':
        self.addMethods()
        # all registration completed
        self.validateSpec()


    def addMethods(self):
        # Add the auxiliary function specs to this Generator's namespace
        for auxfnname in self.funcspec.auxfns:
            fninfo = self.funcspec.auxfns[auxfnname]
            if not hasattr(self, fninfo[1]):
                # user-defined auxiliary functions
                # (built-ins are provided explicitly)
                if self._solver:
                    fnstr = fninfo[0].replace(self._solver.name, 'ds._solver')
##                    self._funcreg[self._solver.name] = self._solver
                else:
                    fnstr = fninfo[0]
                try:
                    exec fnstr
                except:
                    print 'Error in supplied auxiliary function code'
                self._funcreg[fninfo[1]] = ('self', fnstr)
                setattr(self, fninfo[1], types.MethodType(locals()[fninfo[1]],
                                                           self,
                                                           self.__class__))
            # bind all the auxfns here
            if HAVE_PSYCO:
                psyco.bind(getattr(self, fninfo[1]))
        # Add the spec function to this Generator's namespace
        fninfo = self.funcspec.spec
        if self._solver:
            fnstr = fninfo[0].replace(self._solver.name, 'ds._solver')
        else:
            fnstr = fninfo[0]
        try:
            exec fnstr
        except:
            print 'Error in supplied functional specification code'
            raise
        self._funcreg[fninfo[1]] = ('self', fnstr)
        setattr(self, fninfo[1], types.MethodType(locals()[fninfo[1]],
                                                           self,
                                                           self.__class__))
        if HAVE_PSYCO and not self._solver:
            psyco.bind(getattr(self, fninfo[1]))
        # Add the auxiliary spec function (if present) to this
        # Generator's namespace
        if self.funcspec.auxspec != '':
            fninfo = self.funcspec.auxspec
            if self._solver:
                fnstr = fninfo[0].replace(self._solver.name, 'ds._solver')
            else:
                fnstr = fninfo[0]
            try:
                exec fnstr
            except:
                print 'Error in supplied auxiliary variable code'
                raise
            self._funcreg[fninfo[1]] = ('self', fnstr)
            setattr(self, fninfo[1], types.MethodType(locals()[fninfo[1]],
                                                           self,
                                                           self.__class__))
            if HAVE_PSYCO and not self._solver:
                psyco.bind(getattr(self, fninfo[1]))


    def haveJacobian_pars(self):
        """Report whether generator has an explicit user-specified Jacobian
        with respect to pars associated with it."""
        return 'Jacobian_pars' in self.funcspec.auxfns


    def haveJacobian(self):
        """Report whether map system has an explicit user-specified Jacobian
        associated with it."""
        return 'Jacobian' in self.funcspec.auxfns


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
        """Set map system parameters"""

        validKeys = ['globalt0', 'xdomain', 'tdata', 'tdomain', 'checklevel',
                     'ics', 'pars', 'inputs', 'pdomain']
        assert remain(kw.keys(), validKeys) == [], "Invalid keys in argument"
        if 'globalt0' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, globalt0=kw['globalt0'])
        if 'checklevel' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, checklevel=kw['checklevel'])
        # optional keys for this call are ['pars', 'tdomain', 'ics',
        #   'algparams', 'tdata', 'xdomain', 'inputs', 'pdomain']
        if 'ics' in kw:
            for k_temp, v in kw['ics'].iteritems():
                k = k_temp.replace(NAMESEP, "_")
                if k in self.funcspec.vars+self.funcspec.auxvars:
                    self._xdatadict[k] = v
                else:
                    raise ValueError, 'Illegal variable name, %s'%k
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
            self._tdata[0] = self._tdomain[0]
        if self._tdomain[1] < self._tdata[1]:
            print 'tdata cannot be specified above largest '\
                  'value in tdomain\n (possibly due to uncertain bounding).'\
                  ' It has been automatically adjusted from\n ', \
                  self._tdomain[1], 'to', \
                  self._tdomain[1], '(difference of', \
                  self._tdata[1]-self._tdomain[1], ')'
            self._tdata[1] = self._tdomain[1]
        self.indepvariable.depdomain.set(self._tdata)
        if 'xdomain' in kw:
            for k_temp, v in kw['xdomain'].iteritems():
                k = k_temp.replace(NAMESEP, "_")
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
                    self._xdomain[k] = v
                else:
                    raise ValueError, 'Illegal variable name'
                try:
                    self.variables[k].depdomain.set(v)
                except TypeError:
                    raise TypeError, ('xdomain must be a dictionary of variable'
                                      ' names -> valid interval 2-tuples or '
                                      'singletons')
                for ev in self.eventstruct.events.values():
                    ev._xdomain[k] = v
        if 'pdomain' in kw:
            for k_temp, v in kw['pdomain'].iteritems():
                k = k_temp.replace(NAMESEP, "_")
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
            # pass on parameter changes to embedded system, if present
            if self._solver:
                shared_pars = intersect(kw['pars'].keys(), self._solver.pars)
                if shared_pars != []:
                    self._solver.set(pars=filteredDict(kw['pars'], shared_pars))
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
            if _inputs:
                for i in _inputs:
                    assert i in self.inputs, 'Incorrect input name provided'
                    self.inputs[i] = _inputs[i]
                # re-calc inputs ixmaps
                self._generate_ixmaps('inputs')


    # for backwards compatibility
    setPars = set


    def validateSpec(self):
        discGen.validateSpec(self)


    def compute(self, trajname, ics=None):
        assert self.funcspec.targetlang == 'python', \
               ('Wrong target language for functional specification. '
                'Python needed for this class')
        assert isinstance(self.funcspec, RHSfuncSpec), ('Map system '
                                    'requires RHSfuncSpec type to proceed')
        if self.defined:
            self.validateSpec()
##            assert remain(self.initialconditions.keys(),
##                          self._xdatadict.keys()+self.funcspec.auxvars) == [],\
##                    ('mismatching entries between initial conditions and '
##                     'declared variable names')
            self.clearWarnings()
            self.clearErrors()
        if ics is not None:
            self.set(ics=ics)
        xnames = self._var_ixmap  # ensures correct order
        # wrap up each dictionary initial value as a singleton list
        alltData = [self.indepvariable.depdomain.get(0)]
        allxDataDict = dict(zip(xnames, map(listid,
                                   sortedDictValues(self.initialconditions,
                                                    self.funcspec.vars))))
        rhsfn = eval("self."+self.funcspec.spec[1])
        # Check i.c.'s are well defined (finite)
        self.checkInitialConditions()
        self.setEventICs(self.initialconditions, self.globalt0)
        ic = sortedDictValues(self.initialconditions, self.funcspec.vars)
        plist = sortedDictValues(self.pars)
        extralist = copy(plist)
        ilist = []
        if self.inputs:
            # inputVarList is a list of Variables
            listend = self.numpars + len(self.inputs)
            inputVarList = sortedDictValues(self.inputs)
            try:
                for f in inputVarList:
                    f.clearWarnings
                    ilist.append(f(alltData[0], self.checklevel))
            except AssertionError:
                print 'External input call has t out of range: t = ', \
                    self.indepvariable.depdomain.get(0)
                print 'Maybe checklevel is 3 and initial time is not', \
                            'completely inside valid time interval'
                raise
            except ValueError:
                print 'External input call has value out of range: t = ', \
                      self.indepvariable.depdomain.get(0)
                for f in inputVarList:
                    if len(f.warnings):
                        print 'External input %s out of range:' % f.name
                        print '   t = ', repr(f.warnings[-1][0]), ', ', \
                              f.name, ' = ', repr(f.warnings[-1][1])
                raise
        else:
            listend = self.numpars
            inputVarList = []
        extralist.extend(ilist)
        precevents = self.eventstruct.query(['precise'])
        if precevents != []:
            raise PyDSTool_ValueError('precise events are not valid for map systems')
        eventslist = self.eventstruct.query(['highlevel', 'active',
                                             'notvarlinked'])
        termevents = self.eventstruct.query(['term'], eventslist)
                # initialize event info dictionaries
        Evtimes = {}
        Evpoints = {}
        for (evname, ev) in eventslist:
            Evtimes[evname] = []
            Evpoints[evname] = []
        if eventslist != []:
            self.eventstruct.resetHighLevelEvents(self.indepvariable.depdomain.get(0),
                                                  eventslist)
            self.eventstruct.validateEvents(self.funcspec.vars + \
                                            self.funcspec.auxvars + \
                                            ['t'], eventslist)

        # per-iteration storage of variable data (initial values are irrelevant)
        xDataDict = {}
        # storage of all auxiliary variable data
        allaDataDict = {}
        anames = self.funcspec.auxvars
        avals = apply(eval("self."+self.funcspec.auxspec[1]),
                      [self.indepvariable.depdomain.get(0),
                       sortedDictValues(self.initialconditions,
                                        self.funcspec.vars),
                       extralist])
        for aix in range(len(anames)):
            aname = anames[aix]
            allaDataDict[aname] = [avals[aix]]
        # temp storage of first time at which terminal events found
        # (this is used for keeping the correct end point of new mesh)
        first_found_t = None
        tmesh = self.indepvariable.depdomain.uniformSample(self.tstep, strict=False,
                                        avoidendpoints=self.checklevel>2)
        # Main loop
        breakwhile = False
        success = False
        x = ic
        notdone = True
        # did i=0 for initial condition already
        i = 1
        while notdone:
            t = tmesh[i]
            ## COMPUTE NEXT STATE y from x
            try:
                y = rhsfn(t, x, extralist)
            except:
                print "Error in calling right hand side function:"
                self.showSpec()
                raise
            for xi in xrange(self.dimension):
                xDataDict[xnames[xi]] = y[xi]
                if not self.contains(self.variables[xnames[xi]].depdomain,
                                 y[xi], self.checklevel):
                    self.warnings.append((W_TERMSTATEBD,
                                    (t, xnames[xi], y[xi],
                                     self.variables[xnames[xi]].depdomain)))
                    breakwhile = True
                    break  # for loop
            if breakwhile:
                notdone = False
                continue
            if eventslist != []:
                dataDict = copy(xDataDict)
                dataDict['t'] = t
                evsflagged = self.eventstruct.pollHighLevelEvents(None,
                                                            dataDict,
                                                            self.pars,
                                                            eventslist)
                termevsflagged = filter(lambda e: e in evsflagged, termevents)
                nontermevsflagged = filter(lambda e: e not in termevsflagged,
                                           evsflagged)
                # register any non-terminating events in the warnings list
                if len(nontermevsflagged) > 0:
                    evnames = [ev[0] for ev in nontermevsflagged]
                    self.warnings.append((W_NONTERMEVENT,
                                 (t, evnames)))
                    for evname in evnames:
                        Evtimes[evname].append(t)
                        Evpoints[evname].append(y)
                if termevsflagged != []:
                    # active terminal event flagged at this time point
                    # register the event in the warnings
                    evnames = [ev[0] for ev in termevsflagged]
                    self.warnings.append((W_TERMEVENT, \
                                             (t, evnames)))
                    for evname in evnames:
                        Evtimes[evname].append(t)
                        Evpoints[evname].append(y)
                    notdone = False
                    continue
            alltData.append(t)
            for xi in range(self.dimension):
                allxDataDict[xnames[xi]].append(y[xi])
            avals = apply(eval("self."+self.funcspec.auxspec[1]), [t,
                            sortedDictValues(xDataDict),
                            extralist])
            for aix in range(len(anames)):
                aname = anames[aix]
                allaDataDict[aname].append(avals[aix])
            try:
                extralist[self.numpars:listend] = [apply(f,
                                                [t, self.checklevel]) \
                                              for f in inputVarList]
            except ValueError:
                print 'External input call caused value out of range error:', \
                      't = ', t
                for f in inputVarList:
                    if len(f.warnings):
                        print 'External input variable %s out of range:' % f.name
                        print '   t = ', repr(f.warnings[-1][0]), ', ', \
                              f.name, ' = ', repr(f.warnings[-1][1])
                raise
            except AssertionError:
                print 'External input call caused t out of range error: t = ', t
                raise
            if i >= len(tmesh) - 1:
                notdone = False
            else:
                i += 1
                x = y
        # update success flag
        success = not notdone
        # Check that any terminal events found terminated the code correctly
        if first_found_t is not None:
            assert self.warnings[-1][0] == W_TERMEVENT, ("Event finding code "
                                        "for terminal event failed")
        # Package up computed trajectory in Variable variables
        # Add external inputs warnings to self.warnings, if any
        for f in inputVarList:
            for winfo in f.warnings:
                self.warnings.append((W_NONTERMSTATEBD,
                                     (winfo[0], f.name, winfo[1],
                                      f.depdomain)))
        # check for non-unique terminal event
        termcount = 0
        for (w,i) in self.warnings:
            if w == W_TERMEVENT or w == W_TERMSTATEBD:
                termcount += 1
                if termcount > 1:
                    self.errors.append((E_NONUNIQUETERM, (alltData[-1], i[1])))
##                print 'Time interval adjusted according to %s: %s' % \
##                      (self._warnmessages[w], str(i[0])+", "+ str(i[1]))
        # Create variables (self.variables contains no actual data)
        variables = copyVarDict(self.variables)
        # build event pointset information (reset previous trajectory's)
        self.trajevents = {}
        for (evname,  ev) in eventslist:
            evpt = Evpoints[evname]
            if evpt == []:
                self.trajevents[evname] = None
            else:
                evpt = transpose(array(evpt))
                self.trajevents[evname] = Pointset({'coordnames': xnames,
                                           'indepvarname': 't',
                                           'coordarray': evpt,
                                           'indepvararray': Evtimes[evname],
                                           'indepvartype': self.variables[xnames[0]].indepvartype})
        for x in xnames:
            if len(alltData) > 1:
                variables[x] = Variable(Pointset({'coordnames': [x],
                               'coordarray': allxDataDict[x],
                               'coordtype': self.variables[x].coordtype,
                               'indepvarname': 't',
                               'indepvararray': alltData,
                               'indepvartype': self.variables[x].indepvartype}), 't', x, x)
            else:
                raise PyDSTool_ValueError, "Fewer than 2 data points computed"
        for a in anames:
            if len(alltData) > 1:
                variables[a] = Variable(Pointset({'coordnames': [a],
                               'coordarray': allaDataDict[a],
                               'coordtype': self.variables[a].coordtype,
                               'indepvarname': 't',
                               'indepvararray': alltData,
                               'indepvartype': self.variables[a].indepvartype}), 't', a, a)
            else:
                raise PyDSTool_ValueError, "Fewer than 2 data points computed"

        if success:
            self.validateSpec()
            self.defined = True
            return Trajectory(trajname, variables.values(),
                              self.globalt0, self.checklevel)
        else:
            print 'Trajectory computation failed'
            self.errors.append((E_COMPUTFAIL, (t, self._errorcodes[errcode])))
            self.defined = False


    computeTraj = compute


    def Rhs(self, t, xdict, pdict):
        x = sortedDictValues(xdict)
        p = sortedDictValues(pdict)
        i = sortedDictValues(dict([(n,i(t)) for n,i in self.inputs.iteritems()]))
        return apply(eval("self."+self.funcspec.spec[1]), [t, x, p+i])


    def Jacobian(self, t, xdict, pdict):
        if self.haveJacobian():
            x = sortedDictValues(xdict)
            p = sortedDictValues(pdict)
            i = sortedDictValues(dict([(n,i(t)) for n,i in self.inputs.iteritems()]))
            return apply(eval("self."+self.funcspec.auxfns["Jacobian"][1]), \
                         [t, x, p+i])
        else:
            raise PyDSTool_ExistError, "Jacobian not defined"


    def JacobianP(self, t, xdict, pdict):
        if self.haveJacobian_pars():
            x = sortedDictValues(xdict)
            p = sortedDictValues(pdict)
            i = sortedDictValues(dict([(n,i(t)) for n,i in self.inputs.iteritems()]))
            return apply(eval("self."+self.funcspec.auxfns["Jacobian_pars"][1]), \
                        [t, x, p+i])
        else:
            raise PyDSTool_ExistError, "Jacobian w.r.t. parameters not defined"


    def AuxVars(self, t, xdict, pdict):
        x = sortedDictValues(xdict)
        p = sortedDictValues(pdict)
        i = sortedDictValues(dict([(n,i(t)) for n,i in self.inputs.iteritems()]))
        return apply(eval("self."+self.funcspec.auxspec[1]), [t, x, p+i])


    def __del__(self):
        discGen.__del__(self)



# Register this Generator with the database

symbolMapDict = {}
# in future, provide appropriate mappings for libraries math,
# random, etc. (for now it's left to FuncSpec)
theGenSpecHelper.add(MapSystem, symbolMapDict, 'python')
