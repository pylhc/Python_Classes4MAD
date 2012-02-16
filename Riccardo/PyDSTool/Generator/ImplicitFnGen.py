# Implicit function generator
from __future__ import division

from allimports import *
from baseclasses import ctsGen, theGenSpecHelper
from PyDSTool.utils import *
from PyDSTool.common import *
from PyDSTool.Symbolic import ensureStrArgDict
from PyDSTool.parseUtils import NAMESEP

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
# Version History
#
# 17 July 2005
#  Added support for singleton xdomains
#
# -----------------------------------------------------------------------------

class ImplicitFnGen(ctsGen):
    """Implicitly defined functional-form trajectory generator.
    """

    def __init__(self, kw):
        ctsGen.__init__(self, kw)
        self.needKeys.extend(['varspecs', 'ics'])
        self.optionalKeys.extend(['tdomain', 'pars', 'pdomain', 'xdomain',
                                  'auxvars', 'vars', 'events',
                                  'algparams', 'fnspecs', 'tdata'])
        fs_args = {}
        if 'varspecs' in kw:
            self.foundKeys += 1
            varspecs = ensureStrArgDict(kw['varspecs'])
        else:
            raise PyDSTool_KeyError("Keyword 'varspecs' missing in argument")
        if 'tdomain' in kw:
            self._tdomain = kw['tdomain']
            self.foundKeys += 1
        else:
            self._tdomain = [-Inf, Inf]
        if 'tdata' in kw:
            self._tdata = kw['tdata']
            # _tdata is made into list to be consistent with
            # other uses of it in other Dynamical Systems...
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
        if 'inputs' in kw:
            raise PyDSTool_KeyError, 'inputs option invalid for ImplicitFnGen class'
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
        if 'ics' in kw:
            self._xdatadict = {}
            for k, v in dict(kw['ics']).iteritems():
                self._xdatadict[str(k)] = v
            # initial condition for implicit function solver algorithm
            self.initialconditions = self._xdatadict
            for name in remain(varspecs.keys(),
                               self._xdatadict.keys()):
                self.initialconditions[name] = NaN
            self.foundKeys += 1
        if 'auxvars' in kw:
            assert 'vars' not in kw, \
                   "Cannot use both 'auxvars' and 'vars' keywords"
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
            vars = remain(self._xdomain.keys(), auxvars)
        if auxvars != []:
            fs_args['auxvars'] = auxvars
        fs_args.update({'vars': vars,
                   'varspecs': varspecs,
                   'name': self.name
                   })
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
        fs_args['targetlang'] = 'python'
        self.funcspec = ImpFuncSpec(fs_args)
        for s in self.funcspec.spec[0]:
            if s.find('x[') > -1:
                raise ValueError('Variable values cannot depend on '
                            'other variables in implicit function specs -- '
                            'in function:\n'+s)
        if 'algparams' in kw:
            self._algparams = copy(kw['algparams'])
            self.foundKeys += 1
        else:
            self._algparams = {}
        if 'solvemethod' in self._algparams:
            if self._algparams['solvemethod'] not in _implicitSolveMethods:
                raise PyDSTool_ValueError, 'Invalid implicit solver type'
        # Holder and interface for events
        self.eventstruct = EventStruct()
        if 'events' in kw:
            raise PyDSTool_ValueError('ImplicitFnGen does not presently' \
                                      ' support events')
##            self._addEvents(kw['events'])
##            assert self.eventstruct.getLowLevelEvents() == [], \
##               "Can only pass high level events to ImplicitFnGen objects"
##            assert self.eventstruct.query(['highlevel', 'varlinked']) == [], \
##               "Only non-variable linked events are valid for this class"
##            self.foundKeys += 1
        self.checkArgs(kw)
        if self.pars: # is defined
            self._register(self.pars)
            self.numpars = len(self.pars)
        else:
            self.numpars = 0
        assert self.funcspec.targetlang == 'python', \
               ('Wrong target language for functional specification. '
                'Python needed for this class')
        self.newTempVars()
        self._generate_ixmaps()
        self.dimension = len(self.funcspec.vars)


    def newTempVars(self):
        self.indepvariable = Variable(listid, Interval('t_domain', 'float',
                                        self._tdomain, self._abseps),
                                      Interval('t', 'float', self._tdata,
                                        self._abseps), 't')
        if not self.defined:
            self._register(self.indepvariable)
        for x in self.funcspec.vars + self.funcspec.auxvars:
            try:
                xinterval=Interval(x, 'float', self._xdomain[x], self._abseps)
            except KeyError, e:
                raise PyDSTool_KeyError, ('Mismatch between declared variables'
                                 ' and xspecs: ' + str(e))
            # placeholder variable so that this class can be
            # copied before it is defined (listid function is a dummy)
            self.variables[x] = Variable(None, self.indepvariable.depdomain,
                                         xinterval, x)


    def compute(self, trajname, ics=None):
        """Attach specification functions to callable interface."""

        assert self.funcspec.targetlang == 'python', \
               ('Wrong target language for functional specification. '
                'Python needed for this class')
        assert isinstance(self.funcspec, ImpFuncSpec), ('ImplicitFnGen'
                                    ' requires ImpFuncSpec type to proceed')
        # repeat this check (made in __init__) in case events were added since
        assert self.eventstruct.getLowLevelEvents() == [], \
               "Can only pass high level events to ImplicitFnGen objects"
        assert self.eventstruct.query(['highlevel', 'varlinked']) == [], \
               "Only non-variable linked events are valid for this class"
        if ics is not None:
            self.set(ics=ics)
        # set some defaults for implicit function
        if 'solvemethod' in self._algparams:
            if self._algparams['solvemethod'] not in _implicitSolveMethods:
                raise PyDSTool_ValueError, 'Invalid implicit solver type'
        else:
            self._algparams['solvemethod'] = 'fsolve'
        if self._algparams['solvemethod'] in _1DimplicitSolveMethods and \
           self.dimension > 1:
            raise PyDSTool_TypeError, 'Inappropriate implicit solver for non-scalar system'
        if 'atol' not in self._algparams:
            self._algparams['atol'] = 1e-8
        if 'maxnumiter' not in self._algparams:
            self._algparams['maxnumiter'] = 100
        if self.defined:
            # reset variables
            self.newTempVars()
        self.setEventICs(self.initialconditions, self.globalt0)
        tempfs = deepcopy(self.funcspec)
        tempvars = copyVarDict(self.variables)
        # make unique fn for this trajectory
        tempspec = makeUniqueFn(copy(tempfs.spec[0]), 7, self.name)
        tempfs.spec = tempspec
        # test supplied code
        try:
            exec tempspec[0] in globals()
        except:
            print 'Error in supplied functional specification code'
            raise
        # set up implicit function: utils.makeImplicitFunction gets
        # called finally in Variable.addMethods() method.
        tempfs.algparams.update(self._algparams)
        if self.haveJacobian():
            tempfs.algparams['jac'] = self.funcspec.auxfns['Jacobian']
        else:
            tempfs.algparams['jac'] = None
        tempfs.algparams['pars'] = sortedDictValues(self.pars)
        if self.dimension == 1 and \
           self._algparams['solvemethod'] in _1DimplicitSolveMethods:
            tempfs.algparams['x0'] = sortedDictValues(self.initialconditions,
                                                  self.funcspec.vars)[0]
        else:
            tempfs.algparams['x0'] = sortedDictValues(self.initialconditions,
                                                  self.funcspec.vars)
        tempfs.algparams['impfn_name'] = "impfn_" + timestamp(7)
##        tempfs.algparams['impfn_name'] = "impfn"
        # create wrapper functions around implicit function, for each variable
        for x in self.funcspec.vars:
            x_ix = self.funcspec.vars.index(x)
            funcname = "_mapspecfn_" + x + "_" + timestamp(7)
            funcstr = "def " + funcname + "(self, t):\n\treturn " \
                    + tempfs.algparams['impfn_name'] + "(self, t)"
            if self.dimension > 1:
                funcstr += "[" + str(x_ix) + "]\n"
            else:
                funcstr += "\n"
            # make output fn of each variable the entry in the output from the
            # same implicit function.
            # initial conditions aren't needed beyond algparams['x0']
            tempvars[x].setOutput((funcname,funcstr), tempfs,
                                  self.globalt0, self._var_namemap)
        if self.funcspec.auxvars != []:
            # make unique fn for this trajectory
            tempauxspec = makeUniqueFn(copy(tempfs.auxspec[0]), 7, self.name)
            tempfs.auxspec = tempauxspec
        for a in self.funcspec.auxvars:
            a_ix = self.funcspec.auxvars.index(a)
            funcname = "_mapspecfn_" + a + "_" + timestamp(7)
            funcstr = "def " + funcname + "(self, t):\n\treturn "
            if len(self.funcspec.auxvars) == 1:
                # we'll only go through this once!
                funcstr += tempauxspec[1] + "(self, t, [v(t) " \
                      + "for v in self._refvars], " \
                      + repr(sortedDictValues(self.pars)) \
                      + ")[0]\n"
            else:
                funcstr += tempauxspec[1] + "(self, t, [v(t) " \
                      + "for v in self._refvars], " \
                      + repr(sortedDictValues(self.pars)) \
                      + ")[" + str(a_ix) + "]\n"
            # initial conditions aren't needed beyond algparams['x0']
            tempvars[a].setOutput((funcname, funcstr), tempfs,
                                        self.globalt0, self.funcspec.auxvars,
                                        None,
                                        sortedDictValues(tempvars,
                                                         self.funcspec.vars))
        self.clearWarnings()
        self.clearErrors()
        if self.eventstruct.getHighLevelEvents():
            raise PyDSTool_ValueError('ImplicitFnGen does not presently' \
                                      ' support events')
##        # Find any events in tdomain, and adjust tdomain in case they
##        # are terminal
##        eventslist = self.eventstruct.query(['highlevel', 'active',
##                                             'notvarlinked'])
##        termevents = self.eventstruct.query(['term'], eventslist)
##        if eventslist != []:
##            Evts = []
##            for evix in xrange(len(eventslist)):
##                (evname, ev) = eventslist[evix]
##                evsfound = ev.searchForEvents(self.indepvariable.depdomain.get(),
##                                              parDict=self.pars,
##                                              vars=tempvars,
##                                              checklevel=self.checklevel)
##                Evts.append([evinfo[0] for evinfo in evsfound])
##        if eventslist != []:
##            self.eventstruct.resetHighLevelEvents(self.indepvariable.depdomain.get(0),
##                                                  eventslist)
##            self.eventstruct.validateEvents(self.funcspec.vars + \
##                                            self.funcspec.auxvars + \
##                                            ['t'], eventslist)
##            termevtimes = {}
##            nontermevtimes = {}
##            for evix in xrange(len(eventslist)):
##                numevs = shape(Evts[evix])[-1]
##                if numevs == 0:
##                    continue
##                if eventslist[evix][1].activeFlag:
##                    if numevs > 1:
##                        print "Event info:", Evts[evix]
##                    assert numevs <= 1, ("Internal error: more than one "
##                                     "terminal event of same type found")
##                    # For safety, we should assert that this event
##                    # also appears in termevents, but we don't
##                    evname = eventslist[evix][0]
##                    if Evts[evix][0] in termevtimes.keys():
##                        # append event name to this warning
##                        warning_ix = termevtimes[Evts[evix][0]]
##                        self.warnings[warning_ix][1][1].append(evname)
##                    else:
##                        # make new termevtime entry for the new warning
##                        termevtimes[Evts[evix][0]] = len(self.warnings)
##                        self.warnings.append((W_TERMEVENT,
##                                         (Evts[evix][0],
##                                         [eventslist[evix][0]])))
##                else:
##                    for ev in range(numevs):
##                        if Evts[evix][ev] in nontermevtimes.keys():
##                            # append event name to this warning
##                            warning_ix = nontermevtimes[Evts[evix][ev]]
##                            self.warnings[warning_ix][1][1].append(evname)
##                        else:
##                            # make new nontermevtime entry for the new warning
##                            nontermevtimes[Evts[evix][ev]] = \
##                                                        len(self.warnings)
##                            self.warnings.append((W_NONTERMEVENT,
##                                             (Evts[evix][ev],
##                                              [eventslist[evix][0]])))
##        termcount = 0
##        earliest_termtime = self.indepvariable.depdomain.get(1)
##        for (w,i) in self.warnings:
##            if w == W_TERMEVENT or w == W_TERMSTATEBD:
##                termcount += 1
##                if i[0] < earliest_termtime:
##                    earliest_termtime = i[0]
##        # now delete any events found after the earliest terminal event, if any
##        if termcount > 0:
##            warn_temp = []
##            for (w,i) in self.warnings:
##                if i[0] <= earliest_termtime:
##                    warn_temp.append((w,i))
##            self.warnings = warn_temp
##        self.indepvariable.depdomain.set([self.indepvariable.depdomain.get(0),
##            earliest_termtime])
##        for v in tempvars.values():
##            v.indepdomain.set(self.indepvariable.depdomain.get())
####                print 'Time interval adjusted according to %s: %s' % \
####                      (self._warnmessages[w], str(i[0])+", "+ str(i[1]))
        if not self.defined:
            self._register(self.variables)
        self.validateSpec()
        self.defined = True
        return Trajectory(trajname,
                          tempvars.values(),
                          self.globalt0, self.checklevel, parameterized=False)


    # for backwards compatibility
    computeTraj = compute


    def haveJacobian_pars(self):
        """Report whether generator has an explicit user-specified Jacobian
        with respect to pars associated with it."""
        return 'Jacobian_pars' in self.funcspec.auxfns

    def haveJacobian(self):
        """Report whether generator has an explicit user-specified Jacobian
        associated with it."""
        return 'Jacobian' in self.funcspec.auxfns


    def set(self, **kw):
        """Set ImplicitFnGen parameters"""

        validKeys = ['globalt0', 'xdomain', 'tdata', 'tdomain', 'checklevel',
                     'ics', 'pars', 'algparams', 'pdomain']
        assert remain(kw.keys(), validKeys) == [], "Invalid keys in argument"
        if 'globalt0' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, globalt0=kw['globalt0'])
        if 'checklevel' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, checklevel=kw['checklevel'])
        # optional keys for this call are ['pars', 'tdomain', 'xdomain', 'pdomain']
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
                    elif type(v) in [float, int]:
                        pass
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
        if 'tdata' in kw:
            self._tdata = kw['tdata']
        if 'tdomain' in kw:
            self._tdomain = kw['tdomain']
            self.indepvariable.indepdomain.set(self._tdomain)
        if self._tdomain[0] > self._tdata[0]:
            print 'tdata cannot be specified below smallest '\
                  'value in tdomain\n (possibly due to uncertain bounding).'\
                  ' It has been automatically adjusted from\n ', self._tdata[0], \
                  'to', self._tdomain[0], '(difference of', \
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
                    raise TypeError, ('pdomain must be a dictionary of parameter'
                                      ' names -> valid interval 2-tuples or '
                                      'singletons')
        if 'ics' in kw:
            for k_temp, v in kw['ics'].iteritems():
                k = k_temp.replace(NAMESEP, "_")
                if k in self.funcspec.vars+self.funcspec.auxvars:
                    self._xdatadict[k] = v
                else:
                    raise ValueError, 'Illegal variable name'
            self.initialconditions.update(self._xdatadict)
        if 'pars' in kw:
            if not self.pars:
                raise ValueError, ('No pars were declared for this object'
                                   ' at initialization.')
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
                    raise PyDSTool_AttributeError, 'Illegal parameter name'
        if 'algparams' in kw:
            for k, v in kw['algparams'].iteritems():
                self._algparams[k] = v
        if 'solvemethod' in self._algparams:
            if self._algparams['solvemethod'] not in _implicitSolveMethods:
                raise PyDSTool_ValueError, 'Invalid implicit solver type'


    # for backwards compatibility
    setPars = set


    def validateSpec(self):
        ctsGen.validateSpec(self)
        try:
            for v in self.variables.values():
                assert isinstance(v, Variable)
            assert not self.inputs
        except AssertionError:
            print 'Invalid system specification'
            raise


    def __del__(self):
        ctsGen.__del__(self)




# Register this Generator with the database

symbolMapDict = {}
# in future, provide appropriate mappings for libraries math,
# random, etc. (for now it's left to FuncSpec)
theGenSpecHelper.add(ImplicitFnGen, symbolMapDict, 'python', 'ImpFuncSpec')
