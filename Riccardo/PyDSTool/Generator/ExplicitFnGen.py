# Explicit function generator
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
# 21 June 2006
#  Added reuseterms support.
#
# 22 January 2006
#  Added trajevents structure.
#
# 17 July 2005
#  Added support for singleton xdomains.
#
# -----------------------------------------------------------------------------

class ExplicitFnGen(ctsGen):
    """Explicit functional form specifying a trajectory.

    E.g. for an external input. This class allows parametric
    forms of the function, but with no dependence on x or its
    own external inputs."""

    def __init__(self, kw):
        ctsGen.__init__(self, kw)
        self.needKeys.extend(['varspecs'])
        self.optionalKeys.extend(['tdomain', 'pars', 'pdomain', 'xdomain',
                                  'ics', 'auxvars', 'vars', 'events',
                                  'fnspecs', 'tdata', 'enforcebounds',
                                  'activatedbounds', 'reuseterms'])
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
            raise PyDSTool_KeyError, 'inputs option invalid for ExplicitFnGen class'
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
        if 'reuseterms' in kw:
            fs_args['reuseterms'] = kw['reuseterms']
            self.foundKeys += 1
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
        self.funcspec = ExpFuncSpec(fs_args)
        for s in self.funcspec.spec[0]:
            if s.find('x[') > -1:
                raise ValueError('Variable values cannot depend on '
                            'other variables in explicit function specs -- '
                            'in function:\n'+s)
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
            assert self.eventstruct.getLowLevelEvents() == [], \
               "Can only pass high level events to ExplicitFnGen objects"
            assert self.eventstruct.query(['highlevel', 'varlinked']) == [], \
               "Only non-variable linked events are valid for this class"
            self.foundKeys += 1
        self.checkArgs(kw)
        if self.pars: # is defined
            self._register(self.pars)
            self.numpars = len(self.pars)
        else:
            self.numpars = 0
        assert self.funcspec.targetlang == 'python', \
               ('Wrong target language for functional specification. '
                'Python needed for this class')
        self.indepvariable = Variable(listid, Interval('t_domain', 'float',
                                              self._tdomain, self._abseps),
                             Interval('t', 'float', self._tdata,
                                      self._abseps), 't')
        self._register(self.indepvariable)
        for x in self.funcspec.vars + self.funcspec.auxvars:
            try:
                xinterval=Interval(x, 'float', self._xdomain[x], self._abseps)
            except KeyError, e:
                raise PyDSTool_KeyError, ('Mismatch between declared variables and '
                                 'xspecs: ' + str(e))
            # placeholder variable so that this class can be
            # copied before it is defined (listid function is a dummy)
            self.variables[x] = Variable(None, self.indepvariable.depdomain,
                                         xinterval, x)
        self._generate_ixmaps()
        self.dimension = len(self.funcspec.vars)


    def compute(self, trajname, ics=None):
        """Attach specification functions to callable interface."""

        assert self.funcspec.targetlang == 'python', \
               ('Wrong target language for functional specification. '
                'Python needed for this class')
        assert isinstance(self.funcspec, ExpFuncSpec), ('ExplicitFnGen'
                                    ' requires ExpFuncSpec type to proceed')
        # repeat this check (made in __init__) in case events were added since
        assert self.eventstruct.getLowLevelEvents() == [], \
               "Can only pass high level events to ExplicitFnGen objects"
        assert self.eventstruct.query(['highlevel', 'varlinked']) == [], \
               "Only non-variable linked events are valid for this class"
##        icdict_local = copy(self.initialconditions)
##        t0 = self.indepvariable.depdomain.get(0)
##        icdict_local['t'] = t0
##        have_aux = len(self.funcspec.auxvars)>0
##        for a in self.funcspec.auxvars:
##            # these functions are intended to be methods in their target
##            # Variable object, so expect a first argument 'self'
##            exec self.funcspec.auxspec[0]
##        if have_aux:
##            self.initialconditions.update(dict(zip(self.funcspec.auxvars,
##                    apply(locals()[self.funcspec.auxspec[1]],
##                                  (None, t0, sortedDictValues(icdict_local),
##                                       sortedDictValues(self.pars))) )))
        if ics is not None:
            self.set(ics=ics)
        self.setEventICs(self.initialconditions, self.globalt0)
        tempfs = deepcopy(self.funcspec)
        tempvars = copyVarDict(self.variables)
        # make unique fn for this trajectory: function definition gets executed
        # finally in Variable.addMethods() method
        tempspec = makeUniqueFn(copy(tempfs.spec[0]), 7, self.name)
        tempfs.spec = tempspec
        for x in self.funcspec.vars:
            x_ix = self.funcspec.vars.index(x)
            funcname = "_mapspecfn_" + x + "_" + timestamp(7)
            funcstr = "def " + funcname + "(self, t):\n\treturn "
            if len(self.funcspec.vars) == 1:
                # this clause is unnecessary if [0] is ever dropped
                # i.e. if spec would return plain scalar in 1D case
                funcstr += tempfs.spec[1] + "(self, t, [0], " \
                       + repr(sortedDictValues(self.pars)) + ")[0]\n"
            else:
                funcstr += tempfs.spec[1] + "(self, t, [0], " \
                       + repr(sortedDictValues(self.pars)) + ")[" \
                       + str(x_ix) + "]\n"
            tempvars[x].setOutput((funcname, funcstr), tempfs,
                                        self.globalt0, self._var_namemap,
                                        copy(self.initialconditions))
        if self.funcspec.auxvars != []:
            # make unique fn for this trajectory
            tempauxspec = makeUniqueFn(copy(tempfs.auxspec[0]), 7, self.name)
            tempfs.auxspec = tempauxspec
        for a in self.funcspec.auxvars:
            a_ix = self.funcspec.auxvars.index(a)
            funcname = "_mapspecfn_" + a + "_" + timestamp(7)
            funcstr = "def " + funcname + "(self, t):\n\treturn "
            if len(self.funcspec.auxvars) == 1:
                # this clause is unnecessary if [0] is ever dropped
                # i.e. if auxspec would return plain scalar in 1D case
                funcstr += tempfs.auxspec[1] + "(self, t, [v(t) " \
                      + "for v in self._refvars], " \
                      + repr(sortedDictValues(self.pars)) \
                      + ")[0]\n"
            else:
                funcstr += tempfs.auxspec[1] + "(self, t, [v(t) " \
                      + "for v in self._refvars], " \
                      + repr(sortedDictValues(self.pars)) \
                      + ")[" + str(a_ix) + "]\n"
            tempvars[a].setOutput((funcname, funcstr), tempfs,
                                        self.globalt0, self.funcspec.auxvars,
                                        copy(self.initialconditions),
                                        sortedDictValues(tempvars,
                                                         self.funcspec.vars))
        self.clearWarnings()
        self.clearErrors()
        # Find any events in tdomain, and adjust tdomain in case they
        # are terminal
        eventslist = self.eventstruct.query(['highlevel', 'active',
                                             'notvarlinked'])
        termevents = self.eventstruct.query(['term'], eventslist)
        Evtimes = {}
        Evpoints = {}
        for (evname, ev) in eventslist:
            Evtimes[evname] = []
            Evpoints[evname] = []
        if eventslist != []:
            for evix in xrange(len(eventslist)):
                (evname, ev) = eventslist[evix]
                evsfound = ev.searchForEvents(self.indepvariable.depdomain.get(),
                                              parDict=self.pars,
                                              vars=tempvars,
                                              checklevel=self.checklevel)
                tvals = sortedDictValues(tempvars)
                for evinfo in evsfound:
                    Evtimes[evname].append(evinfo[0])
                    Evpoints[evname].append(array([v(evinfo[0]) for v in tvals]))
            self.eventstruct.resetHighLevelEvents(self.indepvariable.depdomain.get(0),
                                                  eventslist)
            self.eventstruct.validateEvents(self.funcspec.vars + \
                                            self.funcspec.auxvars + \
                                            ['t'], eventslist)
            termevtimes = {}
            nontermevtimes = {}
            for evix in xrange(len(eventslist)):
                (evname, ev) = eventslist[evix]
                numevs = shape(Evtimes[evname])[-1]
                if numevs == 0:
                    continue
                if eventslist[evix][1].activeFlag:
                    if numevs > 1:
                        print "Event info:", Evtimes[evname]
                    assert numevs <= 1, ("Internal error: more than one "
                                     "terminal event of same type found")
                    # For safety, we should assert that this event
                    # also appears in termevents, but we don't
                    if Evtimes[evname][0] in termevtimes.keys():
                        # append event name to this warning
                        warning_ix = termevtimes[Evtimes[evname][0]]
                        self.warnings[warning_ix][1][1].append(evname)
                    else:
                        # make new termevtime entry for the new warning
                        termevtimes[Evtimes[evname][0]] = len(self.warnings)
                        self.warnings.append((W_TERMEVENT,
                                         (Evtimes[evname][0],
                                         [eventslist[evix][0]])))
                else:
                    for ev in range(numevs):
                        if Evtimes[evname][ev] in nontermevtimes.keys():
                            # append event name to this warning
                            warning_ix = nontermevtimes[Evtimes[evname][ev]]
                            self.warnings[warning_ix][1][1].append(evname)
                        else:
                            # make new nontermevtime entry for the new warning
                            nontermevtimes[Evtimes[evname][ev]] = \
                                                        len(self.warnings)
                            self.warnings.append((W_NONTERMEVENT,
                                             (Evtimes[evname][ev],
                                              [eventslist[evix][0]])))
        termcount = 0
        earliest_termtime = self.indepvariable.depdomain.get(1)
        for (w,i) in self.warnings:
            if w == W_TERMEVENT or w == W_TERMSTATEBD:
                termcount += 1
                if i[0] < earliest_termtime:
                    earliest_termtime = i[0]
        # now delete any events found after the earliest terminal event, if any
        if termcount > 0:
            warn_temp = []
            for (w,i) in self.warnings:
                if i[0] <= earliest_termtime:
                    warn_temp.append((w,i))
            self.warnings = warn_temp
#        self.indepvariable.depdomain.set([self.indepvariable.depdomain.get(0),
#                                          earliest_termtime])
        for v in tempvars.values():
            v.indepdomain.set(self.indepvariable.depdomain.get())
##                print 'Time interval adjusted according to %s: %s' % \
##                      (self._warnmessages[w], str(i[0])+", "+ str(i[1]))
        # build event pointset information (reset previous trajectory's)
        self.trajevents = {}
        for (evname, ev) in eventslist:
            evpt = Evpoints[evname]
            if evpt == []:
                self.trajevents[evname] = None
            else:
                evpt = transpose(array(evpt))
                self.trajevents[evname] = Pointset({'coordnames': sortedDictKeys(tempvars),
                                           'indepvarname': 't',
                                           'coordarray': evpt,
                                           'indepvararray': Evtimes[evname],
                                           'indepvartype': Float})
        if not self.defined:
            self._register(self.variables)
        self.validateSpec()
        self.defined = True
        return Trajectory(trajname, tempvars.values(),
                          self.globalt0, self.checklevel)


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
        """Set ExplicitFnGen parameters"""

        validKeys = ['globalt0', 'xdomain', 'tdata', 'tdomain',
                     'ics', 'pars', 'checklevel', 'pdomain']
        assert remain(kw.keys(), validKeys) == [], "Invalid keys in argument"
        if 'globalt0' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, globalt0=kw['globalt0'])
        if 'checklevel' in kw:
            # pass up to generic treatment for this
            ctsGen.set(self, checklevel=kw['checklevel'])
        # optional keys for this call are
        #   ['pars', 'tdomain', 'xdomain', 'pdomain']
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
                    raise TypeError, ('pdomain must be a dictionary of parameter'
                                      ' names -> valid interval 2-tuples or '
                                      'singletons')
                for ev in self.eventstruct.events.values():
                    ev._pdomain[k] = v
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
theGenSpecHelper.add(ExplicitFnGen, symbolMapDict, 'python', 'ExpFuncSpec')
