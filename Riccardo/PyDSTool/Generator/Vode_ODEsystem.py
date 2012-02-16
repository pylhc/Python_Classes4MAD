# VODE system
from __future__ import division

from allimports import *
from PyDSTool.Generator import ODEsystem as ODEsystem
from baseclasses import Generator, theGenSpecHelper
from PyDSTool.utils import *
from PyDSTool.common import *
# PyDSTool overrides
from PyDSTool.scipy_ode import ode

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue, sign, all, any
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, less_equal, transpose, shape
import math, random
from copy import copy, deepcopy
import os, platform, shutil, sys, gc

# PSYCO FUNCTIONS DON'T WORK AS VODE CALLBACK FUNCTIONS!
#try:
#    # use pscyo JIT byte-compiler optimization, if available
#    import psyco
#    HAVE_PSYCO = True
#except ImportError:
#    HAVE_PSYCO = False
HAVE_PSYCO = False

# -----------------------------------------------------------------------------
# VERSION HISTORY
#
# 13 Sep 2006
#   Changed lastevtime to only be initialized from eventstruct when continuing
#    integration (bug fix), and made event reset state for regular forward
#    integration option None (by leaving state keyword to default) rather than
#    setting to 'off'.
#   Added support for external inputs in events.
#
# 02 Mar 2006
#   Added support for bisectlimit in event detection.
#   Improved event preparation and checking, especially involving
#     continued integrations, and fixed minor bugs.
#   Changed time mesh refinement factor to find events to 1/5 from 1/3.
#
# 22 Jan 2006
#   Added trajevents structure.
#
# October 2005
#  Minor bug fixes.
#  Overhaul of terminal event tracking.
#
# August 2005
#  Trapped occasional error when avals does not get created properly by
#    the eval statement in computeTraj.
#  Added working support for explicit Jacobian function.
#  Renamed _sys to _solver
#
# 17 July 2005
#  Fixed empty warnings added for non-terminal events.
#
# 08 July 2005
#  Moved calculation of auxiliary variable values earlier so that they could
#    be incorporated into values sent to events.
#  Events can depend on values of auxiliary variables.
#  Terminal event times & states are now included in the trajectory.
#
# -----------------------------------------------------------------------------



class Vode_ODEsystem(ODEsystem):
    """Wrapper for VODE, from SciPy.

    Uses Python functional specifications only."""

    def __init__(self, kw):
        ODEsystem.__init__(self, kw)
        assert self.funcspec.targetlang == 'python', \
               ('Wrong target language for functional specification. '
                'Python needed for this class')
        assert isinstance(self.funcspec, RHSfuncSpec), ('Vode '
                                    'requires RHSfuncSpec type to proceed')
        self._paraminfo = {'init_step': 'Fixed step size for time mesh.',
                           'strictdt': 'Boolean determining whether to evenly space time mesh (default=False), or to use exactly dt spacing.',
                           'stiff': 'Boolean to activate the BDF method, otherwise Adams method used. Default False.'}
        self._errorcodes = {0: 'Unrecognized error code returned (see stderr output)',
                -1: 'Excess work done on this call. (Perhaps wrong method MF.)',
                -2: 'Excess accuracy requested. (Tolerances too small.)',
                -3: 'Illegal input detected. (See printed message.)',
                -4: 'Repeated error test failures. (Check all input.)',
                -5: 'Repeated convergence failures. (Perhaps bad'
                    ' Jacobian supplied or wrong choice of method MF or tolerances.)',
                -6: 'Error weight became zero during problem. (Solution'
                    ' component i vanished, and ATOL or ATOL(i) = 0.)'
                }
        self._outputinfo = {'errorStatus': 'Error status on completion.'}
        # note: VODE only supports array atol, not rtol.
        algparams_def = {'rtol': 1e-9,
                         'atol': [1e-12 for dimix in xrange(self.dimension)],
                         'stiff': False,
                         'max_step': 0.0,
                         'min_step': 0.0,
                         'init_step': 0.01,
                         'max_pts': 1000000,
                         'strictdt': False
                         }
        for k, v in algparams_def.iteritems():
            if k not in self._algparams:
                self._algparams[k] = v
        self.outputstats = {}


    def addMethods(self):
        # override to add _solver function
        ODEsystem.addMethods(self, usePsyco=False)
        if self.haveJacobian():
            self._solver = ode(getattr(self,self.funcspec.spec[1]),
                            getattr(self,self.funcspec.auxfns["Jacobian"][1]))
            self._funcreg['_solver'] = ('self',
                 'ode(getattr(self,self.funcspec.spec[1]),' \
                   + 'getattr(self,self.funcspec.auxfns["Jacobian"][1]))')
        else:
            self._solver = ode(getattr(self,self.funcspec.spec[1]))
            self._funcreg['_solver'] = ('self', 'ode(getattr(self,' \
                                                  + 'self.funcspec.spec[1]))')


    def compute(self, trajname, dirn='f'):
        continue_integ = ODEsystem.prepDirection(self, dirn)
        if self._dircode == -1:
            raise NotImplementedError, ('Backwards integration is not implemented')
        # validate spec if there exists a prior trajectory computation
        if self.defined:
            self.validateSpec()
            self.validateICs()
            self.clearWarnings()
            self.clearErrors()
        pnames = sortedDictKeys(self.pars)
        xnames = self._var_ixmap  # ensures correct order
        # Check i.c.'s are well defined (finite)
        self.checkInitialConditions()
        if self._algparams['stiff']:
            methstr = 'bdf'
            methcode = 2
        else:
            methstr = 'adams'
            methcode = 1
        if self.haveJacobian():
            haveJac = 1
        else:
            haveJac = 0
        if isinstance(self._algparams['atol'], list):
            if len(self._algparams['atol']) != self.dimension:
                raise ValueError, 'atol list must have same length as phase dimension'
        else:
            atol = self._algparams['atol']
            self._algparams['atol'] = [atol for dimix in xrange(self.dimension)]
        indepdom0 = self.indepvariable.depdomain.get(0)
        indepdom1 = self.indepvariable.depdomain.get(1)
        if continue_integ:
            if self._tdata[0] != self._solver.t:
                print "Previous end time is %f"%self._solver.t
                raise ValueError, \
                      "Start time not correctly updated for continuing orbit"
            x0 = self._solver.y
            indepdom0 = self._solver.t
        else:
            x0 = sortedDictValues(self.initialconditions,
                                        self.funcspec.vars)
        if self._solver._integrator is None:
            # Banded Jacobians not yet supported
            #
            # start a new integrator, because method may have been
            # switched
            self._solver.set_integrator('vode', method=methstr,
                                 rtol=self._algparams['rtol'],
                                 atol=self._algparams['atol'],
                                 nsteps=self._algparams['max_pts'],
                                 max_step=self._algparams['max_step'],
                                 min_step=self._algparams['min_step'],
                                 first_step=self._algparams['init_step'],
                                 with_jacobian=haveJac)
            # speed up repeated access to solver by making a temp name for it
            solver = self._solver
        else:
            # speed up repeated access to solver by making a temp name for it
            solver = self._solver
            solver.with_jacobian = haveJac
#            self.mu = lband
#            self.ml = uband
            solver.rtol = self._algparams['rtol']
            solver.atol = self._algparams['atol']
            solver.method = methcode
#            self.order = order
            solver.nsteps = self._algparams['max_pts']
            solver.max_step = self._algparams['max_step']
            solver.min_step = self._algparams['min_step']
            solver.first_step = self._algparams['init_step']
        solver.set_initial_value(x0, indepdom0)
##            if self._dircode == 1:
##                solver.set_initial_value(x0, indepdom0)
##            else:
##                solver.set_initial_value(x0, indepdom1)
        # wrap up each dictionary initial value as a singleton list
        alltData = [indepdom0]
        allxDataDict = dict(zip(xnames, map(listid, x0)))
        plist = sortedDictValues(self.pars)
        extralist = copy(plist)
        if self.inputs:
            # inputVarList is a list of Variables
            inames = sortedDictKeys(self.inputs)
            listend = self.numpars + len(self.inputs)
            inputVarList = sortedDictValues(self.inputs)
            ilist = _pollInputs(inputVarList, alltData[0]+self.globalt0,
                               self.checklevel)
        else:
            ilist = []
            inames = []
            listend = self.numpars
            inputVarList = []
        extralist.extend(ilist)
        solver.set_f_params(extralist)
        if haveJac:
            solver.set_jac_params(extralist)
        dt = self._algparams['init_step']
        strict = self._algparams['strictdt']
        # Make t mesh
        if not all(isfinite(self.indepvariable.depdomain.get())):
            print "Time domain was: ", self.indepvariable.depdomain.get()
            raise ValueError, "Ensure time domain is finite"
        if dt == indepdom1 - indepdom0:
            # single-step integration required
            tmesh = [indepdom0, indepdom1]
        else:
            notDone = True
            repeatTol = 10
            count = 0
            while notDone and count <= repeatTol:
                try:
                    tmesh = self.indepvariable.depdomain.uniformSample(dt,
                                              strict=strict,
                                              avoidendpoints=self.checklevel>2)
                    notDone = False
                except AssertionError:
                    count += 1
                    dt = dt/3.0
            if count == repeatTol:
                raise AssertionError, \
                          ("supplied time step is too large for selected time"
                           " interval")
            if len(tmesh)<=2:
                # safety net, in case too few points in mesh
                # too few points unless we can add endpoint
                if tmesh[-1] != indepdom1:
                    # dt too large for tmesh to have more than one point
                    tmesh.append(indepdom1)
            if not strict:  # get actual time step used
                # don't use [0] in case avoided end points
                dt = tmesh[2]-tmesh[1]
        if self.eventstruct.query(['lowlevel']) != []:
            raise ValueError, "Only high level events can be passed to VODE"
        eventslist = self.eventstruct.query(['highlevel', 'active',
                                             'notvarlinked'])
        termevents = self.eventstruct.query(['term'], eventslist)
        # reverse time by reversing mesh doesn't work
##        if self._dircode == -1:
##            tmesh.reverse()
        tmesh.pop(0)  # get rid of first entry for initial condition
        # per-iteration storage of variable data (initial values are irrelevant)
        xDataDict = {}
        xnames = self.funcspec.vars
        # storage of all auxiliary variable data
        allaDataDict = {}
        anames = self.funcspec.auxvars
        avals = apply(getattr(self,self.funcspec.auxspec[1]),
                      [indepdom0, x0, extralist])
        for aix in range(len(anames)):
            aname = anames[aix]
            try:
                allaDataDict[aname] = [avals[aix]]
            except IndexError:
                print "\nVODE generator: There was a problem evaluating " \
                      + "an auxiliary variable"
                print "Debug info: avals (length", len(avals), ") was ", avals
                print "Index out of range was ", aix
                print self.funcspec.auxspec[1]
                print hasattr(self, self.funcspec.auxspec[1])
                print "Args were:", [indepdom0, x0, extralist]
                raise
        # Initialize signs of event detection objects at IC
        dataDict = copy(self.initialconditions)
        dataDict.update(dict(zip(anames, avals)))
        dataDict['t'] = indepdom0
        self.setEventICs(self.initialconditions, self.globalt0)
        if self.inputs:
            parsinps = copy(self.pars)
            parsinps.update(dict(zip(inames,ilist)))
        else:
            parsinps = self.pars
        if eventslist != []:
            evsflagged = self.eventstruct.pollHighLevelEvents(None,
                                                            dataDict,
                                                            parsinps,
                                                            eventslist)
            if len(evsflagged) > 0:
                raise RuntimeError, "Some events flagged at initial condition"
            if continue_integ:
                # revert to prevprevsign, since prevsign changed after call
                self.eventstruct.resetHighLevelEvents(indepdom0, eventslist, 'prev')
            else:
                self.eventstruct.resetHighLevelEvents(indepdom0, eventslist) #, 'off')
                self.eventstruct.validateEvents(self.funcspec.vars + \
                                            self.funcspec.auxvars + \
                                            self.funcspec.inputs + \
                                            ['t'], eventslist)
        # temp storage of first time at which terminal events found
        # (this is used for keeping the correct end point of new mesh)
        first_found_t = None
        # list of precise non-terminal events to be resolved after integration
        nontermprecevs = []
        evnames = [ev[0] for ev in eventslist]
        lastevtime = {}.fromkeys(evnames, None)
        # initialize new event info dictionaries
        Evtimes = {}
        Evpoints = {}
        if continue_integ:
            for evname in evnames:
                try:
                    # these are in global time, so convert to local time
                    lastevtime[evname] = self.eventstruct.Evtimes[evname][-1]-self.globalt0
                except (IndexError, KeyError):
                    # IndexError: Evtimes[evname] was None
                    # KeyError: Evtimes does not have key evname
                    pass
        for evname in evnames:
            Evtimes[evname] = []
            Evpoints[evname] = []
        # temp storage for repeatedly used object attributes (for lookup efficiency)
        depdomains = {}
        for xi in xrange(self.dimension):
            depdomains[xi] = self.variables[xnames[xi]].depdomain
        # Main integration loop
        num_points = 0
        breakwhile = False
        while not breakwhile:
            try:
                new_t = tmesh.pop(0)  # this destroys tmesh for future use
            except IndexError:
                break
            try:
                errcode = solver.integrate(new_t)
            except:
                print "Error calling right hand side function:"
                self.showSpec()
                print "Numerical traceback information (current state, " \
                      + "parameters, etc.)"
                print "in generator dictionary 'traceback'"
                self.traceback = {'vars': dict(zip(xnames,solver.y)),
                                  'pars': dict(zip(pnames,plist)),
                                  'inputs': dict(zip(inames,ilist)),
                                  self.indepvariable.name: new_t}
                raise
            for xi in xrange(self.dimension):
                xDataDict[xnames[xi]] = solver.y[xi]
                if not self.contains(depdomains[xi],
                                 solver.y[xi],
                                 self.checklevel):
                    self.warnings.append((W_TERMSTATEBD,
                                (solver.t,
                                 xnames[xi],solver.y[xi],
                                 depdomains[xi].get())))
                    breakwhile = True
                    break  # for loop
            if breakwhile:
                break
            avals = apply(getattr(self,self.funcspec.auxspec[1]), [new_t,
                            sortedDictValues(xDataDict),
                            extralist])
            # Uncomment the following assertion for debugging
#            assert all([isfinite(a) for a in avals]), \
#               "Some auxiliary variable values not finite"
            aDataDict = dict(zip(anames,avals))
            if eventslist != []:
                dataDict = copy(xDataDict)
                dataDict.update(aDataDict)
                dataDict['t'] = new_t
                if self.inputs:
                    parsinps = copy(self.pars)
                    parsinps.update(dict(zip(inames,
##                                    extralist[self.numpars:listend])))
                              _pollInputs(inputVarList, new_t+self.globalt0,
                                          self.checklevel))))
                else:
                    parsinps = self.pars
                evsflagged = self.eventstruct.pollHighLevelEvents(None,
                                                            dataDict,
                                                            parsinps,
                                                            eventslist)
##                print new_t, evsflagged
#                evsflagged = [ev for ev in evsflagged if solver.t-indepdom0 > ev[1].eventinterval]
                termevsflagged = filter(lambda e: e in evsflagged, termevents)
                nontermevsflagged = filter(lambda e: e not in termevsflagged,
                                           evsflagged)
                # register any non-terminating events in the warnings
                # list, unless they are 'precise' in which case flag
                # them to be resolved after integration completes
                if len(nontermevsflagged) > 0:
                    evnames = [ev[0] for ev in nontermevsflagged]
                    precEvts = self.eventstruct.query(['precise'],
                                                          nontermevsflagged)
                    prec_evnames = [e[0] for e in precEvts]
                    # first register non-precise events
                    nonprec_evnames = remain(evnames, prec_evnames)
                    # only record events if they have not been previously
                    # flagged within their event interval
                    if nonprec_evnames != []:
                        temp_names = []
                        for evname in nonprec_evnames:
                            prevevt_time = lastevtime[evname]
                            if prevevt_time is None:
                                ignore_ev = False
                            else:
                                if solver.t-prevevt_time < e[1].eventinterval:
                                    ignore_ev = True
                                else:
                                    ignore_ev = False
                            if not ignore_ev:
                                temp_names.append(evname)
                                lastevtime[evname] = solver.t
                        self.warnings.append((W_NONTERMEVENT,
                                     (solver.t, temp_names)))
                        for evname in temp_names:
                            Evtimes[evname].append(solver.t)
                            Evpoints[evname].append(solver.y)
                    for e in precEvts:
                        # only record events if they have not been previously
                        # flagged within their event interval
                        prevevt_time = lastevtime[e[0]]
                        if prevevt_time is None:
                            ignore_ev = False
                        else:
                            if solver.t-dt-prevevt_time < e[1].eventinterval:
                                ignore_ev = True
                            else:
                                ignore_ev = False
                        if not ignore_ev:
                            nontermprecevs.append((solver.t-dt, solver.t,e))
                            # be conservative as to where the event is, so
                            # that don't miss any events.
                            lastevtime[e[0]] = solver.t-dt
                        e[1].reset() #e[1].prevsign = None #
                do_termevs = []
                if termevsflagged != []:
                    # only record events if they have not been previously
                    # flagged within their event interval
                    for e in termevsflagged:
                        prevevt_time = lastevtime[e[0]]
##                        print "Event %s flagged."%e[0]
##                        print "  ... last time was ", prevevt_time
##                        print "  ... event interval = ", e[1].eventinterval
##                        print "  ... t = %f, dt = %f"%(solver.t, dt)
                        if prevevt_time is None:
                            ignore_ev = False
                        else:
##                            print "  ... comparison = %f < %f"%(solver.t-dt-prevevt_time, e[1].eventinterval)
                            if solver.t-dt-prevevt_time < e[1].eventinterval:
                                ignore_ev = True
##                                print "VODE ignore ev"
                            else:
                                ignore_ev = False
                        if not ignore_ev:
                            do_termevs.append(e)
                if len(do_termevs) > 0:
                    # >= 1 active terminal event flagged at this time point
                    if all([not ev[1].preciseFlag for ev in do_termevs]):
                        # then none of the events specify greater accuracy
                        # register the event in the warnings
                        evnames = [ev[0] for ev in do_termevs]
                        self.warnings.append((W_TERMEVENT, \
                                             (solver.t, evnames)))
                        for evname in evnames:
                            Evtimes[evname].append(solver.t)
                            Evpoints[evname].append(solver.y)
                        breakwhile = True  # break while loop after appending t, x
                    else:
                        # find which are the 'precise' events that flagged
                        precEvts = self.eventstruct.query(['precise'],
                                                          do_termevs)
                        # these events have flagged once so eventdelay has
                        # been used. now switch it off while finding event
                        # precisely (should be redundant after change to
                        # eventinterval and eventdelay parameters)
                        evnames = [ev[0] for ev in precEvts]
                        if first_found_t is None:
##                            print "first time round at", solver.t
                            numtries = 0
                            first_found_t = solver.t
                            restore_evts = deepcopy(precEvts)
                            minbisectlimit = min([ev[1].bisectlimit for ev in precEvts])
                            for ev in precEvts:
                                ev[1].eventdelay = 0.
                        else:
                            numtries += 1
##                            print "time round: ", numtries
                            if numtries > minbisectlimit:
                                self.warnings.append((W_BISECTLIMIT,
                                                      (solver.t, evnames)))
                                breakwhile = True
                        # find minimum eventtol in precEvts
                        dt_min = min([e[1].eventtol for e in precEvts])
                        # get previous time point
                        if len(alltData)>=1:
                            # take one step back -> told, which will
                            # get dt added back to first new meshpoint
                            # (solver.t is the step *after* the event was
                            # detected)
                            told = solver.t-dt
                        else:
                            raise ValueError("Event %s found too "%evnames[0]+\
                                 "close to local start time: try decreasing "
                                 "initial step size (current size is "
                                 "%f @ t=%f)"%(dt,solver.t+self.globalt0))
                        if dt_min >= dt:
##                            print "Registering event:", dt_min, dt
                            # register the event in the warnings
                            self.warnings.append((W_TERMEVENT,
                                             (solver.t, evnames)))
                            for evname in evnames:
                                Evtimes[evname].append(solver.t)
                                Evpoints[evname].append(solver.y)
                            # Cannot continue -- dt_min no smaller than
                            # previous dt. If this is more than the first time
                            # in this code then have found the event to within
                            # the minimum 'precise' event's eventtol, o/w need
                            # to set eventtol smaller.
                            breakwhile = True  # while loop
                        if not breakwhile:
                            dt_new = dt/5.0
                            # calc new tmesh
                            trangewidth = 2*dt #first_found_t - told
                            numpoints = int(math.ceil(trangewidth/dt_new))
                            # choose slightly smaller dt to fit trange exactly
                            dt = trangewidth/numpoints
                            tmesh = [told + i*dt for i in xrange(1, numpoints+1)]
                            # reset events according to new time mesh,
                            # setting known previous event state to be their
                            # "no event found" state
                            self.eventstruct.resetHighLevelEvents(told,
                                                                  precEvts,
                                                                  state='off')
                            # build new ic with last good values (at t=told)
                            if len(alltData)>1:
                                new_ic = [allxDataDict[xname][-1] \
                                            for xname in xnames]
                            else:
                                new_ic = x0
                            # reset integrator
                            solver.set_initial_value(new_ic, told)
                            extralist[self.numpars:listend] = \
                                    [apply(f, [told+self.globalt0,
                                               self.checklevel]) \
                                                  for f in inputVarList]
                            solver.set_f_params(extralist)
                            # continue integrating over new mesh
                            continue  # while
            alltData.append(solver.t)
            for xi in xrange(self.dimension):
                allxDataDict[xnames[xi]].append(solver.y[xi])
            for aix in xrange(len(anames)):
                aname = anames[aix]
                allaDataDict[aname].append(avals[aix])
            num_points += 1
            if not breakwhile:
                try:
                    extralist[self.numpars:listend] = [apply(f, \
                                                    [solver.t+self.globalt0,
                                                     self.checklevel]) \
                                                  for f in inputVarList]
                except ValueError:
                    print 'External input call caused value out of range error:',\
                          't = ', solver.t
                    for f in inputVarList:
                        if len(f.warnings):
                            print 'External input variable %s out of range:' % f.name
                            print '   t = ', repr(f.warnings[-1][0]), ', ', \
                                  f.name, ' = ', repr(f.warnings[-1][1])
                    raise
                except AssertionError:
                    print 'External input call caused t out of range error: t = ', \
                          solver.t
                    raise
                solver.set_f_params(extralist)
                breakwhile = not solver.successful()
        # Check that any terminal events found terminated the code correctly
        if first_found_t is not None:
            # ... then terminal events were found. Those that were 'precise' had
            # their 'eventdelay' attribute temporarily set to 0. It now should
            # be restored.
            for evname1, ev1 in termevents:
                # restore_evts are copies of the originally flagged 'precise'
                # events
                for evname2, ev2 in restore_evts:
                    if evname2 == evname1:
                        ev1.eventdelay = ev2.eventdelay
            try:
                if self.warnings[-1][0] not in [W_TERMEVENT, W_TERMSTATEBD]:
                    raise RuntimeError("Event finding code for terminal event "
                                       "failed in Generator " + self.name + \
                                       ": try decreasing eventdelay or "
                                       "eventinterval below eventtol")
            except IndexError:
                print "Output stats: ", self.outputstats
                raise RuntimeError("Event finding failed in Generator " + \
                                   self.name + ": try decreasing eventdelay "
                                   "or eventinterval below eventtol")
        # Package up computed trajectory in Variable variables
        # Add external inputs warnings to self.warnings, if any
        for f in inputVarList:
            for winfo in f.warnings:
                self.warnings.append((W_NONTERMSTATEBD,
                                     (winfo[0], f.name, winfo[1],
                                      f.depdomain.get())))
        # check for non-unique terminal event
        termcount = 0
        for (w,i) in self.warnings:
            if w == W_TERMEVENT or w == W_TERMSTATEBD:
                termcount += 1
                if termcount > 1:
                    self.errors.append((E_NONUNIQUETERM, (alltData[-1], i[1])))
        # uncomment the following lines for debugging
#        assert len(alltData) == len(allxDataDict.values()[0]) \
#             == len(allaDataDict.values()[0]), "Output data size mismatch"
#        for val_list in allaDataDict.values():
#            assert all([isfinite(x) for x in val_list])
        # Create variables (self.variables contains no actual data)
        # These versions of the variables are only final if no non-terminal
        # events need to be inserted.
        variables = copyVarDict(self.variables)
        for x in xnames:
            if len(alltData) > 1:
                variables[x] = Variable(interp1d(alltData, allxDataDict[x]),
                                             't', x, x)
            else:
                print "Error in Generator:", self.name
                print "t = ", alltData
                print "x = ", allxDataDict
                raise PyDSTool_ValueError, "Fewer than 2 data points computed"
        for a in anames:
            if len(alltData) > 1:
                variables[a] = Variable(interp1d(alltData, allaDataDict[a]),
                                        't', a, a)
            else:
                print "Error in Generator:", self.name
                print "t = ", alltData
                print "x = ", allxDataDict
                raise PyDSTool_ValueError, "Fewer than 2 data points computed"
        # Resolve non-terminal 'precise' events that were flagged, using the
        # variables created. Then, add them to a new version of the variables.
        ntpe_tdict = {}
        for (et0,et1,e) in nontermprecevs:
            lost_evt = False
            search_dt = max((et1-et0)/5,e[1].eventinterval)
            et_precise_list = e[1].searchForEvents(trange=[et0,et1],
                                              dt=search_dt,
                                              checklevel=self.checklevel,
                                              parDict=self.pars,
                                              vars=variables,
                                              inputs=self.inputs,
                                              abseps=self._abseps,
                                              eventdelay=False,
                                              globalt0=self.globalt0)
            if et_precise_list == []:
                lost_evt = True
            for et_precise in et_precise_list:
                if et_precise[0] is not None:
                    if et_precise[0] in ntpe_tdict:
                        # add event name at this time (that already exists in the dict)
                        ntpe_tdict[et_precise[0]].append(e[0])
                    else:
                        # add event name at this time (when time is not already in dict)
                        ntpe_tdict[et_precise[0]] = [e[0]]
                else:
                    lost_evt = True
            if lost_evt:
                raise PyDSTool_ExistError, \
                      ("Internal error: A non-terminal, 'precise' event '"
                       +e[0]+"' was lost after integration!")
        # add non-terminal event points to variables
        if ntpe_tdict != {}:
            # find indices of times at which event times will be inserted
            tix = 0
            evts = ntpe_tdict.keys()
            evts.sort()
            for evix in xrange(len(evts)):
                evt = evts[evix]
                evnames = ntpe_tdict[evt]
                self.warnings.append((W_NONTERMEVENT, (evt, evnames)))
                xval = [variables[x](evt) for x in xnames]
                for evname in evnames:
                    Evtimes[evname].append(evt)
                    Evpoints[evname].append(array(xval))
                tcond = less_equal(alltData[tix:], evt).tolist()
                try:
                    tix = tcond.index(0) + tix  # lowest index for t > evt
                    do_insert = (alltData[tix-1] != evt)
                except ValueError:
                    # evt = last t value so no need to add it
                    do_insert = False
                if do_insert:
                    alltData.insert(tix, evt)
                    xvaldict = dict(zip(xnames, xval))
                    for x in xnames:
                        allxDataDict[x].insert(tix, xvaldict[x])
                    for a in anames:
                        allaDataDict[a].insert(tix, variables[a](evt))
            for x in xnames:
                variables[x] = Variable(interp1d(alltData, allxDataDict[x]),
                                                 't', x, x)
            for a in anames:
                variables[a] = Variable(interp1d(alltData, allaDataDict[a]),
                                            't', a, a)
        self.outputstats = {'last_step': dt,
                            'num_fcns': num_points,
                            'num_steps': num_points,
                            'errorStatus': errcode
                            }
        if solver.successful():
            self.validateSpec()
            for evname, evtlist in Evtimes.iteritems():
                try:
                    self.eventstruct.Evtimes[evname].extend([et+self.globalt0 \
                                            for et in evtlist])
                except KeyError:
                    self.eventstruct.Evtimes[evname] = [et+self.globalt0 \
                                            for et in evtlist]
            # build event pointset information (reset previous trajectory's)
            self.trajevents = {}
            for (evname, ev) in eventslist:
                evpt = Evpoints[evname]
                if evpt == []:
                    self.trajevents[evname] = None
                else:
                    evpt = transpose(array(evpt))
                    self.trajevents[evname] = Pointset({'coordnames': xnames,
                                               'indepvarname': 't',
                                               'coordarray': evpt,
                                               'indepvararray': Evtimes[evname],
                                               'indepvartype': Float})
            self.defined = True
            return Trajectory(trajname, variables.values(),
                              self.globalt0, self.checklevel)
        else:
            try:
                self.errors.append((E_COMPUTFAIL, (solver.t,
                                              self._errorcodes[errcode])))
            except TypeError:
                # e.g. when errcode has been used to return info list
                print "Error information: ", errcode
                self.errors.append((E_COMPUTFAIL, (solver.t,
                                              self._errorcodes[0])))
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
        ODEsystem.__del__(self)


def _pollInputs(inputVarList, t, checklevel):
    ilist = []
    try:
        for f in inputVarList:
            f.clearWarnings
            ilist.append(f(t, checklevel))
    except AssertionError:
        print 'External input call has t out of range: t = ', t
        print 'Maybe checklevel is 3 and initial time is not', \
                    'completely inside valid time interval'
        raise
    except ValueError:
        print 'External input call has value out of range: t = ', t
        print 'Check beginning and end time of integration'
        for f in inputVarList:
            if len(f.warnings):
                print 'External input %s out of range:' % f.name
                print '   t = ', repr(f.warnings[-1][0]), ', ', \
                      f.name, ' = ', repr(f.warnings[-1][1])
        raise
    return ilist


# Register this Generator with the database

symbolMapDict = {}
# in future, provide appropriate mappings for libraries math,
# random, etc. (for now it's left to FuncSpec)
theGenSpecHelper.add(Vode_ODEsystem, symbolMapDict, 'python')
