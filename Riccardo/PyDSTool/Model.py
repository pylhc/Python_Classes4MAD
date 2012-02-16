#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      11 Sep 2006
#$Revision: 2.2.9 $
"""General purpose (hybrid) model class.

   Robert Clewley, March 2005.

A Model object's hybrid trajectory can be treated as a curve, or as
a mapping. Call the model object with the trajectory name, time(s), and
set the 'asmap' argument to be True to use an integer time to select the
trajectory segment. These are numbered from zero.

A trajectory value in a Model object's 'trajectories' dictionary
attribute is a TrajInfo object, having the following attributes:

    timeInterval is the entire time interval for the trajectory.

    timePartitions is a sequence of time_intervals (for each trajectory
        segment in trajseq), and

    trajSeq is a list of epoch or regime trajectory segments [traj_0, traj_1,
        ..., traj_(R-1)],

        where traj_i is a callable Trajectory object.

    eventStruct is the event structure used to determine that trajectory.

    genEvents is a dictionary of event names -> list of times at which that
      event took place.

    genNames is a list of the generators used for the trajectory (one per
                                                             partition).

NOTE: This file is not intended to be executable as a script -- there are no
       tests or examples at the end.
"""

# -----------------------------------------------------------------------------
# VERSION HISTORY
#
# Version 2.2.9, Sep 2006
#   Fixed sample method to keep evs_gen = None when dt is None.
#   Updated Model initialization to allow keyword arguments.
#
# Version 2.2.8, May 2006
#   Changed 'sampleTraj', 'computeTraj', and 'setParsModel' method names to
#     'sample', 'compute', and 'set', respectively, leaving the former as
#     aliases.
#   Added __getitem__ and __delitem__ methods to extract and delete
#     trajectories using reference notation.
#
# Version 2.2.7, April 2006
#   Updated hybrid curve access so that in uncertain case, trajectory is
#     called with checklevel=1 option -- prevents possible BoundsErrors from
#     interp_1d class.
#
# Version 2.2.6, March 2006
#   Minor updates for numarray 1.4.1 and 1.5.1 compatibility.
#   Updated end-user methods that involve variable and parameter names to
#     more consistently use hierarchical name convention (using '.' separator).
#
# Version 2.2.5, March 2006
#   Changed findInitialGenerator to allow a tolerance of up to gen._abseps
#     in the domain test "sign" value.
#
# Version 2.2.4, February 2006
#   Prevented setting of pars that are not present in any Generator.
#   Added use of _normord for Pointset norm order.
#   Allowed models with single generators to not require specification of
#     initial conditions when calling computeTraj (will fetch them from the
#     generator directly).
#   Renamed 'xdatadict' to 'ics' to reflect fact that Point object is a
#     valid value.
#
# Version 2.2.3, February 2006
#   Added showRegimes method to TrajInfo.
#   Minor bug fixes.
#
# Version 2.2.2, January 2006
#   Added __copy__ and __deepcopy__ methods to TrajInfo and Model classes.
#   Minor bug fixes.
#   Added MassMatrix interfacing function, which has been inadvertently
#     left out.
#   Fixed error in computeTraj when time range is altered after initial
#     trajectory computation: indepvar.depdomain -> indepvar.indepdomain in
#     test of t1_global > indepvar.indepdomain.get(1)+t0
#
# Version 2.2.1, October 2005
#   Corrected output of sampleTraj to have true hierarchically-named
#     Pointsets.
#   Model trajectory attribute is now an object of class TrajInfo.
#   In computeTraj, moved calculation of new initial condition for second and
#     subsequent hybrid traj segments to the front end of the while loop.
#   In computeTraj, changed calculation of new initial condition from old
#     state so that generators that store datapoints are accessed more quickly.
#   Simplified sampleTraj's mesh generation, and default behaviour when dt is
#     not given (defaults to underlying time mesh, when available).
#   Added 'variables' and 'auxvariables' keys to query method.
#
# Version 2.2, October 2005
#   Added TrajInfo object to replace tuple.
#   Adapted setGenEventActive() and setGenEventTerm() to accept lists of
#     generator names.
#   Optimized initial generator setup.
#   Removed DomainTest usage and replaced with direct determination of
#     initial generator to use through evaluation of each generator's
#     active terminal events at the initial conditions.
#   Set up future hybrid trajectory segment termination by multiple terminal
#     events (esp. for use in discrete time systems).
#   Minor bug fixes and checks added for debugging.
#
# Version 2.1, October 2005
#   Changed genInfo to allow an event mapping function per terminal event.
#   Renamed getTrajGenEvents() -> getTrajEvents()
#   Added support for hierarchical names.
#   Added info method.
#   Added 'mspecdict' optional initialization key, to store associated
#     ModelSpec.
#   Moved embed(), makeGenInfoEntry() and makeGenInfo() to ModelConstructor.py
#   Added module_cache feature to prevent copied generators losing loaded
#     and initialized shared library integrators (at least C-based ones).
#   Added showGenInfo() and getGenEventMappings() methods.
#   Added generators dictionary attribute for direct access to Generators.
#   Added Vfield, Jacobian, JacobianP, and AuxVars calls to access generators.
#   Added setPars() to set pars more easily (including support for
#     multiply defined pars such as p0 - p10).
#   Added searchForNames(), searchForVars(), setGenEventTerm(), and
#     setGenEventActive() methods
#
# Version 2.0.3, July 2005
#   Added generator name to trajectory structure, and associated method
#     to extract it (getTrajGenName).
#   findInitialGenerator now includes regular variables in its update of
#     icdict, but includes a continue statement for internal variables
#     that was missing before.
#   Added a trap for multiple events in sampleTraj.
#   Added a trap for empty evtlist in sampleTraj.
#   Added default termination switch rules for variable bounds.
#   Added a trap for empty tmesh in sampleTraj.
#   Changed < to <= in sampleTraj:
#                elif t <= thi:
#                    tmesh.append(thi)
#     and added code to handle evtlist only having a single entry at thi
#
# Version 2.0.2, June 2005
#   In computeTraj(), fetching of Generator's final time interval changed
#     from looking at its indepvariable.depdomain to a variable's indepdomain,
#     as Generators no longer change their own indepdomain after terminal
#     events -- they only pass on the change to their output trajectories.
#   Test for trajectory being parameterized is made. This ensures that the
#     right kind of Generator was used.
#
# Version 2.0.1, June 2005
#   Made method names more consistent with use of camel case
#   Updated sampleTraj to actually use eventData in its output
#
# ----------------------------------------------------------------------------

## PyDSTool imports
import Generator, Events #, MProject
from utils import *
from common import *
from errors import *
from Interval import *
from Trajectory import *
from Variable import *
from Points import *
from ModelSpec import *
from Symbolic import isMultiRef
from parseUtils import isHierarchicalName, NAMESEP, mapNames, symbolMapClass

## Other imports
import math, sys
from scipy import Inf, NaN, isfinite, sometrue, alltrue, sign
from numarray import array, arange, zeros, Int, Float, Complex, Int32, \
     Float64, Complex64, concatenate, transpose, shape
import copy
from time import clock

__all__ = ['Model', 'TrajInfo']

# ----------------------------------------------------------------------------

# later, make this class callable and turn it into a HybridTrajectory class
class TrajInfo(object):
    """Trajectory information container."""
    def __init__(self, timeInterval, timePartitions, trajSeq,
                 eventStruct, genEventTimes, genEvents, genNames,
                 genEventStructs):
        self.timeInterval = timeInterval
        self.timePartitions = timePartitions
        self.trajSeq = trajSeq
        self.eventStruct = eventStruct
        self.genEventTimes = genEventTimes
        self.genEvents = genEvents
        self.genNames = genNames
        self.genEventStructs = genEventStructs

    def __copy__(self):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)

    def __deepcopy__(self, memo=None, _nil=[]):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)

    def showRegimes(self):
        parts = [p[1] for p in self.timePartitions]
        res = zip(self.genNames, parts)
        for g, p in res:
            print "Regime: %s, from t = %.5f"%(g,p)

    def info(self):
        return info(self) 


class Model(object):
    """
    General-purpose Hybrid Model abstract class.
    """

    def __init__(self, *a, **kw):
        self.needKeys = ['name', 'genInfo']
        self.optionalKeys = ['events', 'ics', 'mspecdict', # 'mproject'
                             'FScompatibleNames', 'FScompatibleNamesInv',
                             'norm']
##                                'inputinfo']
        # obsvars and intvars are dictionaries model variable
        # declarations must be of the form `name: valtype`, where
        # valtype is float or int. obsvars specifies the observable
        # variables for this model, which must be present in all
        # trajectory segments (regardless of which other variables are
        # specified in those segments). intvars specifies
        # non-observable (internal) variables for this model, which
        # may or may not be present in trajectory segments, depending
        # on reductions, etc., that are present in the Generators that
        # determine the segment.
        if len(a) > 0:
            if len(a) == 1 and type(a[0])==dict:
                if intersect(a[0].keys(),kw.keys()) != []:
                    raise ValueError, \
                       ("Cannot have initialization keys common to both "
                        "dictionary and keyword arguments")
                kw.update(a[0])
            else:
                raise ValueError, ("Non-keyword arguments must be a single "
                                   "dictionary")
        try:
            self.name = kw['name']
            # dict of tuples: (Generator object, switchrulesDict,
            # epochMapStateDict) needed to build trajectories. dict
            # keys are the Generator object names
            self.genInfo = kw['genInfo']
        except KeyError:
            raise KeyError, 'Necessary keys missing from argument'
        foundKeys = len(self.needKeys)
        # by default, 'observables' are all variables that are common
        # to the Generator objects in self.genInfo, and the
        # 'internals' are those remaining generate obsvars, intvars...
        self.defaultVars()
        self.allvars = self.obsvars + self.intvars
        self.allvars.sort()
        # self.obsdimension: externally observable dimension internal,
        # full dimension (maximum possible for a Generator in this
        # model)
        self.fulldimension = len(self.allvars)
        # trajectories is a dict of trajectory segments (in sequence)
        # ... can only add one at a time!
        self.trajectories = {}
        # collect parameter info from all genInfo objects
        self.generateParamInfo()
        self.generators = {}   # registry of generators
        for infopair in self.genInfo.values():
            # (genInfo already validated in defaultVars())
            # set model tag
            infopair[0].model = self.name
            self.generators[infopair[0].name] = infopair[0]
        # FScompatibleNames is a name mapping from hierarchical names
        # using user-level '.' separator to internal '_' separator
        if 'FScompatibleNames' in kw:
            self._FScompatibleNames = kw['FScompatibleNames']
            foundKeys += 1
        else:
            self._FScompatibleNames = symbolMapClass({})
        if 'FScompatibleNamesInv' in kw:
            self._FScompatibleNamesInv = kw['FScompatibleNamesInv']
            foundKeys += 1
        else:
            self._FScompatibleNamesInv = symbolMapClass({})
        # set initial conditions if specified already (if not, they must
        # be specified before or during when compute() is called)
        self.icdict = {}
        if 'ics' in kw:
            xdatadict = dict(kw['ics'])
            for k, v in xdatadict.iteritems():
                self.icdict[self._FScompatibleNames(k)] = v
            foundKeys += 1
        if 'norm' in kw:
            self._normord = kw['norm']
        else:
            self._normord = 2
        # event structure that can be used to apply global self-consistency
        # conditions for validation of orbits
        self.eventstruct = Events.EventStruct()
##        # If part of a MProject, generate abstraction digraph info
##        self._absdigraphinfo = MProject.AbsDigraphInfo(self)
##        # If part of a MProject, keep a reference to which one
##        # this model belongs to
##        if 'mproject' in kw:
##            self._mprojecttag = kw['mproject']
##            assert isinstance(self._mprojecttag, Project.Project), ('MProject '
##                                            ' tag must be an MProject object')
##            foundKeys += 1
##        else:
##            self._mprojecttag = None
        if 'mspecdict' in kw:
            self._mspecdict = kw['mspecdict']
            foundKeys += 1
        else:
            self._mspecdict = None
        if foundKeys < len(kw):
            raise KeyError, 'Invalid keys found in arguments'
        # If not already created, a True result for
        # self.haveJacobian() means that the Model will have a
        # Jacobian method made during compute(), which
        # references the appropriate _auxfn_Jac function in the Generator
        # objects.


    def generateParamInfo(self):
        """Record parameter info locally, for future queries."""
        self.pars = {}
        for infopair in self.genInfo.values():
            try:
                self.pars.update(infopair[0].pars)
            except AttributeError:
                # no pars present in Generator
                pass


    def _makeDefaultVarNames(self):
        """Return default observable and internal variable names
        from genInfo."""
        obsvars = []
        intvars = []
        for infopair in self.genInfo.values():
            varnames = infopair[0].variables.keys()
##            varnames.extend(infopair[0].funcspec.auxvars)
            if obsvars == []:
                obsvars = varnames
            else:
                obsvars = intersect(obsvars, varnames)
                intvars.extend(remain(varnames, obsvars))
        return (obsvars, intvars)


    def haveJacobian(self):
        """Returns True iff all Generator objects in genInfo have
        defined Jacobians."""
        result = True
        for infopair in self.genInfo.values():
            result = result and infopair[0].haveJacobian()
        return result


    def haveJacobian_pars(self):
        """Returns True iff all Generator objects in genInfo have
        defined Jacobians."""
        result = True
        for infopair in self.genInfo.values():
            result = result and infopair[0].haveJacobian_pars()
        return result


    def compute(self, **kw):
        """Compute a trajectory out of constituent segments."""

        try:
            trajname = kw['trajname']
            tdata = kw['tdata']
        except KeyError:
            raise PyDSTool_KeyError, 'Invalid or missing argument keys in call to compute()'
        foundKeys = 2
        if 'verboselevel' in kw:
            if kw['verboselevel'] in [0,1,2]:
                verboselevel = kw['verboselevel']
            else:
                raise ValueError, "Invalid verbose level value"
            foundKeys += 1
        else:
            verboselevel = 0
        if 'ics' in kw:
            self.icdict = dict(kw['ics'])
            foundKeys += 1
        if 'pars' in kw:
            self.set(pars=kw['pars'])
            foundKeys += 1
        if 'force' in kw:
            force_overwrite = kw['force']
            foundKeys += 1
        else:
            force_overwrite = False
        if len(kw) != foundKeys:
            raise PyDSTool_KeyError, 'Invalid argument keys passed to compute()'
        if len(tdata) == 1:
            assert isinstance(tdata, float) or isinstance(tdata, int), \
                   'tdata must be either a single number or a 2-tuple'
            t0_global = tdata[0]
            t1_global = Inf
        elif len(tdata) == 2:
            t0_global = tdata[0]
            t1_global = tdata[1]
        else:
            raise ValueError, ('tdata argument key may be either a single '
                               'float or a 2-tuple of floats')
        # From t0 (if non-autonomous system), icdict, and switching rules
        # (self-consistency validity conditions), determine which Generator
        # applies for the initial portion of the trajectory
        assert self.genInfo != {}, 'No Generator objects defined for this model'
        if not force_overwrite:
            assert trajname not in self.trajectories, \
                                   'Trajectory name already exists'
        # Ensure all Generators provided to build trajectory share the same
        # observables, and that genInfo switch rules mappings are OK
        self._validateGenInfo(self.obsvars, self.intvars)
        # Check that all Generator pars of the same name are equal
        self._validateParameters()
        # Set initial reason for an epoch to end to None to initiate
        # search for eligible Generator for first epoch
        end_reasons = None
        # initial values
        notDone = True  # finished computing trajectory segment?
        t0 = t0_global
        # check initial condition specification in icdict
        if self.icdict == {}:
            if len(self.genInfo) == 1:
                # just get i.c. from the single generator, making sure it's non-empty
                self.icdict = copy.copy(self.generators.values()[0].initialconditions)
                if self.icdict == {}:  # still empty!
                    raise PyDSTool_ExistError, "No initial conditions specified"
            else:
                raise PyDSTool_ExistError, "No initial conditions specified"
        xdict = {}
        for xname_user, value in self.icdict.iteritems():
            xname = self._FScompatibleNames(xname_user)
            if xname not in self.obsvars+self.intvars:
                raise ValueError, ("Invalid variable name in initial "
                                   "conditions: " + xname)
            try:
                # in case xdict made up of 1D Points
                xdict[xname] = value.toarray()
            except AttributeError:
                # assume just a regular array
                xdict[xname] = value
        # flag for re-use of a generator from one hybrid segment to the next
        genReused = False
        partition_num = 0
        trajseq = []
        genNames = []
        epoch_genEvents = {}
        # reset persistent storage of event times
        for gen in self.generators.values():
            gen.eventstruct.resetEvtimes()
        while notDone:
            # find appropriate Generator to compute trajectory segment
            try:
                if end_reasons is None:
                    # first time in the while loop, so set up xnames
                    assert partition_num == 0, ("end_reasons was None on a "
                                                "non-initial epoch")
                    xnames = self.icdict.keys()
                    xnames.sort()
                    xnames = [self._FScompatibleNames(xname) for xname in xnames]
                    # only for the initial Generator of a trajectory
                    infopair = findInitialGenerator(self.genInfo, t0, xdict,
                                                    self.pars,
                                                    self.intvars,
                                                    verboselevel)
                    gen = infopair[0]
                    genSwRules = infopair[1]
                    assert remain(xnames, gen.variables.keys()) == [], \
                       ('Missing initial condition specifications')
                    for entry, val in xdict.iteritems():
                        if not isfinite(val):
                            print "Warning: %s initial condition for "%gen.name\
                                   + str(entry) + " = " + str(val)
                    if verboselevel == 2:
                        print "\nStarting partition #%i"%partition_num
                        print "Chose initial generator '%s'"%gen.name
                else:
                    # find new Generator using switching rules
                    num_reasons = len(end_reasons)
                    # num_reasons > 1 will occur most often
                    # for discrete time systems
                    if num_reasons > 1:
                        # then all next generators must be the same!
                        cand_nextGen = ""
                        mapinfos = []
                        for ri in range(num_reasons):
                            try:
                                mapinfos.append(genSwRules[end_reasons[ri]])
                            except KeyError:
                                # no rule for this termination reason
                                continue
                            nextGenName = mapinfos[-1][0]
                            if cand_nextGen == "":
                                cand_nextGen = nextGenName
                            else:
                                if cand_nextGen != nextGenName:
                                    raise RuntimeError("The multiple reasons "
                                      "for terminating last trajectory segment"
                                      " did not point to the same generator to"
                                      " continue the calculation.")
                        if len(mapinfos) == 0:
                            # No next generators found
                            nextGenName = 'terminate'
                            genEpochStateMaps = None   # will not be needed
                        else:
                            genEpochStateMaps = [m[1] for m in mapinfos]
                    else:
                        if end_reasons[0] in genSwRules:
                            mapinfo = genSwRules[end_reasons[0]]
                            nextGenName = mapinfo[0]
                            genEpochStateMaps = [mapinfo[1]]
                        else:
                            nextGenName = 'terminate'
                            genEpochStateMaps = None   # will not be needed
                    if nextGenName == 'terminate':
                        if verboselevel > 0:
                            print 'Trajectory calculation for Generator `'+gen.name\
                                  +'` terminated without a DS specified in switch'\
                                  +' rules with which to continue.'
                        notDone = False
                        continue
                    elif nextGenName == gen.name:
                        # Reuse the same generator copy (take advantage of
                        # not having to re-copy the Generator from genInfo)
                        # Currently we can't use the continue integration
                        # option of ODE solvers because of the way that each
                        # new hybrid partition is solved in "local time", i.e.
                        # from t=0. So this would mess up the solvers`
                        # 'continue' option.
                        genReused = True
                    else:
                        gen = self.genInfo[nextGenName][0]
                        genSwRules = self.genInfo[nextGenName][1]
                        genReused = False
                    if verboselevel == 2:
                        print "\nStarting partition #%i"%partition_num
                        print "Chose generator '%s'"%gen.name
            except PyDSTool_ValueError, errinfo:
                print errinfo
                raise PyDSTool_ExistError('No unique eligible Generator found:' 
                  ' cannot continue (check active terminal event definitions)')
            # compute trajectory segment until an event or t1 is reached
            if t1_global > gen.indepvariable.indepdomain.get(1)+t0:
                if isfinite(gen.indepvariable.indepdomain.get(1)):
                    if verboselevel > 0:
                        print "Warning: end time was truncated to max size of specified independent variable domain"
                    t1 = gen.indepvariable.indepdomain.get(1)+t0
                else:
                    t1 = t1_global
            else:
                t1 = t1_global
##            # prepare initial conditions
##            update_ic_vals = {}
##            for xname in remain(self.obsvars+self.intvars,self.icdict.keys()):
##                # fill in unspecified keys from _xdatadict
##                update_ic_vals[xname] = gen._xdatadict[xname]
##            xdict.update(update_ic_vals)
            # For Generators that support setPars() and
            # compute() methods, set the Generator's global time
            # reference (in case it needs it), and verify that t=0 is valid
            # in the Generator's time domain
            gen.set(tdata=[0, t1-t0], globalt0=t0)
            if partition_num > 0:
                # get new initial conditions for this epoch (partition)
                for xname in self.allvars:
                    # update final condition of epoch to be used to generate
                    # next initial condition when we re-enter the loop
                    try:
                        xdict[xname] = traj.variables[xname].output.datapoints[1][-1]
                    except KeyError:
                        # auxiliary var didn't need calling
                        pass
                    except AttributeError:
                        # non-point based output attributes of a Variable need
                        # to be called ...
                        xdict[xname] = traj.variables[xname](ti_1)
                    except PyDSTool_BoundsError:
                        print "Value out of bounds in variable call:"
                        print "  variable '%s' was called at time %f"%(xname, \
                                                             ti_1)
                        raise
                # Uncomment the following line for debugging
#                if not all([isfinite(x) for x in xdict.values()]):
#                    info(xdict, "Final state")
#                    raise RuntimeError("Final state for generator '%s' "%gen.name \
#                                   + "contained non-finite values")
                # xdict holds old Generator end state for
                # *all* variables -- map ending gen's state
                # to new gen state using (previous gen's) epoch mapping
                # and current gen's pars and external inputs
                if gen.inputs != {}:
                    # build dictionary of initial values for the
                    # external inputs
                    extinputs_ic = {}.fromkeys(gen.inputs)
                    try:
                        for k in extinputs_ic:
                            extinputs_ic[k] = gen.inputs[k](0)
                    except PyDSTool_BoundsError:
                        print "Cannot proceed - generator '" + gen.name \
                              + "' had an external input undefined " \
                              + "at t = 0"
                else:
                    extinputs_ic = {}
                genSwRules = self.genInfo[nextGenName][1]
                # insert auxiliary variable values to get correct
                # remapping of variable values
                parsinps = {}
                parsinps.update(gen.pars)
                parsinps.update(extinputs_ic)
                try:
                    xdict.update(getAuxVars(gen, 0, xdict, parsinps))
                except ValueError:
                    print "Problem evaluating auxiliary variables at state:"
                    info(xdict, "")
                    print "Auxiliary variable defs:"
                    print gen.funcspec.auxspec[0]
                    info(gen.funcspec.auxfns, "Auxiliary function defs")
                    raise
                # the following mapping might be the ID fn
                num_maps = len(genEpochStateMaps)
                # apply all maps to the initial condition
                # but check that the ordering doesn't matter (they are
                # commutative), otherwise tell the user we can't safely
                # proceed [FOR NOW, JUST DON'T ALLOW SIMULTANEOUS EVENTS]
                if num_maps > 1:
                    # ensure that all the maps are the same for multiple
                    # simultaneous events that point to the same generator
                    # (the latter was already verified above)
                    for mi in range(1,num_maps):
                        if genEpochStateMaps[0] != genEpochStateMaps[mi]:
                            raise RuntimeError("PyDSTool does not yet allow "
                                 "truly simultaneous events that do not point "
                                 "to the same generator with the same mapping")
                gen_eventstruct = gen.eventstruct
                for genEpochStateMap in genEpochStateMaps:
                    # genEpochStateMap will update xdict if necessary
##                    print "Mapping states from :", xdict
                    genEpochStateMap.evmapping(xdict, gen.pars,
                                    extinputs_ic, gen_eventstruct, t0)
##                    print " ... to :", xdict
                # only retain dynamic variables in xdict for initial conditions
                # unless we have an Explicit or Implicit generator
                if not (isinstance(gen, Generator.ExplicitFnGen) or \
                         isinstance(gen, Generator.ImplicitFnGen)):
                    for auxname in gen.funcspec.auxvars:
                        del xdict[auxname]
            # add remaining pars for system
            setup_pars = {'ics': xdict}
            try:
                if gen._algparams["init_step"] > t1-t0:
                    if verboselevel > 0:
                        print "Warning: time step too large for remaining time"\
                            + " interval. Temporarily reducing time step to " \
                            + "1/10th of its previous value"
                    setup_pars['algparams'] = {'init_step': (t1-t0)/10}
            except (AttributeError, KeyError):
                # system does not support this integration parameter
                pass
            gen.set(**setup_pars)
            gen.clearWarnings()
            gen.clearErrors()
##            print "Entering generator %s at t0=%f, t1=%f"%(gen.name, t0, t1)
            try:
                traj = gen.compute(trajname+'_'+str(partition_num))
            except PyDSTool_ValueError, e:
                print "\nError in Generator:", gen.name
                gen.showWarnings()
                gen.showErrors()
                print "Are the constituent generator absolute epsilon "
                print " tolerances too small? -- this abseps =", gen._abseps
                raise
            except:
                print "\nError in Generator:", gen.name
                self.traceback = {}
                for k,v in gen.traceback.iteritems():
                    if isinstance(v, dict):
                        dentry = {}
                        for vk, vv in v.iteritems():
                            dentry[self._FScompatibleNamesInv(vk)] = vv
                    else:
                        dentry = v
                    self.traceback[k] = dentry
                if self.traceback != {}:
                    print "Traceback dictionary copied to model"
                raise
            if not isparameterized(traj):
                raise ValueError("Generator " + gen.name + " produced a " \
                                 + " non-parameterized Trajectory")
            # update original copy of the generator with new event times during
            # this trajectory, for high level event detection to be able to check
            # eventinterval for terminal events between hybrid steps
            self.generators[gen.name].eventstruct.Evtimes = gen.eventstruct.Evtimes
##            if hasattr(gen, '_solver'):
##                # clean up memory usage after calculations in Dopri and Radau
##                try:
##                    gen._solver.CleanupEvents()
##                    gen._solver.CleanupInteg()
##                except AttributeError:
##                    # incompatible solver, so no need to worry
##                    pass
            genNames.append(gen.name)
            # look at warnings etc. for any terminating events that occurred
            if gen.errors != [] and verboselevel > 0:
                for e in gen.errors:
                    print 'Generator ' + gen.name + \
                                        ' in trajectory segment had errors'
                    gen.showErrors()
            end_reasons = ['time']
            # NB. end reason will always be 'time' for ExplicitFnGen
            # Default time interval after Generator completed,
            # in case it was truncated by a terminal event
            time_interval = traj.indepdomain
            ti_1 = time_interval.get(1)
            if gen.warnings != []:
                if verboselevel >= 1:
                    gen.showWarnings()
                for w in gen.warnings:
                    if w[0] == Generator.W_TERMEVENT or \
                       w[0] == Generator.W_TERMSTATEBD:
                        if verboselevel > 1:
                            print "Time:", ti_1+t0
                            print 'Generator ' + gen.name + \
                                ' had a terminal event. Details (in local ' + \
                                'system time) ...\n ' + str(w[1])
                        # w[1] always has t as first entry, and 
                        # either a state variable name
                        # or a list of terminal event names as second.
                        if type(w[1][1])==list:
                            if len(w[1][1]) > 1:
                                if verboselevel > 0:
                                    print 'Warning: More than one terminal event found.'
                                    print '  Consider altering event pars.'
                            end_reasons = [w[1][1][ri] \
                                   for ri in range(len(w[1][1]))]
                        else:
                            end_reasons = [w[1][1]]
##                        print "End reason: ", end_reasons[0]
            # Must scan Generator for externally-imposed active events, and
            # determine t interval if they occur before current time_interval
            # ends
##            Events.find() ??
            # Update end_reasons if externally-imposed validation events found
##            end_reasons = ??
##            time_interval.set(RESTRICTED_T_INT?)
            partition_num += 1
            if time_interval.atEndPoint(t1_global-t0, 'hi'):
                # we reached end of desired computation interval
                notDone = False
            if ti_1 > t1_global-t0:
                print "Warning: Generator time interval exceeds prescribed limits:"
                print " ... ", ti_1, " > ", t1_global-t0
                notDone = False
            # last Generator time == t1 should hold after successful traj computation
            # otherwise t1 needs to be updated according to last terminating
            # event time point
            if ti_1 < t1-t0:
                if end_reasons[0] not in genSwRules:
                    # there's prescribed time left to compute but no
                    # eligible Generator to transfer control to
                    print 'Trajectory calculation for Generator `'+gen.name+'` ' \
                          +'terminated without a Generator specified in the ' \
                          +'switching rules from which to continue.'
                    print 'Perhaps a variable went out of bounds:'
                    gen.showWarnings()
                    print 'Last reasons for stopping were: ', end_reasons
                    print "Trajectory calculated: "
                    dt_sample = (ti_1-time_interval.get(0))/10.
                    print traj.sample(dt=dt_sample)
                    notDone = False
                    raise PyDSTool_ValueError, \
                          ("Premature termination of trajectory calculation")
                elif genSwRules[end_reasons[0]][0] == 'terminate':
                    if verboselevel > 0:
                        print 'Trajectory calculation for Generator `'+gen.name+'` ' \
                              +'terminated.'
                    notDone = False
                    t1 = ti_1+t0
                else:
                    t1 = ti_1+t0
            # store trajectory segment in the sequence
            trajseq.append(traj)
            # before moving on, update t0 for next epoch / regime
            t0 = t1
            try:
                epoch_genEvents[traj.name] = gen.getEvents(globalindepvar=True)
            except AttributeError:
                # this Generator has no eventstruct
                epoch_genEvents[traj.name] = {}
        # end while loop
        # Take a copy of the events as they were at the time that the
        # trajectory was computed, for future reference
        self._addTraj(trajname, trajseq, copy.copy(self.eventstruct),
                      epoch_genEvents, genNames, force_overwrite)


    computeTraj = compute


    def query(self, querykey=''):
        """Return info about Model set-up.
        Valid query key: 'pars', 'events', 'generators',
         'ics', 'variables', 'auxvariables'"""
        assert isinstance(querykey, str), \
                       ("Query argument must be a single string")
        _keylist = ['pars', 'events', 'generators',
                    'ics', 'variables', 'auxvariables']
        if querykey not in _keylist:
            print 'Valid query keys are:', _keylist
            print "('events' key only queries model-level events, not those"
            print " inside generators)"
            raise TypeError, ('Query key '+querykey+' is not valid')
        if querykey == 'pars':
            pars = self.pars
            result = {}
            for k, v in pars.iteritems():
                result[self._FScompatibleNamesInv(k)] = v
        elif querykey == 'ics':
            icinfo = self.icdict
            result = {}
            for k, v in icinfo.iteritems():
                result[self._FScompatibleNamesInv(k)] = v
        elif querykey == 'events':
            # global model events (not bound to individual generators)
            result = {}
            result.update(self.eventstruct.events)
        elif querykey == 'generators':
            result = {}
            result.update(self.showGenInfo(returnInfo=True, verbosity=1))
        elif querykey == 'variables':
            result = self._FScompatibleNamesInv(self.icdict.keys())
        elif querykey == 'auxvariables':
            result = {}
            for genName, gen in self.generators.iteritems():
                result[genName] = self._FScompatibleNamesInv(gen.funcspec.auxvars)
        return result


    def showGenDef(self, genName=None, type=''):
        """type = 'spec', 'auxspec', 'auxfnspec', 'events', or 'modelspec'
        (leave blank for the first *four* together).

        'spec', 'auxspec' and 'auxfnspec', 'events' refer to the compiled
        target language code for the specifications. 'modelspec' refers to
        the pre-compiled abstract specifications of the model."""
        if self._mspecdict is None:
            raise PyDSTool_ExistError("Cannot use this function for models"
                                      " not defined through ModelSpec")
        if genName is None:
            print "Use showGenInfo() to find names of defined Generators"
            return
        else:
            showAll = type==''
            try:
                if type=='spec' or showAll:
                    self.generators[genName].showSpec()
                if type=='auxspec' or showAll:
                    self.generators[genName].showAuxSpec()
                if type=='auxfnspec' or showAll:
                    self.generators[genName].showAuxFnSpec()
                if type=='events' or showAll:
                    self.generators[genName].showEventSpec()
                if type=='modelspec':
                    info(self._mspecdict[genName]['modelspec'].flattenSpec(\
                                [self.generators[genName].indepvariable.name]))
            except KeyError:
                raise ValueError("Generator named %s is not known"%genName)


    def showGenInfo(self, genName=None, returnInfo=False, verbosity=0):
        """Display information about named generator"""
        if genName is None:
            if returnInfo:
                result = {}
            else:
                print "Generators defined in model " + self.name
            for infopair in self.genInfo.values():
                if returnInfo:
                    result[infopair[0].name] = infopair[0]._infostr(verbosity)
                else:
                    print infopair[0]
            if returnInfo:
                return result
        else:
            # return more information for a single generator
            try:
                result = {genName: self.genInfo[genName][0]._infostr(1)+"\n"+\
                      "Event mapping info:\n" + str(self.genInfo[genName][1])}
            except KeyError:
                raise NameError("Generator %s not found in model"%genName)
            if returnInfo:
                return result
            else:
                info(result)


    def getEventMappings(self, genName):
        try:
            return self.genInfo[genName][1]
        except KeyError:
            raise NameError("Generator %s not found in model"%genName)


    def showGenEventInfo(self, genName, verbosity=1):
        self.generators[genName].eventstruct(verbosity)


    def setGenEventTerm(self, genTarget, eventTarget, flagVal):
        if isinstance(genTarget, list):
            for genName in genTarget:
                self.generators[genName].eventstruct.setTermFlag(eventTarget,
                                                                 flagVal)
        else:
            self.generators[genTarget].eventstruct.setTermFlag(eventTarget,
                                                                 flagVal)


    def setGenEventActive(self, genTarget, eventTarget, flagVal):
        if isinstance(genTarget, list):
            for genName in genTarget:
                self.generators[genName].eventstruct.setActiveFlag(eventTarget,
                                                                 flagVal)
        else:
            self.generators[genTarget].eventstruct.setActiveFlag(eventTarget,
                                                                 flagVal)


    def _validateParameters(self):
        pars_found = {}
        for infopair in self.genInfo.values():
            try:
                pars_recognized = intersect(pars_found.keys(),
                                        infopair[0].pars)
            except AttributeError:
                # Generator doesn't have any pars
                continue  # for loop
            pars_unknown = remain(infopair[0].pars, pars_found.keys())
            for parname in pars_recognized:
                if infopair[0].pars[parname] != pars_found[parname]:
                    print "Inconsistent parameter values between Generator objects inside model"
                    print "probably means that they were changed during trajectory computation"
                    print "at discrete events: you should reset the pars when calling compute()"
                    raise ValueError, ("Inconsistency between parameter values"
                     " of same name ('%s') in different constiutent Generator"%parname +\
                     " objects: reset during call to compute")
            new_vals = dict(zip(pars_unknown, [infopair[0].pars[p] \
                                               for p in pars_unknown]))
            pars_found.update(new_vals)


    def setPars(self, p, val):
        # process multirefs first, then hierarchical names
        # calls itself recursively to resolve all names
        if isMultiRef(p):
            # e.g to change a group of numerically indexed parameter names
            # such as p0 - p9 or comp1.p0.g - comp1.p9.g
            # extract numeric range of pars
            # [ and ] are guaranteed to be present, from isMultiRef()
            lbrace = p.find('[')
            rbrace = p.find(']')
            if rbrace < lbrace:
                raise ValueError("Invalid multiple reference to pars")
            rootname = p[:lbrace]
            rangespec = p[lbrace+1:rbrace].split(',')
            try:
                remainder = p[rbrace+1:]
            except KeyError:
                # no more of p after this multireference
                remainder = ''
            if len(rangespec) != 2:
                raise ValueError("Invalid multiple reference to pars")
            loix = int(rangespec[0])
            hiix = int(rangespec[1])
            if loix >= hiix:
                raise ValueError("Invalid multiple reference to pars")
            # call setPars for each resolved name (these may include further
            # multi references or hierarchical names
            for ix in range(loix, hiix+1):
                self.setPars(rootname+str(ix), val)
        elif isHierarchicalName(p):
            if self._mspecdict is None:
                raise PyDSTool_ExistError("Cannot use this functionality for"
                                    " models not defined through ModelSpec")
            # find out: is root of p a valid 'type' in model spec?
            # find all occurrences of last p
            allFoundNames = []
            for mspecinfo in self._mspecdict.values():
                foundNames = searchModelSpec(mspecinfo['modelspec'], p)
                # don't add duplicates
                allFoundNames.extend(remain(foundNames,allFoundNames))
            # if allFoundNames == [] then either the hierarchical name was
            # a specific reference, or the type matching is invalid.
            # either way, we can just call set(p) and let that resolve the issue
            # (note that all multi-refs will have been dealt with by this point)
            if allFoundNames == []:
                self.set(pars={p: val})
            else:
                self.set(pars={}.fromkeys(allFoundNames, val))
        else:
            self.set(pars={p: val})


    def setICs(self, p, val):
        # process multirefs first, then hierarchical names
        # calls itself recursively to resolve all names
        if isMultiRef(p):
            # e.g to change a group of numerically indexed parameter names
            # such as p0 - p9 or comp1.p0.g - comp1.p9.g
            # extract numeric range of pars
            # [ and ] are guaranteed to be present, from isMultiRef()
            lbrace = p.find('[')
            rbrace = p.find(']')
            if rbrace < lbrace:
                raise ValueError("Invalid multiple reference to initial conditions")
            rootname = p[:lbrace]
            rangespec = p[lbrace+1:rbrace].split(',')
            try:
                remainder = p[rbrace+1:]
            except KeyError:
                # no more of p after this multireference
                remainder = ''
            if len(rangespec) != 2:
                raise ValueError("Invalid multiple reference to initial conditions")
            loix = int(rangespec[0])
            hiix = int(rangespec[1])
            if loix >= hiix:
                raise ValueError("Invalid multiple reference to initial conditions")
            # call setICs for each resolved name (these may include further
            # multi references or hierarchical names
            for ix in range(loix, hiix+1):
                self.setICs(rootname+str(ix), val)
        elif isHierarchicalName(p):
            if self._mspecdict is None:
                raise PyDSTool_ExistError("Cannot use this functionality for"
                                    " models not defined through ModelSpec")
            # find out: is root of p a valid 'type' in model spec?
            # find all occurrences of last p
            allFoundNames = []
            for mspecinfo in self._mspecdict.values():
                foundNames = searchModelSpec(mspecinfo['modelspec'], p)
                # don't add duplicates
                allFoundNames.extend(remain(foundNames,allFoundNames))
            # if allFoundNames == [] then either the hierarchical name was
            # a specific reference, or the type matching is invalid.
            # either way, we can just call set(p) and let that resolve the issue
            # (note that all multi-refs will have been dealt with by this point)
            if allFoundNames == []:
                self.set(ics={p: val})
            else:
                self.set(ics={}.fromkeys(allFoundNames, val))
        else:
            self.set(ics={p: val})


    def set(self, **kw):
        """Set specific parameters of Model. These will get passed on to
        all Generators that support these keys unless the restrictGenList
        argument is set (only applies to keys algparams, checklevel, and abseps).

        Permitted keys: 'pars', 'algparams', 'checklevel', 'abseps',
                        'ics'"""
        allowed_keys = ['pars', 'algparams', 'checklevel', 'abseps',
                        'ics', 'restrictGenList']
        for key in kw:
            assert key in allowed_keys, ("Not a permitted parameter argument: "
                                         "Allowed keys: "+str(allowed_keys))
        # Handle initial conditions here, because compute will pass
        # the values on to the appropriate Generators when they are called.
        if 'restrictGenList' in kw:
            restrictGenList = kw['restrictGenList']
        else:
            restrictGenList = []
        if 'ics' in kw:
            xdatadict_temp = dict(kw.pop('ics'))
            xdatadict = {}
            for k, v in xdatadict_temp.iteritems():
                xdatadict[self._FScompatibleNamesInv(k)] = v
            self.icdict.update(xdatadict)
        if restrictGenList == []:
            restrictGenList = self.generators.keys()
        # For the remaining keys, must propagate parameter changes to all
        # Generators throughout genInfo structure.
        #
        # Changed by WES 10FEB06 to handle problem of 'new' pars being
        # added if the parameter names do not exist in any generators
        genvals = self.genInfo.values()
        numgens = len(genvals)
        # loop over keywords
        for key in kw:
            # keep track of the number of errors on this keyword
            if isinstance(kw[key], dict):
                # keep track of entry errors for this key
                entry_err_attr = {}
                entry_err_val = {}
                for entrykey, entryval in kw[key].iteritems():
                    entry_err_attr[entrykey] = 0
                    entry_err_val[entrykey] = 0
            else:
                entry_err_attr = 0
                entry_err_val = 0
            # try setting pars in each generator
            for infopair in genvals:
                callparsDict = {}
                # select out the ones relevant to this Generator
                if key in infopair[0].needKeys + infopair[0].optionalKeys:
                    if key in ['algparams', 'checklevel', 'abseps'] and \
                        infopair[0].name not in restrictGenList:
                            # only apply these keys to the restricted list
                            continue
                    if isinstance(kw[key], dict):
                        for entrykey, entryval in kw[key].iteritems():
                            try:
                                infopair[0].set(**{key:{entrykey:entryval}})
                            except PyDSTool_AttributeError:
                                entry_err_attr[entrykey] += 1
                            except PyDSTool_ValueError:
                                entry_err_val[entrykey] += 1
                    else:
                        try:
                            infopair[0].set(**{key: kw[key]})
                        except PyDSTool_AttributeError:
                            entry_err_attr += 1
                        except PyDSTool_ValueError:
                            entry_err_val += 1
            # Check that none of the entries in the dictionary caused errors
            # in each generator
            if isinstance(kw[key], dict):
                for entrykey, entryval in kw[key].iteritems():
                    if entry_err_attr[entrykey] == numgens:
                        raise PyDSTool_AttributeError('Parameter does not' +\
                              ' exist in any generator') #, entrykey, entryval
                    if entry_err_val[entrykey] == numgens:
                        raise PyDSTool_ValueError('Parameter value error in' +\
                              ' every generator') #, entrykey, entryval
                    else:
                        # can't think of other ways for this error to crop up
                        pass
            else:
                if entry_err_attr == numgens:
                    raise PyDSTool_AttributeError('Parameter does not exist' +\
                             ' in any generator') #, key
                if entry_err_val == numgens:
                    raise PyDSTool_ValueError('Parameter value error in' +\
                              ' every generator') #, key
            del(entry_err_attr)
            del(entry_err_val)
        self.generateParamInfo()

    # for backwards compatibility
    setParsModel = set

    def __getitem__(self, trajname):
        if trajname in self.trajectories:
            return self.trajectories[trajname]
        else:
            raise ValueError, 'No such trajectory.'

    def __delitem__(self, trajname):
        self.delTraj(trajname)

    def delTraj(self, trajname):
        """Delete a named trajectory from the database."""
        if trajname in self.trajectories:
            del(self.trajectories[trajname])
        else:
            raise ValueError, 'No such trajectory.'


    def _addTraj(self, trajname, trajseq, eventstruct, epoch_genEvents,
                 genNames, force_overwrite=False):
        """Add a computed trajectory to database."""

        if not force_overwrite:
            assert trajname not in self.trajectories, \
                   'Trajectory name already used'
        genEvStructs = dict(zip(genNames,
                                [copy.copy(self.generators[gn].eventstruct) \
                                   for gn in genNames]))
        time_range = [None, None]
        indepvarname = ''
        genEventTimes = {}
        genEvents = {}
        for traj in trajseq:
            assert isinstance(traj, Trajectory), \
                   ('traj must contain trajectories')
            if indepvarname == '':
                indepvarname = traj.indepvarname
##            else:
##                if indepvarname != traj.indepvarname:
##                    print indepvarname, traj.indepvarname
##                    raise ValueError("Independent "
##                                        "variable name mismatch in trajectory")
            for evname in epoch_genEvents[traj.name]:
                if evname in genEvents:
                    if genEvents[evname] is not None:
                        val = mapNames(self._FScompatibleNamesInv,
                                       epoch_genEvents[traj.name][evname])
                        if val is not None:
                            genEvents[evname].append(val)
                            try:
                                genEventTimes[evname].extend(val.indepvararray.tolist())
                            except AttributeError:
                                pass
                else:
                    genEvents[evname] = mapNames(self._FScompatibleNamesInv,
                                          epoch_genEvents[traj.name][evname])
                    try:
                        genEventTimes[evname] = \
                            epoch_genEvents[traj.name][evname].indepvararray.tolist()
                    except AttributeError:
                        pass
            if time_range[0] is None:
                DS_t_numtype = traj.indepdomain.type
                if isinputcts(traj):
                    DS_t_typestr = 'continuous'
                    # continuous time must be float
                    if DS_t_numtype != Float:
                        raise TypeError, \
                         'continuous time inconsistent with non-float type'
                else:
                    DS_t_typestr = 'discrete'
                    # discrete time can be floats or ints
                DS_t_numtypestr = traj.indepdomain.typestr
                assert intersect(self.obsvars, traj.variables.keys()), \
                    'variable name for traj not present in Generator variables'
                # Pick one observable variable to test with -- they
                # must all have the same type -- it is called
                # 'varname'. This variable must be present in all Generator
                # trajectory segments. Other variables can be given in
                # those segments, but they might not be used/needed.
                varname = self.obsvars[0]
                if varname in traj.variables:
                    DS_x_numtype = traj.variables[varname].depdomain.type
                else:
                    raise ValueError, "varname not known"
##                assert DS_x_numtype == self.obsvars[varname], \
##                       ('Mismatch between declared type of variable and '
##                        'that found in trajectory segment')
                time_range = traj.indepdomain.get()
                time_range[0] += traj.globalt0
                time_range[1] += traj.globalt0
                assert len(time_range) == 2, ('time interval of Generator in '
                                            'traj cannot be singleton')
                time_partitions = [(traj.indepdomain, \
                                    traj.globalt0, \
                                    traj.checklevel)]  # initial value
            else:
                if DS_t_numtype != traj.indepdomain.type:
                    raise TypeError, \
                       'Mismatched time types for Generator sequence'
                temp = copy.copy(time_range)
                temp[1] += traj.indepdomain.get(1)
                if not traj.indepdomain.atEndPoint(time_range[1] \
                                                    - traj.globalt0, 'lo'):
                    raise ValueError('Generator sequence time intervals must '
                                     'be contiguous: ' + \
                                     str(traj.indepdomain.get(0)) + ' vs. ' + \
                                     str(time_range[1]-traj.globalt0))
                time_range[1] += traj.indepdomain.get(1)
                time_partitions.append((traj.indepdomain, \
                                        traj.globalt0, \
                                        traj.checklevel))
        # Full time interval of the trajectory
        time_interval = Interval('t', DS_t_numtypestr, time_range)
        # Add to trajectory dictionary, using hierarchical event names for
        # genEvents and genEventTimes (which are called directly by user methods
        # of Model)
        self.trajectories.update({trajname: \
                        TrajInfo(time_interval,
                                 time_partitions, trajseq,
                                 eventstruct,
                                 mapNames(self._FScompatibleNamesInv, genEventTimes),
                                 mapNames(self._FScompatibleNamesInv, genEvents),
                                 genNames,
                                 genEvStructs)
                                  })


    def __call__(self, trajname, t, varnames=None, asmap=False):
        # if asmap == True then t must be an integer in [0, #segments]
        if trajname not in self.trajectories:
            raise ValueError, ("trajectory '"+trajname+"' unknown")
        if varnames is None:
            varnames = self.obsvars   # .keys()
        else:
            if isinstance(varnames, str):
                varnames = [self._FScompatibleNames(varnames)]
            elif isinstance(varnames, list):
                varnames = self._FScompatibleNames(varnames)
            if not [varname in self.obsvars for varname in varnames]:
                raise ValueError('one or more variable names not in observable'
                             ' variable list')
        trajentry = self.trajectories[trajname]
        time_interval = trajentry.timeInterval
        # timePartitions are pairs: (interval, checklevel)
        time_partitions = trajentry.timePartitions
        trajseq = trajentry.trajSeq
        if asmap:
            # hybrid system treated as mapping -- include real time value
            # in returned point / pointset
            num_parts = len(time_partitions)
            if type(t) in [list, tuple, Array, NArray]:
                dim = len(varnames)
                result = []
                for i in range(dim+1):
                    result.append([])
                for tval in t:
                    if not isinstance(tval, int):
                        raise TypError("time values must be of integer type "
                         "when treating system as a mapping")
                    if tval < num_parts:
                        # return first point of partition t
                        part_interval = time_partitions[tval][0]
                        callvals = trajseq[tval](part_interval.get(0),
                                                 varnames, useGlobalTime=False)
                        if dim == 1:
                            result[0].append(callvals)
                        else:
                            for i in range(dim):
                                result[i].append(callvals(varnames[i]))
                        result[dim].append(part_interval.get(0) + \
                                           trajseq[tval].globalt0)
                    elif tval == num_parts:
                        # return final point of final partition
                        part_interval = time_partitions[num_parts-1][0]
                        callvals = trajseq[num_parts-1](part_interval.get(1),
                                        varnames, useGlobalTime=False)
                        if dim == 1:
                            result[0].append(callvals)
                        else:
                            for i in range(dim):
                                result[i].append(callvals(varnames[i]))
                        result[dim].append(part_interval.get(1) + \
                                           trajseq[num_parts-1].globalt0)
                    else:
                        raise ValueError, ("time values not in valid "
                                           "partition range")
                return Pointset({'coordarray': array(result),
                     'coordnames': self._FScompatibleNamesInv(varnames)+['t'],
                     'indepvartype': Int,
                     'indepvararray': t,
                     'indepvarname': 'event',
                     'norm': self._normord
                     })
            else:
                if not isinstance(t, int):
                    raise TypeError("time value must be an integer"
                                    " when treating system as a mapping")
                if t not in range(num_parts+1):
                    raise PyDSTool_BoundsError("time value not in valid"
                                                 " partition range")
                dim = len(varnames)
                if t < num_parts:
                    # return first point of partition t
                    part_interval = time_partitions[t][0]
                    result = trajseq[t](part_interval.get(0), varnames,
                               useGlobalTime=False).toarray().tolist().append( \
                              part_interval.get(0)+trajseq[t].globalt0)
                else:
                    # return final point of final partition
                    part_interval = time_partitions[num_parts-1][0]
                    result = trajseq[num_parts-1](part_interval.get(1),
                        varnames, useGlobalTime=False).toarray().tolist().append( \
                       part_interval.get(1)+trajseq[num_parts-1].globalt0)
                return Point({'coordarray': array(result),
                              'coordnames': \
                                  self._FScompatibleNamesInv(varnames) + ['t'],
                              'norm': self._normord
                             })
        else:
            # hybrid system treated as curve
            if not type(t) in [list, tuple, Array, NArray]:
                t = [t]
            retvals = []
            # lastpart_pos optimizes the re-searching of the time
            # partitions for long ranges of t
            ## CURRENTLY UNUSED
##            lastpart_pos = 0
            tix = 0
            while tix < len(t):
                # list of t vals in a partition, to allow a
                # group-call to the trajectory for that partition
                trel_part = []
                tval = t[tix]
                if time_interval.contains(tval) is notcontained:
                    print "\n** Debugging info: t value, interval, tolerance ="
                    print tval, time_interval.get(), time_interval._abseps
                    raise PyDSTool_BoundsError, \
                       ('time value outside of trajectory`s time interval '
                        'of validity (if checklevel was >=2 then endpoints '
                        'were not included in trajectory)')
                trajseq_pos = 0
                tfound = False
                for (part_interval, gt0, cl) in \
                            time_partitions:
                    trelative = tval - gt0   # a source of rounding error
                    contresult = part_interval.contains(trelative)
                    if contresult is contained:
                        # find more t vals in this partition, if any
                        trel_part.append(trelative)
                        part_done = False
                        while not part_done:
                            tix += 1
                            if tix >= len(t):
                                part_done = True
                                continue  # while loop
                            trel = t[tix] - gt0   # a source of rounding error
                            contresult_sub = part_interval.contains(trel)
                            if contresult_sub is contained:
                                trel_part.append(trel)
                            elif contresult_sub is uncertain:
                                try:
                                    dummy = trajseq[trajseq_pos](trel,
                                                        useGlobalTime=False)
                                    call_ok = True
                                except (ValueError, PyDSTool_BoundsError):
                                    # traj segment didn't accept this t value
                                    # exit if statement and increment trajseq_pos
                                    call_ok = False
                                if call_ok:
                                    trel_part.append(trel)
                                else:
                                    part_done = True
                            else:
                                part_done = True
                        tfound = True
                        break  # for loop
                    elif contresult is uncertain:
                        # first verify that this partition's trajectory
                        # segment will accept this value ... then,
                        # find more t vals in this partition, if any
                        try:
                            dummy = trajseq[trajseq_pos](trelative,
                                                         useGlobalTime=False)
                            call_ok = True
                        except (ValueError, PyDSTool_BoundsError):
                            # traj segment didn't accept this t value
                            # exit if statement and increment trajseq_pos
                            call_ok = False
                        if call_ok:
                            trel_part.append(trelative)
                            part_done = False
                            while not part_done:
                                tix += 1
                                if tix >= len(t):
                                    part_done = True
                                    continue  # while loop
                                trel = t[tix] - gt0   # a source of rounding error
                                contresult_sub = part_interval.contains(trel)
                                if contresult_sub is contained:
                                    trel_part.append(trel)
                                elif contresult_sub is uncertain:
                                    try:
                                        dummy = trajseq[trajseq_pos](trel,
                                                            useGlobalTime=False)
                                        call_ok = True
                                    except (ValueError, PyDSTool_BoundsError):
                                        # traj segment didn't accept this t value
                                        # exit if statement and increment trajseq_pos
                                        call_ok = False
                                    if call_ok:
                                        trel_part.append(trel)
                                    else:
                                        part_done = True
                                else:
                                    part_done = True
                            tfound = True
                            break  # for loop
                    trajseq_pos += 1
                if tfound:
                    # append Pointset of varname values at the t vals
                    val = trajseq[trajseq_pos](trel_part, varnames,
                                               useGlobalTime=False,
                                               checklevel=1).toarray()
                    if len(t) == 1:
                        retvals = val.getflat()
                    else:
                        retvals.append(val)
                else:
                    tinterval = self.getTrajTimeInterval(trajname)
                    print "valid t interval:", \
                          tinterval.get()
                    print "t in interval? -->", tinterval.contains(tval)
                    raise ValueError, ('t = '+str(tval)+' is either not in any '
                                       'time interval defined for this '
                                       'trajectory, or there was an out-of-range'
                                       ' value computed (see above warnings). '
                                       'Decrease step-size / event tolerance or'
                                       ' increasing abseps tolerance.')
                # No need to increment tix here. It was incremented inside
                # the inner while loop the last time a tval was checked.
                # A failed containment there left tix at the next new value
            if len(t) == 1:
                if len(varnames) == 1:
                    return retvals[0]
                else:
                    return Point({'coordarray': retvals,
                          'coordnames': self._FScompatibleNamesInv(intersect( \
                                      trajseq[trajseq_pos].varnames,
                                      varnames)),
                           'norm': self._normord})
            else:
                return Pointset({'coordarray': concatenate(retvals,1),
                         'coordnames': self._FScompatibleNamesInv(intersect( \
                                      trajseq[trajseq_pos].varnames,
                                      varnames)),
                         'indepvartype': trajseq[0].indepvartype,
                         'indepvararray': t,
                         'indepvarname': trajseq[0].indepvarname,
                         'norm': self._normord})


    def getEndPoint(self, trajname, end=1):
        xdict = {}
        endtraj = self.trajectories[trajname].trajSeq[-end]
        for xname in endtraj.varnames:
            try:
                xdict[xname] = endtraj.variables[xname].output.datapoints[1][-end]
            except KeyError:
                # auxiliary var didn't need calling
                pass
            except AttributeError:
                # non-point based output attributes of a Variable need
                # to be called ...
                tend = endtraj.indepdomain.get(end)
                xdict[xname] = endtraj.variables[xname](tend)
            except PyDSTool_BoundsError:
                print "Value out of bounds in variable call:"
                print "  variable '%s' was called at time %f"%(xname, tend)
                raise
        return Point({'coorddict': xdict,
              'coordnames': self._FScompatibleNamesInv(endtraj.varnames),
              'coordtype': Float,
              'norm': self._normord})


    def _validateVarNames(self, names):
        """Check types and uniqueness of variable names."""
        namelist = []  # records those seen so far
        for vname in names:
            assert isinstance(vname, str), 'variable name must be a string'
            assert vname not in namelist, ('variable names must be unique for'
                                           ' model')
            namelist.append(vname)


    def forceObsVars(self, varnames):
        """Force variables to be the only observables in the Model."""
        assert isinstance(varnames, list), 'argument must be a list'
        assert len(intersect(varnames, self.obsvars)) == len(varnames), \
               'All variable names in list must currently be observables'
        self.forceIntVars(remain(self.obsvars, varnames))


    def forceIntVars(self, varnames):
        """Force variables to become internal variables in the Model."""        
        assert isinstance(varnames, list), 'argument must be a list'
        FSvarnames = self._FScompatibleNames(varnames)
        isect_list = intersect(FSvarnames, self.obsvars)
        rem_list = remain(self.obsvars, FSvarnames)
        assert len(isect_list) == len(varnames), \
               'All variable names in list must currently be observables'
        for var in FSvarnames:
            del(self.obsvars[self.obsvars.index(var)])
            self.intvars.append(var)
        self.obsvars.sort()
        self.intvars.sort()
        self.obsdimension = len(self.obsvars)


    def defaultVars(self):
        """(Re)set to default observable and internal variable names."""
        obsvars, intvars = self._makeDefaultVarNames()
        self._validateVarNames(obsvars + intvars)
        self._validateGenInfo(obsvars, intvars)
        # OK to store these permanently after varname validation
        self.obsvars = obsvars
        self.intvars = intvars
        self.obsvars.sort()
        self.intvars.sort()
        self.obsdimension = len(self.obsvars)


    def _infostr(self, verbose=1):
        if verbose > 0:
            outputStr = 'Model '+self.name+" containing generators:"
            for genName, infopair in self.genInfo.iteritems():
                outputStr += "\n  "+infopair[0]._infostr(verbose-1)
        else:
            outputStr = 'Model '+self.name
        return outputStr


    def info(self, verboselevel=1):
        print self._infostr(verboselevel)


    def __repr__(self):
        return self._infostr(verbose=0)

    __str__ = __repr__


    def __copy__(self):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


    def __deepcopy__(self, memo=None, _nil=[]):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


    def getTrajGenName(self, trajname, t=None):
        """Return the named trajectory's associated Generator(s), specific
        to time t if given."""
        try:
            genNames = self.trajectories[trajname].genNames
        except KeyError:
            raise ValueError, 'No such trajectory name'
        if t is None:
            # return list of Generators for associated hybrid trajectory
            return genNames
        else:                
            parts = self.getTrajTimePartitions(trajname)
            pix = 0
            for pinterval in parts:
                if pinterval.contains(t) is not notcontained:
                    return genNames[pix]
                else:
                    pix += 1


    def getTrajEventTimes(self, trajname, events=None):
        """Return the named trajectory's Generator-flagged event times.

        events argument can be singleton string name of an event,
        returning the event data, or events can be a list of event names,
        returning a dictionary of event name -> event data.
        Event names should use hierarchical naming convention, if
        applicable."""
        try:
            if events is None:
                return self.trajectories[trajname].genEventTimes
            elif isinstance(events, str):
                return self.trajectories[trajname].genEventTimes[events]
            elif isinstance(events, list):
                for evname in events:
                    result[evname] = \
                             self.trajectories[trajname].genEventTimes[evname]
                return result
        except KeyError:
            if trajname in self.trajectories:
                raise ValueError, 'No such events found'
            else:
                raise ValueError, 'No such trajectory name'


    def getTrajEvents(self, trajname, events=None):
        """Return the named trajectory's Generator-flagged events.

        events argument can be singleton string name of an event,
        returning the event data, or events can be a list of event names,
        returning a dictionary of event name -> event data.
        Event names should use hierarchical naming convention, if
        applicable."""
        try:
            if events is None:
                return self.trajectories[trajname].genEvents
            elif isinstance(events, str):
                return self.trajectories[trajname].genEvents[events]
            elif isinstance(events, list):
                for evname in events:
                    result[evname] = \
                             self.trajectories[trajname].genEvents[evname]
                return result
        except KeyError:
            if trajname in self.trajectories:
                raise ValueError, 'No such events found'
            else:
                raise ValueError, 'No such trajectory name'


    def getTrajEventStruct(self, trajname):
        """Return the named trajectory's model event structure (representing
        external constraints), as present when it was computed.
        
        (CURRENTLY UNUSED)"""
        try:
            return self.trajectories[trajname].eventStruct
        except KeyError:
            raise ValueError, 'No such trajectory name'


    def getTrajTimeInterval(self, trajname):
        """Return the named trajectory's time domain,
        over which it is defined, in a single interval."""
        try:
            return self.trajectories[trajname].timeInterval
        except KeyError:
            raise ValueError, 'No such trajectory name'


    def numPartitions(self, trajname):
        """Return the number of partitions in a named hybrid trajectory."""
        try:
            return len(self.trajectories[trajname].timePartitions)
        except KeyError:
            raise ValueError, 'No such trajectory name'


    def getTrajTimePartitions(self, trajname):
        """Return the named trajectory's time domain,
        over which it is defined, as a list of time interval partitions."""
        try:
            return self.trajectories[trajname].timePartitions
        except KeyError:
            raise ValueError, 'No such trajectory name'


    def sample(self, trajname, varlist=None, dt=None,
                   tlo=None, thi=None, events=None, precise=False):
        """Uniformly sample the named trajectory over range indicated,
        including any event points. (e.g. use this call for plotting purposes.)

        precise=False attempts to use the underlying mesh of the trajectory
        to return a Pointset more quickly. Currently, this can only be used
        for trajectories that have a single segment.

        If dt is not given, the underlying time mesh is used, if available.
        """

        if varlist is None:
            varlist = self.obsvars
        else:
            varlist = self._FScompatibleNames(varlist)
        assert isinstance(varlist, list), 'varlist argument must be a list'
        t_interval = self.getTrajTimeInterval(trajname)
        if tlo is None:
            tlo = t_interval.get(0)
        if thi is None:
            thi = t_interval.get(1)
        assert tlo < thi, 't start point must be less than t endpoint'
        if dt is None:
            if events is None:
                evs_gen = None
            else:
                evs_gen = self.getTrajEventTimes(trajname, events)
        else:
            assert dt < abs(thi-tlo), 'dt must be smaller than time interval'
            if events is None:
                evs_gen = self.getTrajEventTimes(trajname)
            else:
                evs_gen = self.getTrajEventTimes(trajname, events)
##        evs_mod_struct = self.getTrajEventStruct(trajname)
        # evs_mod = ?
##        if evs_mod or evs_gen:
        if precise:
            if dt is None:
                raise ValueError("Must specify an explicit dt when precise flag"
                                 " is set")
            tmesh = concatenate([arange(tlo, thi, dt), [thi]])
            if evs_gen:
                # combine event data dictionaries
                evs = copy.copy(evs_gen)
##                if evs_mod != {}:
##                    evs.update(dict(evs_mod))
                # sort into ordered list of times
                evtlist = orderEventData(evs, nonames=True)
                if evtlist != []:
                    tmesh = insertInOrder(tmesh, evtlist)
            if len(tmesh) > 0:
                pset = self(trajname,tmesh,varlist)
            else:
                return None
        else:
            #FSvarlist = self._FScompatibleNames(varlist)
            pset = None
            for traj in self.trajectories[trajname].trajSeq:
                # use fact that Trajectory.sample automatically truncates
                # tlo and thi if they are beyond that segment's defined range.
                if pset is None:
                    pset = traj.sample(varlist, dt, tlo, thi, evs_gen,
                                           False)
                else:
                    pset_new = traj.sample(varlist, dt, tlo, thi,
                                                evs_gen, False)
                    pset.append(pset_new, skipMatchingIndepvar=True)
        pset.indepvarname = self._FScompatibleNamesInv(pset.indepvarname)
        pset.coordnames = self._FScompatibleNamesInv(pset.coordnames)
        pset.makeIxMaps()
        return pset


    sampleTraj = sample


    def searchForNames(self, template):
        """Find parameter / variables names matching template in
        defined generators."""
        if self._mspecdict is None:
            raise PyDSTool_ExistError("Cannot use this function for models"
                                      " not defined through ModelSpec")
        result = {}
        for genName, mspecinfo in self._mspecdict.iteritems():
            foundNames = searchModelSpec(mspecinfo['modelspec'], template)
            result[genName] = foundNames
        return result


    def searchForVars(self, template):
        """Find variable and auxiliary variable names that have to be in every
        generator: the result is a list."""
        if self._mspecdict is None:
            raise PyDSTool_ExistError("Cannot use this function for models"
                                      " not defined through ModelSpec")
        # all generators must define the same variables, so need only
        # to look at one
        a_gen_name = self._mspecdict.keys()[0]
        a_gen_mspec = self._mspecdict[a_gen_name]['modelspec']
        foundNames = searchModelSpec(a_gen_mspec, template)
        fspec = self.generators[a_gen_name].funcspec
        return intersect(self._FScompatibleNamesInv(fspec.vars + fspec.auxvars),
                             foundNames)


    def _validateGenInfo(self, obsvars, intvars):
        """Validate Model's genInfo attribute."""

        allgen_names = self.genInfo.keys()
        for genName, infopair in self.genInfo.iteritems():
            # generator for trajectory
            generator = infopair[0]
            # dict of reason -> (gen_name, epmapping) pairs 
            swmapping = infopair[1]
            allvarnames = generator.variables.keys()
            special_reasons = ['time'] + allvarnames
            assert generator.name == genName, ('Generator`s name does not '
                                'correspond to the key in the genInfo dict')
            assert generator.name in allgen_names, ('Generator`s name not in list '
                                        'of all available Generator names!')
            assert len(intersect(obsvars, allvarnames)) ==\
                len(obsvars), 'Generator `'+generator.name+'` did not contain ' \
                                      + ' required observable variables'
            assert len(intersect(intvars, allvarnames)) ==\
                len(intvars), 'Generator '+generator.name+' did not contain ' \
                                      + ' required internal variables'
            try:
                allEndReasonNames = map(lambda (n,o):n,
                                       generator.eventstruct.getTermEvents()) \
                                  + special_reasons
            except AttributeError:
                allEndReasonNames = special_reasons
            seenReasons = []
            assert len(swmapping) == len(allEndReasonNames), ('Incorrect number'
                                            ' of map pairs given in argument')
            for (reason, outcome) in swmapping.iteritems():
                genTargetName = outcome[0]
                epmapping = outcome[1]
                # no checks are made on epmapping here
                assert reason not in seenReasons, ('reasons cannot appear more'
                                            ' than once in map domain')
                seenReasons.append(reason)
                assert reason in allEndReasonNames, ('name `'+reason+'` in map '
                                                    'domain is invalid')
                if genTargetName != 'terminate':
                    assert genTargetName in allgen_names, ('name `'+genTargetName+\
                                            '` in map range is invalid')


    def Rhs(self, genName, t, xdict, pdict, idict={}):
        """Direct access to a generator's Rhs function."""
        try:
            gen = self.genInfo[genName][0]
        except KeyError:
            raise ValueError("No Generator %s was found"%genName)
        xdict_compat = {}
        for k, v in xdict.iteritems():
            xdict_compat[self._FScompatibleNames(k)] = v
        pdict_compat = {}
        for k, v in pdict.iteritems():
            pdict_compat[self._FScompatibleNames(k)] = v        
        idict_compat = {}
        for k, v in idict.iteritems():
            idict_compat[self._FScompatibleNames(k)] = v
        return Point({'coorddict': dict(zip(self._FScompatibleNamesInv(gen.funcspec.vars),
                                gen.Rhs(t, xdict_compat, pdict_compat, idict_compat))),
                      'coordtype': Float,
                      'norm': self._normord})


    def MassMatrix(self, genName, t, xdict, pdict, idict={}):
        """Direct access to a generator's MassMatrix function (if defined)."""
        try:
            gen = self.genInfo[genName][0]
        except KeyError:
            raise ValueError("No Generator %s was found"%genName)
        xdict_compat = {}
        for k, v in xdict.iteritems():
            xdict_compat[self._FScompatibleNames(k)] = v
        pdict_compat = {}
        for k, v in pdict.iteritems():
            pdict_compat[self._FScompatibleNames(k)] = v        
        idict_compat = {}
        for k, v in idict.iteritems():
            idict_compat[self._FScompatibleNames(k)] = v
        return Point({'coorddict': dict(zip(self._FScompatibleNamesInv(gen.funcspec.vars),
                            gen.MassMatrix(t, xdict_compat, pdict_compat, idict_compat))),
                      'coordtype': Float,
                      'norm': self._normord})


    def Jacobian(self, genName, t, xdict, pdict, idict={}):
        """Direct access to a generator's Jacobian function (if defined)."""
        try:
            gen = self.genInfo[genName][0]
        except KeyError:
            raise ValueError("No Generator %s was found"%genName)
        if gen.haveJacobian():
            xdict_compat = {}
            for k, v in xdict.iteritems():
                xdict_compat[self._FScompatibleNames(k)] = v
            pdict_compat = {}
            for k, v in pdict.iteritems():
                pdict_compat[self._FScompatibleNames(k)] = v        
            idict_compat = {}
            for k, v in idict.iteritems():
                idict_compat[self._FScompatibleNames(k)] = v
            # would prefer to be returning the 'D' attached to the lowest name
            # in a hierarchical name
            dvars = ['D'+v for v in self._FScompatibleNamesInv(gen.funcspec.vars)]
            return Pointset({'coorddict': dict(zip(dvars, gen.Jacobian(t,
                                xdict_compat, pdict_compat, idict_compat))),
                             'coordtype': Float,
                             'norm': self._normord})
        else:
            raise PyDSTool_ExistError, "Jacobian not defined"


    def JacobianP(self, genName, t, xdict, pdict, idict={}):
        """Direct access to a generator's JacobianP function (if defined)."""
        try:
            gen = self.genInfo[genName][0]
        except KeyError:
            raise ValueError("No Generator %s was found"%genName)
        if gen.haveJacobian_pars():
            xdict_compat = {}
            for k, v in xdict.iteritems():
                xdict_compat[self._FScompatibleNames(k)] = v
            pdict_compat = {}
            for k, v in pdict.iteritems():
                pdict_compat[self._FScompatibleNames(k)] = v        
            idict_compat = {}
            for k, v in idict.iteritems():
                idict_compat[self._FScompatibleNames(k)] = v
            # would prefer to be returning the 'D' attached to the lowest name
            # in a hierarchical name                
            dpars = ['D'+p for p in self._FScompatibleNamesInv(gen.funcspec.pars)]
            return Pointset({'coorddict': dict(zip(dpars, gen.JacobianP(t,
                                    xdict_compat, pdict_compat, idict_compat))),
                             'coordtype': Float,
                             'norm': self._normord})
        else:
            raise PyDSTool_ExistError, "Jacobian w.r.t. pars not defined"


    def AuxVars(self, genName, t, xdict, pdict, idict={}):
        """Direct access to a generator's auxiliary variables
        definition (if defined)."""
        try:
            gen = self.genInfo[genName][0]
        except KeyError:
            raise ValueError("No Generator %s was found"%genName)
        xdict_compat = {}
        for k, v in xdict.iteritems():
            xdict_compat[self._FScompatibleNames(k)] = v
        pdict_compat = {}
        for k, v in pdict.iteritems():
            pdict_compat[self._FScompatibleNames(k)] = v        
        idict_compat = {}
        for k, v in idict.iteritems():
            idict_compat[self._FScompatibleNames(k)] = v
        return Point({'coorddict': dict(zip(self._FScompatibleNamesInv(gen.funcspec.auxvars),
                                gen.AuxVars(t, xdict_compat, pdict_compat,
                                            idict_compat))),
                      'coordtype': Float,
                      'norm': self._normord})


# -------------------------------------------------------------------------
#### Private functions

def getAuxVars(gen, t, icdict, pardict, keep_t=False):
    # prepare auxiliary variables eval'd at i.c.
    icdict_local = copy.copy(icdict)
    # to restore these to their former values after domain tests are done,
    # (in case events refer to auxiliary function initcond()
    old_gen_icdict = copy.copy(gen.initialconditions)
    # update icdict with auxiliary variable vals eval'd at t, for
    # use in event function calls for domain test.
    gen.initialconditions.update(icdict_local)
    icdict_local['t'] = t
    try:
        icdict_local.update(dict(zip(gen.funcspec.auxvars,
                    apply(getattr(gen, gen.funcspec.auxspec[1]),
                                      (t, sortedDictValues(icdict_local),
                                           sortedDictValues(pardict))) )))
    except AttributeError:
        # no auxiliary variables for this generator
        pass
    # restore generator's original initial conditions
    gen.initialconditions = old_gen_icdict
    if not keep_t:
        del icdict_local['t']
    return icdict_local


def findInitialGenerator(genInfo, t, xdict, pardict, intvars,
                         verboselevel=0):
    """Find eligible Generator to begin computation of trajectory."""
    eligibleGen = []
    if len(genInfo) == 1:
        return genInfo.values()[0]
    outcome = {}
    for infopair in genInfo.values():
        gen = infopair[0]
        outcome[gen.name] = {}
        t_test = gen.contains(gen.indepvariable.indepdomain, t, gen.checklevel)
        allvars = gen.funcspec.vars+gen.funcspec.auxvars
        # if generator has inputs then include initial value of those in pardict
        if gen.funcspec.inputs != []:
            try:
                input_ics = [gen.inputs[in_name](t) for in_name in gen.funcspec.inputs]
            except:
                raise RuntimeError, \
                      "Error in evaluating inputs for generator at initial condition"
            pardict.update(dict(zip(gen.funcspec.inputs, input_ics)))
        icdict = copy.copy(xdict)
        icdict.update(getAuxVars(gen, t, icdict, pardict, keep_t=True))
        # override icdict with any finite-valued generator
        # initial condition (deliberately preset in generator definition)
        for xname in allvars:
            if xname not in xdict or xname in intersect(intvars,xdict):
                if xname in xdict and xname in intvars:
                    # ignore
                    continue
                try:
                    if isfinite(gen.initialconditions[xname]):
                        icdict[xname] = gen.initialconditions[xname]
                except KeyError:
                    # ic (for internal variable) not defined for this
                    # generator, so ignore
                    pass
        x_test = True   # initial value
        for xname in icdict:
            if xname == 't':
                # 't' is a special value inserted into icdict, for use
                # in auxiliary variable evaluation at initial time, but not
                # a regular variable so don't do domain check here.
                continue
            newtest = gen.contains(\
                    gen.variables[xname].depdomain,
                    icdict[xname], 1)  # checklevel = 1 for this
            if verboselevel >=2:
                print "\nx_test for '%s' @ value %f:"%(xname, icdict[xname])
                print "depdomain is ", gen.variables[xname].depdomain.get()
                print "-> test result is ", newtest
            x_test = x_test and newtest
        if verboselevel >= 1:
            print "\nGenerator '%s' tests..."%gen.name
            print " ...for initial time: " + str(t_test)
            print " ...for initial state: " + str(x_test)
        outcome[gen.name]['t_test'] = t_test
        outcome[gen.name]['x_test'] = x_test
        # for all active terminal events in gen:
        # check that directionless zero-crossing events are not 0 at
        # this i.c. check that event conditions w/ directions aren't
        # invalidated at this i.c.
        dom_test = True     # initial value
        if t_test and x_test:
            teps = gen.eventstruct.query(['term', 'active'])
            for tep in teps:
                tep[1].initialconditions = icdict
                fval = tep[1]._fn(icdict, pardict)
                newtest = fval*tep[1].dircode <= gen._abseps
                if verboselevel==2:
                    print " ...domain test for event " \
                      "'%s' (fn val = %.8f): %s"%(tep[0],fval,str(newtest))
                # restore event initial value (not really necessary!)
#                tep[1].initialconditions = {}
                dom_test = dom_test and newtest
            if verboselevel >= 1:
                print " ...for domain inclusion: " + str(dom_test)
            outcome[gen.name]['dom_test'] = dom_test
            if dom_test:
                # this t, xdict is a valid initial state for this Generator
                eligibleGen.append((gen,infopair[1]))
    if eligibleGen == []:
        info(outcome, "\nOutcomes for eligibility tests, by generator")
        raise PyDSTool_ValueError, ('No eligible Generators from'
                           ' which to begin trajectory computation')
    if len(eligibleGen) > 1:
        info(outcome, "\nOutcomes for eligibility tests, by generator")
        raise PyDSTool_ValueError, ('Too many eligible Generators'
                           ' to start trajectory computation')
    # only remaining possibility is a single eligible Generator object
    return eligibleGen[0]

