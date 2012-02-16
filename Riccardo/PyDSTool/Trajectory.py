#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      23 Aug 2006
#$Revision: 1.1.6 $
"""Trajectory class.

   Robert Clewley, June 2005.
"""

# ----------------------------------------------------------------------------
# VERSION HISTORY
#
# Version 1.1.6, Aug 2006
#   Fixed sample method to properly use tlo and thi when precise=False.
#   Fixed underlyingMesh method to return correct shape of array when not
#     using datapoints attribute.
#
# Version 1.1.5, May 2006
#   Changed 'sampleTraj' method name to 'sample', keeping the old name as an
#    alias.
#   Added a "quick and dirty" findApproxPeriod function that calls scipy's
#    all-python newton solver.
#
# Version 1.1.4, Mar 2006
#   Updated sampleTraj to avoid trying to add event data to pointset in case
#    underlying mesh is being used (it will already contain the event data).
#
# Version 1.1.3, Feb 2006
#   Added use of _normord for Pointset norm order.
#
# Version 1.1.2, Jan 2006
#   Updated __call__ to reflect change to Pointset call always returning a
#     Point or a Pointset.
#   Added __deepcopy__ method.
#
# Version 1.1, Oct 2005
#   Changed sampleTraj to access underlying mesh when possible, using default
#     precise=False option. precise=True is still available for previous
#     functionality, but is a lot slower.
#   Simplified sampleTraj's mesh generation.
#   Fixed sampleTraj's default behaviour when dt is not provided, to use
#     underlying mesh where possible.
#   Added useGlobalTime=True option to __call__ and sampleTraj so that by
#     default, global time (defined by self.globalt0) is used for independent
#     variable referencing.
#
# Version 1.0.3, Sep 2005
#   Added support for hierarchical names.
#
# Version 1.0.2, July 2005
#   Added a trap for multiple events in sampleTraj method.
#   Added a trap for empty evtlist in sampleTraj.
#   Changed < to <= in sampleTraj:
#                elif t <= thi:
#                    tmesh.append(thi)
#     and added code to handle evtlist only having a single entry at thi
#
# Version 1.0.1, June 2005
#   Added 'parameterized' argument from ImplicitFnGen.
#   Added automatic variable names alphabetical ordering.
#   Added underlyingMesh method.
#   Changed underscore-separated method names to camel case.
#   Updated sampleTraj to actually use eventData in its output.
#   Added a trap for empty tmesh in sampleTraj.
# ----------------------------------------------------------------------------


# PyDSTool imports
from Variable import *
from Points import *
from utils import *
from common import *
from parseUtils import *
from errors import *

# Other imports
# temporary import of numarray until scipy upgrades to numarray3 (summer '05?),
from scipy import array as scipy_array
from scipy import sometrue, alltrue, any, all, shape
from scipy.optimize import newton
from numarray import array, arange, Float, Int, Complex, Float64, Int32, \
     Complex64, concatenate, zeros
import math
import copy
import sys

# -------------------------------------------------------------

## Exports

__all__ = ['Trajectory']


# -------------------------------------------------------------

class Trajectory(object):
    """Parameterized and non-parameterized trajectory class."""

    def __init__(self, name, varlist, globalt0=0, checklevel=0,
                 parameterized=True, norm=2):
        # parameterized is set False by ImplicitFnGen
        if isinstance(name, str):
            self.name = name
        else:
            raise TypeError, "name argument must be a string"
        self._normord = norm # order for norm of points within trajectory
        if isinstance(varlist, Variable):
            # then a singleton was passed
            varlist = [varlist]
        # order of listing in varnames provides order of coordinates for output
        self.varnames = [v.name for v in varlist]
        if not uniqueList(self.varnames):
            raise ValueError, "Variable names must be unique"
        self.varnames.sort()
        # non-parameterized trajectory?
        indepvarname = varlist[0].indepvarname
        if indepvarname in self.varnames or not parameterized:
            self._parameterized = False
        else:
            self._parameterized = True
        indepvartype = varlist[0].indepvartype
        coordtype = varlist[0].coordtype
        self.depdomain = {}
        if varlist[0].trajirange is None:
            indepdomain = varlist[0].indepdomain
            self.depdomain[varlist[0].coordname] = varlist[0].depdomain
        else:
            indepdomain = varlist[0].trajirange
            self.depdomain[varlist[0].coordname] = varlist[0].trajdrange
        for v in varlist[1:]:
            if isinstance(indepdomain, Array) or isinstance(indepdomain, NArray):
                if v.trajirange is None:
                    if all(v.indepdomain != indepdomain):
                        raise ValueError, ("Some variables in varlist argument have"
                                       " different independent variable domains")
                else:
                    if all(v.trajirange != indepdomain):
                        raise ValueError, ("Some variables in varlist argument have"
                                       " different independent variable domains")
            else:
                if v.trajirange is None:
                    if v.indepdomain != indepdomain:
                        raise ValueError, ("Some variables in varlist argument have"
                                       " different independent variable domains")
                else:
                    if v.trajirange != indepdomain:
                        raise ValueError, ("Some variables in varlist argument have"
                                       " different independent variable domains")
            if v.indepvarname != indepvarname:
                raise ValueError, ("Some variables in varlist "
                  "argument have different independent variable names")
            if v.coordtype != coordtype or v.indepvartype != indepvartype:
                raise ValueError, ("Some variables in varlist "
                  "argument have different coordinate or independent "
                  "variable types")
            if v.trajdrange is None:
                self.depdomain[v.coordname] = v.depdomain
            else:
                self.depdomain[v.coordname] = v.trajdrange
        self.indepvarname = indepvarname
        self.indepvartype = indepvartype
        # self.indepdomain is not checked when traj called, it's just for info
        self.indepdomain = indepdomain
        self.coordtype = coordtype
        self.variables = {}
        for v in varlist:
            assert isinstance(v, Variable)
            assert v.defined, ("Variables must be defined before building "
                               "a Trajectory object")
            self.variables[v.name] = v
        self.dimension = len(varlist)
        self._name_ix_map = invertMap(self.varnames)
        self.globalt0 = globalt0
        self.checklevel = checklevel


    def __call__(self, indepvar, coords=None, checklevel=None,
                 useGlobalTime=True):
        if checklevel is None:
            checklevel = self.checklevel
        if coords is None:
            coordlist = self.varnames
        elif isinstance(coords, int):
            coordlist = [self.varnames[coords]]
        elif isinstance(coords, str):
            coordlist = [coords]
        elif type(coords) in [list, tuple, Array, NArray]:
            if all([isinstance(c, str) for c in coords]):
                coordlist = [v for v in self.varnames if v in coords]
            elif any([isinstance(c, str) for c in coords]):
                raise TypeError, ("Cannot mix string and numeric values in "
                                  "coords")
            else:
                # numeric list assumed
                coordlist = [v for v in self.varnames if v in self._name_ix_map]
        else:
            raise TypeError, "Invalid type for coords argument"
        if type(indepvar) in [list, tuple, Array, NArray]:
            # ensure is an array so that we can subtract globalt0
            if useGlobalTime:
                indepvals = array(indepvar) - self.globalt0
            else:
                indepvals = indepvar
            return Pointset({'coordarray': [v(indepvals, checklevel) \
                                                  for v in \
                                 [self.variables[vn] for vn in coordlist]],
                             'coordnames': coordlist,
                             'coordtype': self.coordtype,
                             'norm': self._normord,
                             'name': self.name + "_sample"})
        else:
            if useGlobalTime:
                indepvar = indepvar - self.globalt0
            if len(coordlist) == 1:
                return self.variables[coordlist[0]](indepvar, checklevel)
            else:
                return Point({'coordarray': [v(indepvar, checklevel) \
                                               for v in \
                                 [self.variables[vn] for vn in coordlist]],
                          'coordnames': coordlist,
                          'coordtype': self.coordtype,
                          'norm': self._normord})


    def underlyingMesh(self, varlist=None):
        """Return a dictionary of the underlying independent variables` meshes,
        where they exist. (Dictionary contains references to the meshes,
        not copies of the meshes.)"""

        # contrast to code to pull out both dependent and independent
        # variables' data points in Variable class' getDataPoints
        if varlist is None:
            varlist = self.varnames
        if not isinstance(varlist, list):
            raise TypeError, 'varlist argument must be a list'
        meshdict = {}
        for varname in varlist:
            try:
                # works if .output is an interp1d instance
                meshdict[varname] = self.variables[varname].output.datapoints
            except AttributeError:
                try:
                    # works if .output is a VarCaller instance (with underlying Pointset)
                    pts = self.variables[varname].output.pts
                    meshdict[varname] = array([pts.indepvararray,
                                               pts.coordarray[0]])
                except AttributeError:
                    meshdict[varname] = None
        return meshdict


    def sample(self, varlist=None, dt=None, tlo=None, thi=None,
                   eventdata=None, precise=False, useGlobalTime=True):
        """Uniformly sample the named trajectory over range indicated.
        If optional 'eventdata' argument included, the event points are included
        in the mesh (e.g. for plotting purposes).

        The order of variable names in varlist is ignored.

        precise=True causes the trajectory position to be evaluated precisely
        at the t values specified, which will invoke slow interpolation.
        precise=False (default) causes the nearest underlying mesh positions
        to be used, if available (otherwise the behaviour is the same as
                                  precise=True)

        If dt is not given, the underlying time mesh is used, if available.
        """

        if varlist is None:
            varlist_sorted = self.varnames
        else:
            assert isinstance(varlist, list), 'varlist argument must be a list'
            varlist_sorted = copy.copy(varlist)
            varlist_sorted.sort()
        tlo_base = self.indepdomain.get(0)
        thi_base = self.indepdomain.get(1)
        if tlo is None:
            tlo = tlo_base
        elif useGlobalTime:
            tlo = tlo - self.globalt0
            if tlo < tlo_base:
                tlo = tlo_base
            if tlo >= thi_base:
                raise ValueError, "tlo too large"
        elif tlo < tlo_base:
            tlo = tlo_base
        elif tlo >= thi_base:
            raise ValueError, "tlo too large"
        if thi is None:
            thi = thi_base
        elif useGlobalTime:
            thi = thi - self.globalt0
            if thi > thi_base:
                thi = thi_base
            if thi <= tlo_base:
                raise ValueError, "thi too small"
        elif thi > thi_base:
            thi = thi_base
        elif thi <= tlo_base:
            raise ValueError, "thi too small"
        assert tlo < thi, 't start point must be less than t endpoint'
        if dt is not None and dt >= abs(thi-tlo):
            if precise:
                print "dt = %f for interval [%f,%f]"%(dt,tlo,thi)
                raise ValueError('dt must be smaller than time interval')
            else:
                dt = (thi-tlo)/10.
##        do_events = eventdata is not None
        if not precise or dt is None:
            # attempt to use underlying mesh for each variable,
            # if available
            meshdict = self.underlyingMesh()
            meshes_ok = False
            if meshdict is not None:
                if not any([meshdict[v] is None for v in varlist_sorted]):
                    # ensure all vars' meshes are identical!
                    firstvar_tmesh = meshdict[varlist_sorted[0]][0]
                    if not type(firstvar_tmesh) in [Array, NArray]:
                        firstvar_tmesh = array(firstvar_tmesh)
                    meshes_ok = True
                    # explicitly checking that the arrays are the same
                    # is very slow -- just catch a resulting error later
#                    try:
#                        for v in varlist_sorted[1:]:
#                            assert all(meshdict[v][0] == firstvar_tmesh)
#                        meshes_ok = True
#                    except AssertionError:
#                        meshes_ok = False
            if meshes_ok:
                loix_a = firstvar_tmesh>=tlo
                loix_a = loix_a.tolist()
                try:
                    loix = loix_a.index(1)
                except ValueError:
                    loix = 0
                hiix_a = firstvar_tmesh>thi
                hiix_a = hiix_a.tolist()
                try:
                    hiix = hiix_a.index(1)
                except ValueError:
                    hiix = len(firstvar_tmesh)
                tmesh = firstvar_tmesh[loix:hiix]
                #do_events = False
        if dt is None:
            if not meshes_ok:
                raise ValueError("Underlying mesh of trajectory is not"
                                    " the same for all variables: dt must "
                                    "be specified in this case")
        else:
            tmesh = concatenate([arange(tlo, thi, dt),[thi]])
            loix = 0
            hiix = len(tmesh)
        if eventdata: #and do_events:
            if isinstance(eventdata, dict):
                assert type(eventdata.values()[0]) in [Array, NArray, list], \
                   "eventdata must be list of times or dictionary of lists of times"
                evtlist = orderEventData(eventdata, nonames=True)
            elif type(eventdata) in [Array, NArray, list]:
                evtlist = eventdata
            else:
                raise TypeError, \
                   "eventdata must be list of times or dictionary of lists of times"
            if evtlist != []:
                tmesh, ins_ixs = insertInOrder(tmesh, evtlist, True)
            else:
                ins_ixs = []
        if len(tmesh) > 0:
            if dt is None:
                # meshes_ok already checked
##                coorddict = {}
##                for v in varlist_sorted:
##                    coorddict[v] = meshdict[v][1]
                coorddict = dict(zip(varlist_sorted,
                                 [m[1][loix:hiix] for m in sortedDictValues(meshdict)]))
                if eventdata: #and do_events:
                    # insert var values at events (SLOW!)
                    # can only insert to lists, so have to convert coorddict arrays
                    # to lists first.
                    for tpos in xrange(len(ins_ixs)):
                        tix = ins_ixs[tpos]
                        t = evtlist[tpos]
                        x = self(t, varlist_sorted)
                        for v in varlist_sorted:
                            try:
                                coorddict[v].insert(tix, x[v])
                            except AttributeError:
                                coorddict[v] = coorddict[v].tolist()
                                coorddict[v].insert(tix, x[v])
                if useGlobalTime and self.globalt0 != 0:
                    tmesh = tmesh + self.globalt0
                return Pointset({'coorddict': coorddict,
                           'coordtype': Float,
                           'coordnames': varlist_sorted,
                           'indepvararray': tmesh,
                           'indepvarname': self.indepvarname,
                           'indepvartype': Float,
                           'norm': self._normord,
                           'name': self.name + "_sample"})
            if not precise:
                # get closest mesh indices corresponding to tmesh
                # times for first var
                try:
                    try:
                        closest = findClosestArray(firstvar_tmesh,tmesh,dt/2)
                        find_ixs = True
                    except UnboundLocalError:
                        # meshes_ok==False and firstvar_tmesh never created
                        # so don't need to find indexes of closest meshpoints,
                        # just use original mesh
                        find_ixs = False
                    coorddict = {}
                    if find_ixs:
                        closest_unique = makeSeqUnique(closest, True, True)
                        if useGlobalTime and self.globalt0 != 0:
                            tmesh = firstvar_tmesh[closest_unique] + self.globalt0
                        else:
                            tmesh = firstvar_tmesh[closest_unique]
                        for v in varlist_sorted:
                            dat = self.variables[v].output.datapoints[1]
                            try:
                                coorddict[v] = dat[closest_unique]
                            except TypeError:
                                dat = array(dat)
                                coorddict[v] = dat[closest_unique]
                    else:
                        for v in varlist_sorted:
                            coorddict[v] = self.variables[v].output.datapoints[1]
                    return Pointset({'coorddict': coorddict,
                                   'coordtype': Float,
                                   'coordnames': varlist_sorted,
                                   'indepvararray': tmesh,
                                   'indepvarname': self.indepvarname,
                                   'indepvartype': Float,
                                   'norm': self._normord,
                                   'name': self.name + "_sample"})
                except (AttributeError, IndexError):
                    # meshes of variables did not match up, so
                    # continue to 'else' case, as if precise=True
                    pass
            # else mesh not available for some variables so we have
            # a mixed-source trajectory, and simple solution is to
            # treat as if precise==True
            # Cases: either precise == True or could not extract
            # underlying meshes for some variables, or the meshes were not
            # identical.
            pset = Pointset({'coordarray': self(tmesh, varlist_sorted,
                                                useGlobalTime=False).toarray(),
                     'coordnames': varlist_sorted,
                     'indepvararray': array(tmesh),
                     'indepvarname': self.indepvarname,
                     'name': self.name + "_sample"})
            # if varlist was not sorted, call to self could potentially return
            # an array in a different order to that of varlist.
            if useGlobalTime and self.globalt0 != 0:
                pset.indepvararray += self.globalt0
                pset.makeIxMaps()
            return pset
        else:
            return None


    sampleTraj = sample


    def findApproxPeriod(self, t0, t1_guess=None, T_guess=None, v=None,
                   ttol=1e-5, rtol=1e-2):
        """findPeriod does not ensure minimum period (it is up to the user to
        select a good initial guess."""
        if t1_guess is None and T_guess is None:
            raise ValueError, "Must supply guess for either t1 or period, T"
        if t1_guess is None:
            if T_guess <= 0:
                raise ValueError, "Guess for period T must be positive"
            t1_guess = t0+T_guess
        if T_guess is None:
            T_guess = t0+t1_guess
        if t1_guess <= t0:
            raise ValueError, "t1 guess must be greater than t0"
        if v is None:
            v = self.varnames[0]
        p0 = self.__call__(t0)
        p0v = p0[v]
        def f(t):
            try:
                return self.__call__(t)[v]-p0v
            except PyDSTool_BoundsError:
                return 10000.*(t-t0)
            except:
                print "Error at t=%f: "%t, sys.exc_info()[0], \
                            sys.exc_info()[1]
                raise
        result = newton(f, t1_guess, tol=ttol)
        val = self.__call__(result)-p0
        rval = [abs(val[vn]/p0[vn]) for vn in self.varnames]
        if max(rval)<rtol:
            return abs(result-t0)
        else:
            print "Did not converge. The endpoint difference was:\n", \
              repr(rval), \
              "\nwith infinity-norm %f > %f tolerance.\n"%(max(rval),rtol), \
              "Try a different starting point,", \
              "a different test variable, or reduce relative tolerance."
        

    def __copy__(self):
        pickledself = pickle.dumps(self)
        c = pickle.loads(pickledself)
##        c.indepdomain = copy.copy(self.indepdomain)
##        c.depdomain = copy.copy(self.depdomain)
##        for vname, v in self.variables.iteritems():
##            c.variables[vname] = copy.copy(v)
        return c


    def __deepcopy__(self, memo=None, _nil=[]):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


    def __getstate__(self):
        d = copy.copy(self.__dict__)
        # remove reference to Cfunc types by converting them to strings
        d['indepvartype'] = _pytypename[self.indepvartype]
        d['coordtype'] = _pytypename[self.coordtype]
        return d


    def __setstate__(self, state):
        self.__dict__.update(state)
        # reinstate Cfunc types
        self.indepvartype = _pynametype[self.indepvartype]
        self.coordtype = _pynametype[self.coordtype]


    def __del__(self):
        # necessary for __del__ methods of variables to be accessed
        for v in self.variables.values():
            v.__del__()


    def __repr__(self):
        return self._infostr(verbose=0)


    __str__ = __repr__


    def _infostr(self, verbose=1):
        if verbose == 0:
            outputStr = "Trajectory " + self.name
        elif verbose >= 1:
            outputStr = "Trajectory " + self.name + "\n  of variables: " + \
                    str(self.variables.keys())
            outputStr += "\n  over domains: " + str(self.depdomain)
            if verbose == 2:
                outputStr += joinStrs(["\n"+v.info(1) for v in \
                                            self.variables.values()])
        return outputStr


    def info(self, verboselevel=1):
        print self._infostr(verboselevel)


# ---------------------------------------------------------------------


if __name__ == '__main__':
    v1 = Variable(Pointset({'coordarray': array(range(10), Float)*0.1,
                            'indepvararray': array(range(10), Float)*0.5
                            }), name='v1')
    v2 = Variable(Pointset({'coordarray': array(range(10), Float)*0.25+1.0,
                            'indepvararray': array(range(10), Float)*0.5
                            }), name='v2')
    traj = Trajectory('test1', [v1,v2])
    print "print traj(0.5, checklevel=2) => ", traj(0.5, checklevel=2)
    print "print traj([0., 0.5]) => ", traj([0., 0.5], 'v1')
    print "traj(0.4, 0, checklevel=2) =>"
    try:
        traj(0.4, 0, checklevel=2)
    except ValueError, e:
        print " ... raised error: ", e
    print "Tests passed."
