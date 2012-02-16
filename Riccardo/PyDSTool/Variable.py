#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      27 January 2006
#$Revision: 2.0.5 $
"""Variable is a one-dimensional discrete and continuous real variable class.

   Robert Clewley, July 2005
"""

# ----------------------------------------------------------------------------
# VERSION HISTORY
#
# Version 2.0.5, January 2006
#   Updated __call__ and setOutput to reflect change to Pointset call always
#     returning a Point or a Pointset. Now uses Varcaller wrapper for
#     output consistent with other underlying types of 'outputdata'.
#   Added __deepcopy__ method.
#
# Version 2.0.4, October 2005
#   Minor bug fixes and additional exception traceback information.
#
# Version 2.0.3, July 2005
#   Changed __del__ to not delete 'output' method -- it was causing
#     IF_squarespike_model and pest_test2 tests to fail
#     Moved indepvar_ok = True in __call__, case checklevel == 1 otherwise
#       only checklevel == 2 uncertain cases were getting through.
#     Removed support for tuples in vectorized calls.
#     Fixed __call__ for checklevels 1 & 2 to catch multiple uncertainties
#       in indepvar and depvar, and also to adjust out-of-bounds but
#       uncertain indepvars that probably occurred due to rounding error.
#     showWarnings method added
#
# Version 2.0.2, June 2005
#   Added support for ImplicitFnGen and saving of objects using Pickle
#   Using a fixed version of Python pickle to get around Inf and Nan load
#     bug
#   Reverted to pickle version of copy method, because self._refvars was
#     not copying properly
#
# Version 2.0.1, June 2005
#   Made method names more consistent with use of camel case
#
# ----------------------------------------------------------------------------


# PyDSTool imports
from utils import *
from common import *
from errors import *
from Points import *
from Interval import *
from FuncSpec import ImpFuncSpec

# Other imports temporary import of numarray until scipy upgrades
from scipy import Inf, NaN, isfinite, sometrue, alltrue, any, all
from scipy import array as scipy_array
from numarray import array, Float, Int, Complex, Float64, Int32, Complex64

import copy
import types, math, random

__all__ = ['Variable', 'OutputFn', 'isinputcts', 'isinputdiscrete',
           'isoutputcts', 'isoutputdiscrete',
           'iscontinuous', 'isdiscrete']


# ------------------------------------------------------------------


class Variable(object):
    """One-dimensional discrete and continuous real variable class.
    """

    def __init__(self, outputdata=None, indepdomain=None, depdomain=None,
                 name='noname'):
        # funcreg stores function data for dynamically created methods
        # to allow a Variable to be copied using pickling
        self._funcreg = {}
        if isinstance(name, str):
            self.name = name
        else:
            raise TypeError, "name argument must be a string"
        # defaults for empty 'placeholder' Variables
        if outputdata is None or isinstance(outputdata, Pointset):
            if indepdomain is None:
                indepdomain = 't'
            if depdomain is None:
                depdomain = 'x'
        # set some initial values so that can test if what changed
        # after calling setOutput()
        self._vectorizable = True
        self.defined = False
        self.indepdomain = None
        self.indepvartype = None
        self.indepvarname = None
        self.depdomain = None
        self.coordtype = None
        self.coordname = None
        self._refvars = None   # for use with ExplicitFnGen
        # Ranges covered by the current trajectory held (if known)
        self.trajirange = None
        self.trajdrange = None
        # independent variable domain
        self.setIndepdomain(indepdomain)
        # output function
        self.setOutput(outputdata)
        # dependent variable domain
        self.setDepdomain(depdomain)
        assert self.coordname != self.indepvarname, ("Independent variable "
                                "name and coordinate name must be different")
        self.warnings = []


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


    def addMethods(self, funcspec):
        """Add dynamically-created methods to Veriable object"""

        # Add the auxiliary function specs to this Variable's namespace
        for auxfnname in funcspec.auxfns:
            fninfo = funcspec.auxfns[auxfnname]
            if not hasattr(Variable, fninfo[1]):
                # user-defined auxiliary functions
                # (built-ins are provided explicitly)
                try:
                    exec fninfo[0] in globals()
                except:
                    print 'Error in supplied auxiliary function code'
                    raise
                self._funcreg[fninfo[1]] = ('Variable', fninfo[0])
                setattr(Variable, fninfo[1], eval(fninfo[1]))
        # Add the spec function to this Variable's namespace
        fninfo_spec = funcspec.spec
        if not hasattr(Variable, fninfo_spec[1]):
            try:
                exec fninfo_spec[0] in globals()
            except:
                print 'Error in supplied functional specification code'
                raise
            self._funcreg[fninfo_spec[1]] = ('Variable', fninfo_spec[0])
            setattr(Variable, fninfo_spec[1], eval(fninfo_spec[1]))
        # Add the auxiliary spec function (if present) to this Var's namespace
        if funcspec.auxspec:
            fninfo_auxspec = funcspec.auxspec
            if not hasattr(Variable, fninfo_auxspec[1]):
                try:
                    exec fninfo_auxspec[0] in globals()
                except:
                    print 'Error in supplied auxiliary variable code'
                    raise
                self._funcreg[fninfo_auxspec[1]] = ('Variable', fninfo_auxspec[0])
                setattr(Variable, fninfo_auxspec[1], eval(fninfo_auxspec[1]))
        # For implicit functions
        if isinstance(funcspec, ImpFuncSpec):
            impfn_name = funcspec.algparams['impfn_name']
            if funcspec.algparams['jac']:
                jac_str = "fprime=funcspec.algparams['jac'],"
            else:
                jac_str = ""
            # Wrap spec fn like this as it has been set up as a
            # method, but want to call as regular function
            # *** ALSO *** spec fn has signature (ds, t, x, p)
            # but implicit function solvers expect
            # (x, t, p), so have to switch 1st and 2nd args here
            # after 'ds' filled with None
            if len(funcspec.vars) == 1:
                # dimension == 1, so have to make list output from spec
                # into a scalar
                # Also, scalar a1 needs to be put into list form for
                # acceptance as x in spec fn
                specfn_str = "lambda a1, a2, a3: " \
                  + fninfo_spec[1] \
                  + "(None, a2, [a1], a3)[0]"
            else:
                # for dimension > 1 a1 will already be an array / list
                specfn_str = "lambda a1, a2, a3: " \
                  + fninfo_spec[1] \
                  + "(None, a2, a1, a3)"
            this_scope = globals()   # WE CHANGE GLOBALS()
            this_scope.update({'funcspec': locals()['funcspec'],
                               'fninfo_spec': locals()['fninfo_spec']})
            impfn_str = impfn_name + \
                " = makeImplicitFunc(" + specfn_str + "," \
                + jac_str + """x0=funcspec.algparams['x0'],
                               extrafargs=(funcspec.algparams['pars'],),
                               xtolval=funcspec.algparams['atol'],
                               maxnumiter=funcspec.algparams['maxnumiter'],
                               solmethod=funcspec.algparams['solvemethod'],
                               standalone=False)"""
            try:
                exec impfn_str in this_scope
            except:
                print 'Error in supplied implicit function code'
                raise
            # record special reference to the implicit fn,
            # as its a method of Variable (for delete method).
            self._funcreg['_impfn'] = (impfn_name, impfn_str)
            # In previous versions setattr was to self, not the Variable class
            setattr(Variable, impfn_name, eval(impfn_name))
            # clean up globals() afterwards
            del this_scope['funcspec']
            del this_scope['fninfo_spec']


    def setOutput(self, outputdata, funcspec=None,
                  globalt0=0, var_namemap=None, ics=None, refvars=None):
        """Dynamically create 'output' method of Variable"""

        self.globalt0 = globalt0
        if type(outputdata) in [types.FunctionType,
                                types.BuiltinFunctionType,
                                types.MethodType]:
            # Variable generated from function, given in closed form
            self.output = outputdata
            assert ics is None, "Invalid option for this type of output"
            if outputdata != noneFn:
                self.defined = True
        elif isinstance(outputdata, tuple):
            # For ExplicitFnGen or ImplicitFnGen types, whose functional forms
            # may need to access these at call time.
            assert len(outputdata) == 2, "Incorrect size of outputdata tuple"
            if funcspec is not None:
                self.addMethods(funcspec)
                self._var_namemap = var_namemap
                self._funcreg['funcspec'] = (None, funcspec)
            else:
                raise ValueError, 'funcspec missing in setOutput'
            # Add the specific mapping functions for Ex/ImplicitFnGen objects
            try:
                exec outputdata[1] in globals()
            except:
                print 'Internal Error in _mapspecfn code'
                raise
            has_op = hasattr(self, 'output')
            # have to define this function in here because use of lambda
            # won't allow me to pickle the Variable object
            if not has_op or (has_op and self.output is noneFn):
                def wrap_output(arg):
                    return eval(outputdata[0])(self, arg)
                setattr(self, 'output', wrap_output)
            self._funcreg['outputdata'] = (None, outputdata)
            t0 = self.indepdomain.get(0)
            if ics is None and not isinstance(funcspec, ImpFuncSpec):
                try:
                    self.initialconditions = {self.coordname: self.output(t0)}
                except ValueError:
                    self.initialconditions = {self.coordname: NaN}
                except TypeError:
                    print "Debugging info: self.output = ", self.output
                    raise
            else:
                self.initialconditions = ics
            self._vectorizable = False
            self._refvars = refvars
            self.defined = True
        elif type(outputdata) in [types.InstanceType, types.TypeType, \
                                 OutputFn] or isinstance(outputdata, interpclass):
            # Variable generated by callable object that generates values over
            # mesh points that it holds, e.g. by interpolation
            # (InstanceType and TypeType are for backwards compatibility, e.g.
            # for old SciPy interpolate code that uses Classic Classes)
            assert ics is None, "Invalid option for this type of output"
            assert '__call__' in dir(outputdata), "Must provide callable object"
            self.output = outputdata
            if hasattr(outputdata, 'datapoints'):
                if hasattr(outputdata, 'types'):
                    deptypestr = _pytypename[outputdata.types[0]]
                    indeptypestr = _pytypename[outputdata.types[1]]
                else:
                    # default
                    deptypestr = indeptypestr = 'float'
                if isinstance(outputdata.datapoints[0], Interval):
                    assert outputdata.types[0] == \
                           outputdata.datapoints[0].type, \
                           "Inconsistent type with Interval bounds"
                    self.trajirange = outputdata.datapoints[0]
                else:
                    self.trajirange = Interval('traj_indep_bd', indeptypestr,
                                              extent(outputdata.datapoints[0]))
                if isinstance(outputdata.datapoints[1], Interval):
                    assert outputdata.types[1] == \
                           outputdata.datapoints[1].type, \
                           "Inconsistent type with Interval bounds"
                    self.trajdrange = outputdata.datapoints[1]
                else:
                    self.trajdrange = Interval('traj_dep_bd', deptypestr,
                                              extent(outputdata.datapoints[1]))
            self.defined = True
        elif isinstance(outputdata, Pointset):
            # Variable generated from a pointset (without interpolation)
            assert ics is None, "Invalid option for this type of output"
            assert isparameterized(outputdata), ("Must only pass parameterized"
                                                 " pointsets")
            if outputdata.dimension == 1:
                self.coordname = copy.copy(outputdata.coordnames[0])
                self.indepvarname = outputdata.indepvarname
                self.output = Varcaller(outputdata)
                self.coordtype = outputdata.coordtype
                self.indepvartype = outputdata.indepvartype
                if self.indepdomain is not None:
                    for v in outputdata[self.indepvarname]:
                        if not v in self.indepdomain:
                            raise ValueError, ("New Pointset data violates "
                               "independent variable domain already specified")
                if self.depdomain is not None:
                    for v in outputdata[self.coordname]:
                        if not v in self.depdomain:
                            raise ValueError, ("New Pointset data violates "
                               "dependent variable domain already specified")
                self.trajirange = Interval('traj_indep_bd',
                                            _pytypename[self.indepvartype],
                                            extent(outputdata.indepvararray))
                self.trajdrange = Interval('traj_dep_bd',
                                           _pytypename[self.coordtype],
                                           extent(outputdata[self.coordname]))
                self.defined = True
            else:
                raise ValueError, \
                      "Pointset data must be 1D to create a Variable"
        elif outputdata is None:
            # placeholder for an unknown output type
            assert ics is None, "Invalid option when outputdata argument is None"
            self.output = noneFn
            self.defined = False
        else:
            raise TypeError, ("Invalid type for data argument: " \
                              +str(type(outputdata)))


    def setIndepdomain(self, indepdomain):
        if isinstance(indepdomain, str):
            self.indepvarname = indepdomain
            if self.indepdomain is not None:
                # If indepdomain already set and indepvarname is none then
                # name won't get put in place unless we force it here
                self.indepvarname = indepdomain
                self.indepdomain.name = indepdomain
            else:
                self.indepdomain = Interval(self.indepvarname, 'float',
                                           [-Inf, Inf])
            self.indepvartype = Float
        else:
            if isinstance(indepdomain, Interval):
                if self.trajirange:
                    if indepdomain.contains(self.trajirange) is notcontained:
                           raise ValueError, ("Cannot set independent variable"
                                          " domain inside current trajectory's"
                                          " range")
                self.indepdomain = indepdomain
                self.indepvarname = indepdomain.name
                self.indepvartype = _pynametype[indepdomain.typestr]
            elif isinstance(indepdomain, dict):
                assert len(indepdomain) == 1, "Dictionary must have only 1 entry"
                d = indepdomain.values()[0]
                assert all(isfinite(d)), "Values must be finite"
                assert isincreasing(d), "Values must be increasing"
                if self.trajirange:
                    assert self.trajirange.get(0) in d
                    assert self.trajirange.get(1) in d
                self.indepvarname = indepdomain.keys()[0]
                if isinstance(d, list):
                    self.indepdomain = array(d)
                elif isinstance(d, Array) or isinstance(d, NArray):
                    self.indepdomain = d
                else:
                    raise TypeError, ("Invalid type for independent "
                                      "variable domain")
                self.indepvartype = self.indepdomain.type()
            else:
                print "Independent variable argument domain was:", indepdomain
                raise TypeError, ("Invalid type for independent variable "
                                  "domain")            


    def setDepdomain(self, depdomain):
        if isinstance(depdomain, str):
            self.coordname = depdomain
            if self.depdomain is None:
                if self.coordtype is None:
                    self.depdomain = Interval(self.coordname, 'float',
                                                 [-Inf, Inf])
                    self.coordtype = Float
                else:
                    self.depdomain = Interval(self.coordname,
                                                 _pytypename[self.coordtype],
                                                 _nummaxmin[self.coordtype])
            else:
                # If interp functions supplied then don't have a name for
                # Interval yet, so update it.
                if isinstance(self.output, interpclass) and \
                   isinstance(self.depdomain, Interval):
                    self.depdomain.name = depdomain
                else:
                    assert isinstance(self.output, Pointset)
                    print "Warning: Dependent variable already named. " + \
                          "Ignoring user-supplied name."
        else:
            if isinstance(depdomain, Interval):
                if self.trajdrange:
                    if depdomain.contains(self.trajdrange) is notcontained:
                        raise ValueError, ("Cannot set dependent variable "
                                          "domain inside current trajectory's "
                                          "range")
                self.depdomain = depdomain
                self.coordname = depdomain.name
                if self.coordtype is None:
                    self.coordtype = depdomain.type
                elif self.coordtype == depdomain.type:
                    pass
                else:
                    raise TypeError, ("Mismatch between type of depdomain "
                                      "argument and Pointset coord data")
            elif isinstance(depdomain, dict):
                assert len(depdomain) == 1, "Dictionary must have only 1 entry"
                d = depdomain.values()[0]
                if self.trajdrange:
                    assert self.trajdrange.get(0) in d
                    assert self.trajdrange.get(1) in d
##                if not isincreasing(d):
##                    d.sort()
                assert all(isfinite(d)), "Values must be finite"
                self.coordname = depdomain.keys()[0]
                if isinstance(d, list):
                    if self.coordtype is not None:
                        self.depdomain = array(d, self.coordtype)
                    else:
                        self.depdomain = array(d)
                elif isinstance(d, Array) or isinstance(d, NArray):
                    da = array(d)
                    if self.coordtype is not None and \
                       self.coordtype != da.type():
                        raise TypeError, ("Mismatch between type of depdomain "
                                          "argument and Pointset coord data")
                    else:
                        self.depdomain = da
                else:
                    raise TypeError, ("Invalid type for dependent variable "
                                      "domain")
                self.coordtype = self.depdomain.type()
            else:
                print "Dependent variable domain argument was:", depdomain
                raise TypeError, "Invalid type for dependent variable domain"
            if isinstance(self.output, Pointset):
                assert self.coordname == self.output.coordnames[0], \
                       "Mismatch between Pointset coord name and declared name"
                assert self.indepvarname == self.output.indepvarname, \
                       ("Mismatch between Pointset independent variable name "
                        "and declared name")


    def __call__(self, indepvar, checklevel=0):
        if checklevel == 0:
            # level 0 -- no depvar bounds checking at all
            # (no need to check for indepvar as list case, which output
            # should know how to handle)
            try:
                if not self._vectorizable and type(indepvar) in \
                   [list, Array, NArray]:
                    return [self.output(ival) for ival in indepvar]
                else:
                    return self.output(indepvar)
            except OverflowError:
                print "Overflow error at independent variable =", indepvar
                raise
        elif checklevel in [1,2]:
            if self.trajirange is None:
                idep = self.indepdomain
            else:
                # use known bounds on indep variable imposed by self.output
                idep = self.trajirange
            indepvar_ok = True
            # level 1 -- ignore uncertain cases (treat as contained)
            # level 2 -- warn on uncertain (treat as contained)
            if type(indepvar) in [list, Array, NArray]:
                vectorizable = self._vectorizable
                for d in indepvar:
                    # use 'in' so that this is compatible with
                    # interval, array and index indeps
                    try:
                        contresult = d in idep
                    except PyDSTool_UncertainValueError:
                        contresult = True
                        # adjust for rounding error so that interpolator
                        # does not barf on out-of-range values
                        if d < idep.get(0):
                            try:
                                # list
                                dix = indepvar.index(d)
                            except AttributeError:
                                # array
                                dix = indepvar.tolist().index()
                            indepvar[dix] = idep.get(0)
                        elif d > idep.get(1):
                            try:
                                # list
                                dix = indepvar.index(d)
                            except AttributeError:
                                # array
                                dix = indepvar.tolist().index()
                            indepvar[dix] = idep.get(1)
                        if checklevel == 2:
                            self.warnings.append((d, None))
                    if not contresult:
                        indepvar_ok = False
                        break
            else:
                vectorizable = True
                try:
                    indepvar_ok = indepvar in idep
                except PyDSTool_UncertainValueError, errinfo:
                    # adjust for rounding error so that interpolator
                    # does not barf on out-of-range values
                    if indepvar < idep.get(0):
                        indepvar = idep.get(0)
                    elif indepvar > idep.get(1):
                        indepvar = idep.get(1)
                    if checklevel == 2:
                        self.warnings.append((indepvar, None))
            # continue to get dependent variable value, unless indep
            # value was not OK
            if not indepvar_ok:
##                print "*** Debug info for variable: ", self.name
##                print "Interval rounding tolerance was", idep._abseps
                if vectorizable:
                    raise ValueError, ('Independent variable value(s) '
                                       'out of range in Variable call')
                else:
                    raise ValueError, ('Independent variable value '+\
                                   str(indepvar) + ' out of '
                                   'range in Variable call')
            try:
                if vectorizable:
                    depvar = self.output(indepvar)
                else:
                    depvar = [self.output(ival) for ival in indepvar]
                depvar_ok = True
            except PyDSTool_BoundsError, errinfo:
                depvar_ok = False
            # Now check that all computed values were in depdomain
            if depvar_ok:
                # no need to use self.trajdrange instead of
                # self.depdomain because we trust that self.output
                # generated the output within its own bounds!
                if type(depvar) in [list, Array, NArray, Pointset]:
                    if isinstance(depvar, Pointset):
                        dv = depvar.toarray()
                    else:
                        dv = depvar
                    for d in dv:
                        # use 'in' so that this is compatible with
                        # interval, array and index indeps
                        try:
                            contresult = d in self.depdomain
                        except PyDSTool_UncertainValueError, errinfo:
                            contresult = True
                            if checklevel == 2:
                                # find which indepvar was the cause of
                                # the uncertain value
                                try:
                                    # list
                                    depix = dv.index(d)
                                except AttributeError:
                                    # array
                                    depix = dv.tolist().index(d)
                                self.warnings.append((indepvar[depix], errinfo.val))
                        if not isfinite(d):
                            print "Warning: Return value for independent " + \
                                   "variable %f was not finite"%indepvar + \
                                   "for Variable '%s'"%self.name
                        if not contresult:
                            depvar_ok = False
                            break
                elif depvar is None:
                    print "*** Debug info for variable: ", self.name
                    print "Independent variable domain: ", self.indepdomain._infostr(1)
                    print "Dependent variable domain: ", self.depdomain._infostr(1)
                    raise ValueError, ("Cannot compute a return value for "
                                          "this independent variable value: "
                                          + str(indepvar))
                else:
                    if isinstance(depvar, Point):
                        dv = depvar[0]
                    else:
                        dv = depvar
                    try:
                        depvar_ok = dv in self.depdomain
                    except PyDSTool_UncertainValueError, errinfo:
                        if checklevel == 2:
                            self.warnings.append((indepvar, errinfo.varval))
                    if not isfinite(dv):
                        print "Warning: Return value for independent variable", indepvar
                        print "was not finite for Variable '%s'"%self.name
            # return value if depvar in bounds
            if depvar_ok:
                return dv
            else:
                print "Variable '%s' -"%self.name, "dependent var domain: ", \
                      self.depdomain._infostr(1)
                self.showWarnings()
                if vectorizable:
                    raise ValueError, ('Computed value(s) outside'
                                   ' validity range in Variable call')                    
                else:
                    raise ValueError, ('Computed value outside'
                                   ' validity range in Variable call')
        else:
            # level 3 -- exception will be raised for uncertain case
            indepvar_ok = False
            try:
                # don't trap uncertain case exception from
                # Interval.__contains__
                if type(indepvar) in [list, Array, NArray]:
                    vectorizable = self._vectorizable
                    indepvar_ok = all([i in self.indepdomain for i in \
                                           indepvar])
                else:
                    vectorizable = True
                    indepvar_ok = indepvar in self.indepdomain
                if not indepvar_ok:
                    raise ValueError, ('Independent variable '+str(indepvar)+\
                                       ' out of range in Variable call')
            except TypeError, e:
                raise TypeError, ('Something messed up with the Variable '
                                  'initialization: '+str(e))
            # Don't need 'if indepvar_ok' because exception would have
            # been raised.
            # For this checklevel, don't trap uncertain case exception from
            # Interval.__contains__
            try:
                if vectorizable:
                    depvar = self.output(indepvar)
                else:
                    depvar = [self.output(ival) for ival in indepvar]
                if type(depvar) in [list, Array, NArray]:
                    depvar_ok = all([d in self.depdomain for d in \
                                         depvar])
                    if depvar_ok:
                        return depvar
                    else:
                        if vectorizable:
                            raise ValueError, ('Computed value(s) '
                                   'outside validity range in Variable call')
                        else:
                            raise ValueError, ('Computed value '+str(depvar)+\
                                   'outside validity range in Variable call')
                else:
                    if depvar in self.depdomain:
                        return depvar
                    else:
                        if vectorizable:
                            raise ValueError, ('Computed value(s) '
                                    'outside validity range in Variable call')
                        else:
                            raise ValueError, ('Computed value '+str(depvar)+\
                                    'outside validity range in Variable call')
            except PyDSTool_BoundsError, e:
                raise ValueError, ("Cannot compute a return value for "
                                      "this independent variable value: "
                                      + str(e))
            except PyDSTool_TypeError:
                if not self.defined:
                    print "Variable '%s' not fully defined."%self.name
                    return None
                else:
                    raise


    def getDataPoints(self):
        """Reveal underlying mesh and values at mesh points, provided
        Variable is based on a mesh (otherwise None is returned)"""
        if isinstance(self.output, Varcaller):
            return {self.indepvarname: self.output.pts.indepvararray, \
                    self.coordname: self.output.pts.toarray()}
        elif hasattr(self.output, 'datapoints'):
            datapoints = self.output.datapoints
            return {self.indepvarname: datapoints[0], \
                    self.coordname: datapoints[1]}
        else:
            return None


    def __repr__(self):
        return self._infostr(verbose=0)


    __str__ = __repr__


    def _infostr(self, verbose=1):
        if verbose == 0:
            return "Variable "+self.coordname+"("+self.indepvarname+")"
        else:
            try:
                if isinputcts(self):
                    ipstr = "continuous"
                else:
                    ipstr = "discrete"
            except ValueError:
                ipstr = "not defined"
            outputStr = "Variable:\n  Independent variable '" \
                        + self.indepvarname + "' [" + ipstr + "]\n"
            try:
                if isoutputcts(self):
                    opstr = "continuous"
                else:
                    opstr = "discrete"
            except ValueError:
                opstr = "not defined"
            outputStr += "    defined in domain  " + str(self.indepdomain)
            if verbose == 2:
                if self.trajirange is None:
                    outputStr += "\n    ranges not known for this trajectory"
                else:
                    outputStr += "\n    trajectory ranges  "+str(self.trajirange)
            outputStr += "\nDependent variable '" + self.coordname + \
                        "' [" + opstr + "]\n    defined in domain  "
            if not isinstance(self.depdomain, Interval):
                outputStr += _pytypename[self.coordtype]+": "
            outputStr += str(self.depdomain)
            if verbose == 2:
                if self.trajdrange is None:
                    outputStr += "\n    ranges not known for this trajectory"
                else:
                    outputStr += "\n    trajectory ranges  "+str(self.trajdrange)
            return outputStr


    def info(self, verboselevel=1):
        print self._infostr(verboselevel)


    def clearWarnings(self):
        self.warnings = []


    def showWarnings(self):
        if self.warnings != []:
            print "Warnings:"
            for (i,d) in self.warnings:
                if d is None:
                    print "Independent variable value %s was out of "% i + \
                          "bounds"
                else:
                    print "Dependent variable value %s was out of " % s + \
                          "bounds at independent variable value %s" % i


    def __copy__(self):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


    def __deepcopy__(self, memo=None, _nil=[]):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


    def __getstate__(self):
        d = copy.copy(self.__dict__)
        # remove reference to Cfunc types by converting to strings
        d['indepvartype'] = _pytypename[self.indepvartype]
        d['coordtype'] = _pytypename[self.coordtype]
        if 'funcspec' in self._funcreg:
            # then self is Imp/ExplicitFnGen and 'output' could not
            # be put in _funcreg because it relies on wrap_output
            # function that's not in the global namespace (so pickle fails
            # to find it)
            del d['output']
        for fname, finfo in self._funcreg.iteritems():
            if finfo[0] == 'self':
                try:
                    del d[fname]
                except KeyError:
                    pass
            # else it's a Variable class method which won't get pickled
            # anyway, and will be restored to any class not in possession
            # of it if this object is unpickled
        return d


    def __setstate__(self, state):
        self.__dict__.update(state)
        #print self.name, "- setstate: self.depdomain = ", self.depdomain.get()
        # reinstate Cfunc types
        self.indepvartype = _pynametype[self.indepvartype]
        self.coordtype = _pynametype[self.coordtype]
        # reinstate dynamic methods / functions
        for fname, finfo in self._funcreg.iteritems():
            if finfo[0] == 'self' and not hasattr(eval(finfo[0]), fname):
                # avoids special entry for 'outputdata'
                setattr(eval(finfo[0]), fname, finfo[1])
        if 'funcspec' in self._funcreg:
            # Add the specific mapping functions for Ex/ImplicitFnGen objects
            funcspec = self._funcreg['funcspec'][1]
            outputdata = self._funcreg['outputdata'][1]
            if hasattr(self, '_var_namemap'):
                var_namemap = self._var_namemap
            else:
                var_namemap = None
            if hasattr(self, 'initialconditions'):
                ics = copy.copy(self.initialconditions)
            else:
                ics = None
            if hasattr(self, '_refvars'):
                if self._refvars is not None and self._refvars != []:
                    refvars = [copy.copy(v) for v in self._refvars]
                else:
                    refvars = None
            else:
                refvars = None
            # if refvars in dictionary then just leave them there!
            self.setOutput(outputdata, funcspec,
              self.globalt0, var_namemap, ics, refvars)


    def __del__(self):
        # delete object-specific class methods etc. before deleting
        # to avoid crowding namespace
##        if hasattr(self, 'output'):
##            del self.output
        for fname, finfo in self._funcreg.iteritems():
            # Treat special cases first
            if finfo[0] is None:
                # don't want to eval(None) below
                continue
            elif fname == '_impfn':
                exec_str = 'del Variable.' + finfo[0]
                try:
                    exec exec_str
                except AttributeError:
                    # Uncertain why the name appears multiple times for their
                    # to be multiple attempts to delete it (which of course
                    # fail after the first successful attempt)
                    pass
            elif fname is 'funcspec':
                # doesn't refer to any dynamically-created methods
                # so ignore
                pass
            elif fname is 'outputdata':
                # doesn't refer to any dynamically-created methods
                # so ignore
                pass
            elif hasattr(eval(finfo[0]), fname):
                exec_str = 'del '+ finfo[0] + '.' + fname
                try:
                    exec exec_str
                except RuntimeError:
                    # sometimes get these when objects improperly delted
                    # and new objects with the same name created
                    pass
        if hasattr(self, '_refvars'):
                if self._refvars is not None and self._refvars != []:
                    for v in self._refvars:
                        v.__del__()


class OutputFn(object):
    """One-dimensional function wrapper."""
    
    def __init__(self, fn, datapoints=None, numtypes=(Float,Float)):
        assert isinstance(fn, types.FunctionType) or \
               isinstance(fn, types.BuiltinFunctionType), \
               ("fn argument must be a regular Python function")
        self.fn = fn
        # datapoints can be exhaustive list of known values for fn or
        # a Interval range for continuous-valued functions
        if datapoints is None:
            datapoints = (Interval('indepvardom', _pynametype[numtypes[0]],
                                      [-Inf, Inf]),
                          Interval('depvardom', _pynametype[numtypes[1]],
                                      [-Inf, Inf]))
        try:
            self.datapoints = (datapoints[0], datapoints[1])
        except TypeError:
            raise TypeError, ("datapoints argument must be a 2-tuple or list "
                              "of 2-tuples or lists")
        try:
            self.types = (numtypes[0], numtypes[1])
        except TypeError:
            raise TypeError, ("numtypes argument must be a 2-tuple or list "
                              "of 2-tuples or lists")


    def __call__(self, arg):
        if type(arg) in [list, Array, NArray]:
            try:
                return self.fn(arg)
            except:
                return array([self.fn(v) for v in arg])
        else:
            return self.fn(arg)


    def __getstate__(self):
        d = copy.copy(self.__dict__)
        # remove reference to Cfunc types by converting to strings
        d['types'] = (_pytypename[self.types[0]],_pytypename[self.types[1]])
        return d


    def __setstate__(self, state):
        self.__dict__.update(state)
        # reinstate Cfunc types
        self.types = (_pynametype[self.types[0]],_pytypename[self.types[1]])


# ---------------------------------------------------------------------


def isinputcts(obj):
    if isinstance(obj, Variable):
        if obj.defined:
            if obj.indepvartype == Float:
                return isinstance(obj.indepdomain, Interval) and not \
                       isinstance(obj.output, Pointset)
            elif obj.indepvartype == Int:
                return False
            else:
                raise TypeError, "Unsupported independent variable type for Variable"
        else:
            raise ValueError, "Variable is not fully defined"
    else:
        # provide support for e.g. Trajectories. Cannot use Trajectory class
        # name explicitly here because will run into an infinite import loop
        # between Variable and Trajectory!
        if obj.indepvartype == Float:
            return isinstance(obj.indepdomain, Interval)


def isinputdiscrete(var):
    return not isinputcts(var)

##def isinputdiscrete(var):
##    if var.indepvartype == Float:
##        return type(var.indepdomain) in [NArray, Array] or \
##               isinstance(var.output, Pointset)
##    elif var.indepvartype == Int:
##        return True
##    else:
##        raise TypeError, "Unsupported independent variable type for Variable"


def isoutputcts(var):
    assert isinstance(var, Variable), "Argument must be a Variable"
    if var.defined:
        if var.coordtype == Float:
            return isinstance(var.depdomain, Interval) and not \
                   isinstance(var.output, Pointset)
        elif var.coordtype == Int:
            return False
        else:
            raise TypeError, "Unsupported dependent variable type for Variable"
    else:
        raise ValueError, "Variable is not fully defined"


def isoutputdiscrete(obj):
    return not isoutputcts(obj)


def iscontinuous(var):
    """Determine if variable is continuously defined on its input and
    output domains."""
    
    assert isinstance(var, Variable), "Argument must be a Variable"
    return isinputcts(var) and isoutputcts(var)


def isdiscrete(var):
    """Determine if variable is discretely defined on its input and
    output domains."""

    return not (isinputcts(var) and isoutputcts(var))



# ---------------------------------------------------------------------

if __name__ == '__main__':
    from Points import Pointset

    w_pts = Pointset({'coordarray': array([4.456, 2.34634, 7.3431, 5.443], Float),
                  'indepvararray': array([0.0, 1.0, 2.0, 3.0], Float)})
    w_var = Variable(w_pts)
    print "w_var(0.0) => ", w_var(0.0)

    print "\n"
    f2 = interp1d([0., 1., 2.], [5,6,7])
    v_var = Variable(f2, 't', 'x')
    print "Use optional 'checklevel' argument to specify degree of bounds/domain checking"
    print "This is useful mainly during initial computation of a variable"
    print "v_var(0.01) => ", v_var(0.01, 2)
    try:
        print "\nv_var(array([0.24,2.566]), 2) =>\n", v_var(array([0.24,2.566]), 2)
    except ValueError, e:
        print " ",e
    print "\nv_var =>\n", v_var
    print "v_var.getDataPoints() => ", v_var.getDataPoints()

    # ------------------------------

    print """Test of domain checking... By making depdomain to be a
    set of integers, and the coordarray to be non-integers, this
    object can only be called for indepvar values in [0, 0.5, 1, 1.5,
    ..., 4.5]"""
    
    print "\n"
    v_int = Variable(Pointset({'coordarray': array(range(10), Float)*0.1,
                            'indepvararray': array(range(10), Float)*0.5
                            }))
    print "v_int(0.5) => ", v_int(0.5, 2)
    print "v_int(0.4) => "
    try:
        v_int(0.4, 2)
    except ValueError, e:
        print "Successfully checked domain validation with v_int:"
        print " ... error was "+str(e)
    print "v_int(0) => ", v_int(0)

    # ------------------------------

    print "\nTest simple functions in Variable object"
    exp_str = """exp_var = Variable(math.exp, 'x', Interval('y', 'float', [0,Inf]))"""
    print exp_str
    exec(exp_str)
    print "exp_var(0.5) => ", exp_var(0.5)

    print "\nTest wrapped functions in OutputFn class"
    print """The optional Intervals specify the "trajectory" range, but are for
    informational purposes only! They are not checked anywhere."""
    sin_str = """sin_opfunc = OutputFn(math.sin, (Interval('t', 'float', [0,Inf]),
                                     Interval('x', 'float', [-1.,1.])))"""
    print sin_str
    exec(sin_str)
    print """\nThese Intervals specify the valid domains of the indep and dep var
Deliberately only allow +ve angles (unless specified here they won't be checked)"""
    sin_str2 = """sin_var = Variable(sin_opfunc, Interval('t', 'float', [0,Inf]),
                    Interval('x', 'float', [-1,1]))"""
    print sin_str2
    exec(sin_str2)
    print "sin_var(math.pi) => ", sin_var(math.pi, 2)
    print "sin_var(-math.pi/2) =>"
    try:
        sin_var(-math.pi/2, 2)
    except ValueError, e:
        print " ", e

    print "sin_var.getDataPoints() => ", sin_var.getDataPoints()
    print "sin_var([0., 0.5*math.pi, math.pi]) => ", sin_var([0., 0.5*math.pi,
                                                              math.pi])

    # ------------------------------

    from utils import makeImplicitFunc
    print "\nTest implicit function routine on half-circle of radius 2"
    cf = """def circ_formula(x,y):
        return x*x+y*y-4"""

    print cf
    exec(cf)

    nf_str = """newton_halfcirc_fn = makeImplicitFunc(circ_formula, x0=0.75, solmethod='newton')"""
    print nf_str
    exec(nf_str)
    
    impl_str = """implicit_halfcirc = OutputFn(newton_halfcirc_fn, (Interval('t',
    'float', (0,2)), Interval('x', 'float', (-2,2))))"""
    print impl_str
    exec(impl_str)

    print "tval = -1.3"
    tval = -1.3
    xval = implicit_halfcirc(tval)
    print "xval = implicit_halfcirc(tval) => ", xval
    print "math.sqrt(xval*xval + tval*tval) => ", math.sqrt(xval*xval + \
                                                            tval*tval), " = radius"
    print "As it stands, the OutputFn doesn't understand the bounds on the variables:"
    print "implicit_halfcirc(3.) =>"
    try:
        implicit_halfcirc(3.)
    except RuntimeError,e:
        print " ... returns error: ", e
    except AssertionError,e:
        print " ... returns error: ", e
    print """\nSo we can embed this OutputFn into a variable for x as
    a function of t, with enforcable bounds (using checklevel > 0 second call argument)"""

    implicit_hc_varstr = """implicit_halfcirc_var = Variable(implicit_halfcirc, Interval('t', 'float', (0,2)),
    Interval('x', 'float', (-2,2)))"""

    print implicit_hc_varstr
    exec(implicit_hc_varstr)

    print "implicit_halfcirc_var(tval) => ", implicit_halfcirc_var(tval), ", as before"
    print "implicit_halfcirc_var(3., 2) =>"
    try:
        implicit_halfcirc_var(3., 2)
    except ValueError,e:
        print " ... returns error: ", e
    print "Tests passed."
