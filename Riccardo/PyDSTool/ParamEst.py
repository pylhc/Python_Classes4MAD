#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      17 Sept 2005
#$Revision: 1.1 $
"""Parameter estimation classes for ODEs.

   Robert Clewley, February 2005.
"""

# PyDSTool imports
from Model import Model
from utils import *
from common import *
from errors import *
import Redirector as redirc

try:
    # use pscyo JIT byte-compiler optimization, if available
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

# Other imports
from scipy.optimize import minpack, optimize
from numarray import array, arange, zeros, sum, power
import math, types
from copy import copy


# ------------------------------------------------------------------------
# VERSION HISTORY:
#
# v1.1, 15 Sept 2005
#   Added support for Psyco JIT byte-compiler optimization, if available.
#     (x86 platforms only).
#   Changed dynamic method creation to use types.MethodType so method is
#     added to instance, not underlying class.
#
# ------------------------------------------------------------------------

__all__ = ['LMpest', 'BoundMin']

rout = redirc.Redirector(redirc.STDOUT)
rerr = redirc.Redirector(redirc.STDERR)


# This default function is not used as a regular method
# It may be assigned as the function bound to an object as a method
def residual_fn_LMpest_default_Ndim(s, p, tmesh):
    # p comes in as numarray
    for i in xrange(s.numFreePars):
        s.modelArgs[s.parTypeStr[i]][s.freeParNames[i]] = p[i]
    s.testModel.compute(**s.modelArgs)
    for i in xrange(s._depvar_len):
        vec = s.goalTraj(tmesh,s.depVarNames[i]) - s.testModel(s.trajname,
                                                          tmesh,
                                                          s.depVarNames[i])
        s._memory[0+i*s._tmesh_len:s._tmesh_len+i*s._tmesh_len] = vec.toarray()
    return s._memory


def residual_fn_LMpest_default_1dim(s, p, tmesh):
    # p comes in as numarray
    # in 1D case, s.depVarNames is of length 1
    for i in xrange(s.numFreePars):
        s.modelArgs[s.parTypeStr[i]][s.freeParNames[i]] = p[i]
    s.testModel.compute(**s.modelArgs)
    vec = s.goalTraj(tmesh,s.depVarNames) - s.testModel(s.trajname, tmesh,
                                                          s.depVarNames)
    return vec.toarray()


def residual_fn_BoundMin_default(s, p, tmesh):
    # uses L2 norm of error (1-D only)
    s.modelArgs[s.parTypeStr][s.freeParNames[0]] = p
    s.testModel.compute(**s.modelArgs)
    vec = s.goalTraj(tmesh,s.depVarNames) - s.testModel(s.trajname, tmesh,
                                                          s.depVarNames)
    return math.sqrt(reduce(float.__add__, power(vec.toarray(),2)))
    

# ----------------------------------------------------------------------------


class ParamEst(Utility):
    """General-purpose parameter estimation class."""

    def __init__(self, **kw):
        self.needKeys = ['freeParams', 'depVarNames', 'testModel',
                         'goalTraj', 'trange']
        self.optionalKeys = ['trajName', 'usePsyco']
        try:
            self.goalTraj = kw['goalTraj']
            assert isinstance(kw['depVarNames'], list), \
                           ("'depVarNames' argument must be a list")
            self.depVarNames = kw['depVarNames']
            assert isinstance(kw['freeParams'], list), ("'freeParams' argument"
                                                " must be a list")
            self.freeParNames = kw['freeParams']
            self.numFreePars = len(self.freeParNames)
            self.testModel = kw['testModel']
            assert isinstance(kw['testModel'], Model), \
                   "testModel argument must be a Model instance"
            self._algParamsSet = False
            self.trange = kw['trange']
        except KeyError:
            raise PyDSTool_KeyError, 'Incorrect argument keys passed'
        self.foundKeys = len(self.needKeys)
        if 'trajName' in kw:
            self.trajname = kw['trajName']
            self.foundKeys += 1
        else:
            print "Using default trajectory name 'par_est'"
            self.trajname = 'par_est'
        if 'usePsyco' in kw:
            if HAVE_PSYCO and kw['usePsyco']:
                self.usePsyco = True
            else:
                self.usePsyco = False
            self.foundKeys += 1
        else:
            self.usePsyco = False
        if self.trajname in self.testModel.trajectories:
            if self.trajname != 'par_est':
                print "Warning: Specified trajectory name already exists" \
                      + " and will be overwritten."
            self.testModel.del_traj(self.trajname)


    def setJacobian(self):
        # Set up Jacobian, if present
        if self.testModel.haveJacobian():
            self.Jacobian = self.testModel.Jacobian
##            if self.usePsyco:
##                psyco.bind(self.Jacobian)
            return True
        else:
            self.Jacobian = None
            return False


    def setAlgParams(self, *args):
        """Set algorithm pars."""
        
        raise NotImplementedError, "This is only an abstract function definition"


    def setResidualFn(self, fn):
        setattr(self, '_residual_fn', types.MethodType(fn, self,
                                                       self.__class__))
        if self.usePsyco:
            psyco.bind(self._residual_fn)


    def run(self):
        """Returns a dictionary:

            'success' -> boolean

            'pars_fit' -> fitted values of pars

            'pars_orig' -> original values of optimized pars
            
            'sys_fit' -> trajectory of best fit Model trajectory

            'alg_results' -> all other algorithm information (list)
        """
        
        raise NotImplementedError, "This is only an abstract function definition"


    def get_tmesh(self):
        """Default tmesh generator (20 steps over time domain of goal
        trajectory)."""
        
        tlimits = self.goalTraj.indepdomain.get()
        return arange(tlimits[0], tlimits[1], 0.05*(tlimits[1]-tlimits[0]))



class LMpest(ParamEst):
    """Unconstrained least-squares parameter and initial condition optimizer
    for DS trajectories. N-dimensional in parameter space.

    Uses MINPACK Levenberg-Marquardt algorithm wrapper from SciPy.
    """

    def __init__(self, **kw):
        ParamEst.__init__(self, **kw)
        # useful quantity for inside residual functions
        self._depvar_len = len(self.depVarNames)
        self.optionalKeys.extend(['residual_fn'])
        if 'residual_fn' in kw:
            self.setResidualFn(kw['residual_fn'])
            self.foundKeys += 1
        else:
            if self._depvar_len == 1:
                self.setResidualFn(residual_fn_LMpest_default_1dim)
            else:
                self.setResidualFn(residual_fn_LMpest_default_Ndim)
        if self.foundKeys < len(kw):
            raise PyDSTool_KeyError, 'Incorrect argument keys passed'
        # Set up model arguments (parameter value will be set before needed)
        self.modelArgs = {'trajname': self.trajname,
                          'force': True,
                          'tdata': self.trange}
        self.parTypeStr = []
        for i in xrange(self.numFreePars):
            if self.freeParNames[i] in self.testModel.obsvars:
                # for varying initial conditions
                self.parTypeStr.append('ics')
                if 'ics' not in self.modelArgs:
                    self.modelArgs['ics'] = {}
            elif self.freeParNames[i] in self.testModel.pars:
                # for varying regular pars
                self.parTypeStr.append('pars')
                if 'pars' not in self.modelArgs:
                    self.modelArgs['pars'] = {}
            else:
                raise ValueError, ("free parameter '"+self.freeParNames[i]+"'"\
                                   " not found in test model")
            # initialize model argument to None
            self.modelArgs[self.parTypeStr[i]][self.freeParNames[i]] = None
        found_jac = self.setJacobian()


    def setAlgParams(self, changed_parDict={}):
        parDict = {
                 'residuals': None,
                 'p_start': None,
                 'args'      : None,
                 'Dfun'      : None,
                 'full_output' : 1,
                 'col_deriv'   : 0,
                 'ftol'        : 5e-5,
                 'xtol'        : 5e-5,
                 'gtol'        : 0.0,
                 'maxfev'      : 50,
                 'epsfcn'      : 0.0,
                 'factor'      : 100,
                 'diag'        : None
                 }

        parDict.update(copy(changed_parDict))
        assert len(parDict) == 13, "Incorrect param dictionary keys used"

        self.__residuals = parDict['residuals']
        self.__p_start   = parDict['p_start']
        self.__args      = parDict['args']
        self.__Dfun        = parDict['Dfun']
        self.__full_output = parDict['full_output']
        self.__col_deriv   = parDict['col_deriv']
        self.__ftol        = parDict['ftol']
        self.__xtol        = parDict['xtol']
        self.__gtol        = parDict['gtol']
        self.__maxfev      = parDict['maxfev']
        self.__epsfcn      = parDict['epsfcn']
        self.__factor      = parDict['factor']
        self.__diag        = parDict['diag']
        # flag for run() to start
        self._algParamsSet = True


    def run(self, parDict={}, tmesh=None, verbose=False):
        parDict_new = copy(parDict)
        if tmesh is None:
            tmesh = self.get_tmesh()
        else:
            if len(tmesh) < 4 + self.numFreePars:
                print "Warning: supplied tmesh contains too few points for "\
                      + "algorithm to run accurately -- it may fail."
        if 'args' in parDict:
            parDict_new['args'] = (tmesh,)+parDict['args']
        else:
            parDict_new['args'] = (tmesh,)
        parsOrig = []
        for i in xrange(self.numFreePars):
            parsOrig.append(self.testModel.query(self.parTypeStr[i])\
                            [self.freeParNames[i]])
        parsOrig = array(parsOrig)
        parDict_new['p_start'] = parsOrig
        parDict_new['Dfun'] = self.Jacobian
        if 'residuals' not in parDict:
            parDict_new['residuals'] = self._residual_fn

        # No need to compute the trajectory here. It gets done in
        # the residual function computation.

        # Setting default minimizer pars
        if not self._algParamsSet:
            self.setAlgParams(parDict_new)

        # store these for use within residual functions
        self._tmesh_len = len(tmesh)
        self._memory = zeros(self._tmesh_len*self._depvar_len, 'Float64')
        
        # perform least-squares fitting
        rout.start()
        rerr.start()
        try:
            results = minpack.leastsq(self.__residuals,
                               self.__p_start,
                               args   = self.__args,
                               Dfun   = self.__Dfun,
                               full_output = self.__full_output,
                               col_deriv   = self.__col_deriv,
                               ftol   = self.__ftol,
                               xtol   = self.__xtol,
                               gtol   = self.__gtol,
                               maxfev = self.__maxfev,
                               epsfcn = self.__epsfcn,
                               factor = self.__factor,
                               diag   = self.__diag)
        except:
            print "Calculating residual failed for pars:", \
                  parsOrig
            raise
        out = rout.stop()
        err = rerr.stop()

        # build return information
        success = results[2] == 1
        if isinstance(results[0], float):
            res_par_list = [results[0]]
            orig_par_list = [parsOrig[0]]
        else:
            res_par_list = results[0].tolist()
            orig_par_list = parsOrig.tolist()
        pestReturn = {'success': success,
                      'pars_fit': dict(zip(self.freeParNames,
                                           res_par_list)),
                      'pars_orig': dict(zip(self.freeParNames,
                                            orig_par_list)),
                      'alg_results': results[1],
                      'sys_fit': self.testModel
                      }

        if verbose:
            if success or results[3].find('at most') != -1:
                if success:
                    print 'Solution of ', self.freeParNames, ' = ', results[0]
                else:
##                    parvals = [self.testModel.pars[p] for p in \
##                               self.freeParNames]
                    print 'Closest values of ', self.freeParNames, ' = ', \
                          results[0]
##                          parvals
                print 'Original values = ', parsOrig
                print 'Number of mesh points = ', len(tmesh)
                print 'Number of fn evals = ', results[1]["nfev"], "(# iterations)"
                if not success:
                    print 'Solution not found: '+results[3]
            else:
                print 'Solution not found: '+results[3]
        return pestReturn



class BoundMin(ParamEst):
    """Bounded minimization parameter and initial condition optimizer
    for DS trajectories. 1 parameter only.

    Uses SciPy.optimize fminbound algorithm.
    """

    def __init__(self, **kw):
        assert len(kw['freeParams']) == 1, ("Only one free parameter can "
                                            "be specified for this class")
        assert len(kw['depVarNames']) == 1, \
            ("Only one dependent variable can be specified for this class")
        ParamEst.__init__(self, **kw)
        self.optionalKeys.extend(['residual_fn'])
        if 'residual_fn' in kw:
            self.setResidualFn(kw['residual_fn'])
            self.foundKeys += 1
        else:
            self.setResidualFn(residual_fn_BoundMin_default)
        if self.foundKeys < len(kw):
            raise PyDSTool_KeyError, 'Incorrect argument keys passed'
        if self.freeParNames[0] in self.testModel.obsvars:
            # for varying initial conditions
            self.parTypeStr = 'ics'
        elif self.freeParNames[0] in self.testModel.pars:
            # for varying regular pars
            self.parTypeStr = 'pars'
        else:
            raise ValueError, 'free parameter name not found in test model'
        # Set up model arguments (parameter value will be set before needed)
        self.modelArgs = {'trajname': self.trajname,
                          'force': True,
                          'tdata': self.trange,
                          self.parTypeStr: {self.freeParNames[0]: None}}


    def run(self, parConstraints, xtol=5e-5, maxiter=500,
            tmesh=None, extra_args=(), verbose=False):
        parsOrig = self.testModel.query(self.parTypeStr)[self.freeParNames[0]]

        # default time mesh if none supplied (t range split into 20)
        if tmesh is None:
            tmesh = self.get_tmesh()

        full_output = 1
        rout.start()
        rerr.start()
        results = optimize.fminbound(self._residual_fn, parConstraints[0],
                        parConstraints[1], (tmesh,) + extra_args, xtol, maxiter,
                                     full_output,
                                     int(verbose))
        out = rout.stop()
        err = rerr.stop()

        # build return information
        success = results[2] == 0
        pestReturn = {'success': success,
                      'pars_fit': {self.freeParNames[0]: results[0]},
                      'pars_orig': {self.freeParNames[0]: parsOrig},
                      'alg_results': results[3],
                      'sys_fit': self.testModel
                      }

        if verbose:
            if success:
                print 'Solution of ', self.freeParNames[0], ' = ', results[0]
                print 'Original value = ', parsOrig
                print 'Number of mesh points = ', len(tmesh)
                print 'Number of fn evals = ', results[3], "(# iterations)"
                print 'Error tolerance = ', xtol
            else:
                print 'No convergence of BoundMin'
                print results
        return pestReturn

# ----------------------------------------------------------------------------

