# Radau ODE system (no DAE support)
from __future__ import division

from allimports import *
from PyDSTool.Generator import ODEsystem as ODEsystem
from baseclasses import Generator, theGenSpecHelper, genDB
from PyDSTool.utils import *
from PyDSTool.common import *
from PyDSTool.ModelSpec import QuantSpec
from PyDSTool.integrator import integrator
from PyDSTool.parseUtils import addArgToCalls, wrapArgInCall
import PyDSTool.Redirector as redirc

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue, any, all
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, transpose, shape
import math, random
from copy import copy, deepcopy
import os, platform, shutil, sys, gc
# Use scipy_distutils that supports fortran compiler
import numpy.distutils
from numpy.distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc

# path to the installation
import PyDSTool
_pydstool_path = PyDSTool.__path__[0]

# -----------------------------------------------------------------------------
# VERSION HISTORY
#
# 19 September 2006
#   Added stdout redirection for calls to minpack functions.
#   Added support for user-modified bounds on state variables.
#   Made _inputVarList and _inputTimesList deep copied in passing to integrator
#     in order to avoid weird side effect where they are lengthened by a factor
#     of # inputs, and the remaining entries filled with zero. This means that
#     if the generator gets called again then the integrator will complain that
#     the input times are not properly ordered.
#   Also made _inputTimesList relative to self.globalt0 so that external
#     inputs work correctly in hybrid model situations.
#
# 25 February 2006
#   Added generator registry to avoid DLL re-usage conflicts.
#
# 8 February 2006
#   Fixed backwards integration.
#
# 22 January 2006
#   Added trajevents structure.
#
# 15 November 2005, 12 December 2005
#   New version of integrator.
#
# 14 October 2005
#   New version of Dopri wrapper, makes CleanUp() work for better memory
#     handling.
# -----------------------------------------------------------------------------

rout = redirc.Redirector(redirc.STDOUT)
rerr = redirc.Redirector(redirc.STDERR)


class radau(integrator):
    """Radau 5 specialization of the basic integrator class."""

    def __init__(self, modname, rhs='default_name', phaseDim=0, paramDim=0,
                 nAux=0, nEvents=0, nExtInputs=0, hasJac=0, hasJacP=0,
                 hasMass=0, extraSpace=0, defaultBound=1e8):
        integrator.__init__(self, rhs=rhs, phaseDim=phaseDim, paramDim=paramDim,
                            nAux=nAux, nEvents=nEvents, nExtInputs=nExtInputs, hasJac=hasJac,
                            hasJacP=hasJacP, hasMass=hasMass, extraSpace=extraSpace,
                            defaultBound=defaultBound)
        self.modname = modname
        try:
            self._integMod = __import__(modname, globals())
        except:
            print "Error in importing compiled vector field and integrator."
            print "Did you compile the RHS C code?"
            raise
        # check module's directory
        assert 'Integrate' in dir(self._integMod), \
               "radau library does not contain Integrate()"

        self.safety = []
        self.jacRecompute = []
        self.newtonStop = []
        self.stepChangeLB = []
        self.stepChangeUB = []
        self.stepSizeLB = []
        self.stepSizeUB = []
        self.hessenberg = []
        self.maxNewton = []
        self.newtonStart = []
        self.index1dim = []
        self.index2dim = []
        self.index3dim = []
        self.stepSizeStrategy = []
        self.DAEstructureM1 = []
        self.DAEstructureM2 = []
        
        retval = self._integMod.InitBasic(self.phaseDim, self.paramDim, self.nAux,
                                          self.nEvents, self.nExtInputs, self.hasJac,
                                          self.hasJacP, self.hasMass, self.extraSpace)
 
        if retval[0] != 1:
            raise InitError, 'Call to InitBasic failed! (radau)'

        self.initBasic = True


    def Run(self, hinit=0, hmax=1.0, checkAux=0, calcSpecTimes=0, verbose=0,
            safety=0.9, jacRecompute=0.001, newtonStop=-1, stepChangeLB=1,
            stepChangeUB=1.2, stepSizeLB=0.2, stepSizeUB=8.0, hessenberg=0,
            maxNewton=7, newtonStart=0, index1dim=-1, index2dim=0, index3dim=0,
            stepSizeStrategy=1, DAEstructureM1=0, DAEstructureM2=0, useJac=0, useMass=0):
        if not self.initBasic:
            raise InitError, 'initBasic is False (radau)'
        if not self.initEvents:
            raise InitError, 'initEvents is False (radau)'
        if not self.initIntegrate:
            raise InitError, 'initInteg is False (radau)'
        if not self.setParams:
            raise InitError, 'setParams is False (radau)'
        if self.nExtInputs > 0 and not self.initExtInputs:
            raise InitError, 'initExtInputs is False (radau)'

        self.setRadauParams(hinit=hinit, hmax=hmax, checkAux=checkAux,
                            calcSpecTimes=calcSpecTimes,
                            verbose=verbose, safety=safety,
                            jacRecompute=jacRecompute, newtonStop=newtonStop,
                            stepChangeLB=stepChangeLB, stepChangeUB=stepChangeUB,
                            stepSizeLB=stepSizeLB, stepSizeUB=stepSizeUB,
                            hessenberg=hessenberg,maxNewton=maxNewton,
                            newtonStart=newtonStart, index1dim=index1dim,
                            index2dim=index2dim, index3dim=index3dim,
                            stepSizeStrategy=stepSizeStrategy,
                            DAEstructureM1=DAEstructureM1,
                            DAEstructureM2=DAEstructureM2,
                            useJac=useJac,useMass=useMass)

        # For a run, we want to ensure indices are set to 0
        self.Reset()
        T, P, A, Stats, H, Err, EvtT, EvtP = self._integMod.Integrate(self.ic,
                                                      self.t0,
                                                      self.hinit,
                                                      self.hmax,
                                                      self.safety,
                                                      self.jacRecompute,
                                                      self.newtonStop,
                                                      self.stepChangeLB,
                                                      self.stepChangeUB,
                                                      self.stepSizeLB,
                                                      self.stepSizeUB,
                                                      self.hessenberg,
                                                      self.maxNewton,
                                                      self.newtonStart,
                                                      self.index1dim,
                                                      self.index2dim,
                                                      self.index3dim,
                                                      self.stepSizeStrategy,
                                                      self.DAEstructureM1,
                                                      self.DAEstructureM2,
                                                      self.useJac,
                                                      self.useMass,
                                                      self.verbose,
                                                      self.checkAux,
                                                      self.calcSpecTimes)
        self.points = P
        self.times = T
        self.auxPoints = A
        self.eventTimes = EvtT
        self.eventPoints = EvtP
        self.errors = Err
        self.stats = Stats
        self.step = H
        
        try:
            self.lastTime = self.times[-1]
            self.lastPoint = [self.points[i][-1] for i in range(self.phaseDim)]
            self.lastStep = self.step
        except IndexError:
            self.lastTime = self.t0
            self.lastPoint = self.ic
            self.lastStep = self.hinit
        self.numRuns += 1
        self.canContinue = True

        return T, P, A, Stats, H, Err, EvtT, EvtP

 
    def Continue(self, tend, params=[], calcSpecTimes=0, verbose=0, extInputChanged=False,
                 extInputVals=[], extInputTimes=[], bounds=[]):
        if not self.initBasic:
            raise InitError, 'initBasic is False (radau)'
        if not self.initEvents:
            raise InitError, 'initEvents is False (radau)'
        if not self.initIntegrate:
            raise InitError, 'initInteg is False (radau)'
        if not self.setParams:
            raise InitError, 'setParams is False (radau)'
        if self.nExtInputs > 0 and not self.initExtInputs:
            raise InitError, 'initExtInputs is False (radau)'

        if not self.canContinue:
            raise ContError, 'Unable to continue trajectory -- have you run the integrator and reset events, etc?'

        self.setContParams(tend=tend, params=copy(params),
                           calcSpecTimes=calcSpecTimes, verbose=verbose, extInputChanged=extInputChanged,
                           extInputVals=copy(extInputVals), extInputTimes=copy(extInputTimes),
                           bounds=copy(bounds))

        # For a continue, we do not set indices to 0
        T, P, A, Stats, H, Err, EvtT, EvtP = self._integMod.Integrate(self.lastPoint,
                                                                      self.lastTime,
                                                                      self.lastStep, self.hmax,
                                                                      self.safety,
                                                                      self.jacRecompute,
                                                                      self.newtonStop,
                                                                      self.stepChangeLB,
                                                                      self.stepChangeUB,
                                                                      self.stepSizeLB,
                                                                      self.stepSizeUB,
                                                                      self.hessenberg,
                                                                      self.maxNewton,
                                                                      self.newtonStart,
                                                                      self.index1dim,
                                                                      self.index2dim,
                                                                      self.index3dim,
                                                                      self.stepSizeStrategy,
                                                                      self.DAEstructureM1,
                                                                      self.DAEstructureM2,
                                                                      self.useJac,
                                                                      self.useMass,
                                                                      self.verbose,
                                                                      self.checkAux,
                                                                      self.calcSpecTimes)
        
        self.points = P
        self.times = T
        self.auxPoints = A
        self.eventTimes = EvtT
        self.eventPoints = EvtP
        self.errors = Err
        self.stats = Stats
        self.step = H
        
        try:
            self.lastTime = self.times[-1]
            self.lastPoint = [self.points[i][-1] for i in range(self.phaseDim)]
            self.lastStep = self.step
        except IndexError:
            self.lastTime = self.t0
            self.lastPoint = self.ic
            self.lastStep = self.hinit
        self.numRuns += 1
        self.numContinues += 1
        self.canContinue = True

        return T, P, A, Stats, H, Err, EvtT, EvtP


    def setRadauParams(self, hinit, hmax, checkAux, calcSpecTimes,
                       verbose, safety, jacRecompute, newtonStop,
                       stepChangeLB, stepChangeUB, stepSizeLB, stepSizeUB,
                       hessenberg, maxNewton, newtonStart, index1dim,
                       index2dim, index3dim, stepSizeStrategy,
                       DAEstructureM1, DAEstructureM2, useJac, useMass):
        useJac = int(useJac)
        useMass = int(useMass)
        checkAux = int(checkAux)
        calcSpecTimes = int(calcSpecTimes)
        hessenberg = int(hessenberg)

        if type(hinit) not in [int, float]:
            raise TypeError("hinit must be int, float")
        if hinit < 0:
            raise ValueError("hinit must be non-negative")
            
        if type(hmax) not in [int, float]:
            raise TypeError("hmax must be int, float")
        if hmax <=0:
            raise ValueError("hmax must be positive")

        if hinit > hmax:
            raise ValueError("hinit must be less than hmax")

        if type(checkAux) is not int:
            raise TypeError("checkAux must be 0 or 1")
        if checkAux not in [0,1]:
            raise TypeError("checkAux must be 0 or 1")
        if checkAux == 1 and self.nAux <= 0:
            raise ValueError("checkAux cannot be 1 if nAux is 0")
            
        if type(verbose) is not int:
            raise TypeError("verbose must be 0 or 1")
        if verbose not in [0,1]:
            raise TypeError("verbose must be 0 or 1")

        if type(calcSpecTimes) is not int:
            raise TypeError("calcSpecTimes must be 0 or 1")
        if calcSpecTimes not in [0,1]:
            raise TypeError("calcSpecTimes must be 0 or 1")
        if calcSpecTimes == 1 and len(self.specTimes) <= 0:
            raise ValueError("calcSpecTimes cannot be 1 if specTimes is empty")

        if safety < 0:
            raise ValueError("safety must be non-negative")
        if jacRecompute <= 0.0:
            raise ValueError("jacRecompute must be positive")
        if newtonStop <= 0:
            newtonStop = min(0.03,sqrt(self.rtol[0]))
        if stepChangeLB <= 0:
            raise ValueError("stepChangeLB must be positive")
        if stepChangeUB <= 0:
            raise ValueError("stepChangeUB must be positive")
        if stepSizeLB <= 0:
            raise ValueError("stepSizeLB must be positive")
        if stepSizeUB <= 0:
            raise ValueError("stepSizeUB must be positive")

        if stepChangeLB >= stepChangeUB:
            raise ValueError("stepChangeLB must be less than stepChangeUB")
        if stepSizeLB >= stepSizeUB:
            raise ValueError("stepSizeLB must be less than stepSizeUB")

        if hessenberg not in [0,1]:
            raise ValueError("hessenberg must be 0 or 1")
        if hessenberg == 1 and useMass != 0:
            raise ValueError("hessenberg form cannot be used for implicit systems (mass matrix)")
        if type(maxNewton) is not int:
            raise TypeError("maxNewton must be int")
        if maxNewton <= 0:
            raise ValueError("maxNewton must be positive")

        if newtonStart not in [0,1]:
            raise ValueError("newtonStart must be 0 or 1")

        if index1dim <= 0:
            index1dim = self.phaseDim
        if index2dim != 0:
            raise ValueError("Currently index2dim must be 0")
        if index3dim != 0:
            raise ValueError("Currently index3dim must be 0")

        if stepSizeStrategy not in [1,2]:
            raise ValueError("stepSizeStrategy must be 1 or 2")

        if DAEstructureM1 != 0:
            raise ValueError("Currently DAEstructureM1 must be 0")
        if DAEstructureM2 != 0:
            raise ValueError("Currently DAEstructureM2 must be 0")

        if useJac not in [0,1]:
            raise ValueError("useJac must be 0 or 1")
        if useMass not in [0,1]:
            raise ValueError("useMass must be 0 or 1")

        if useJac == 1 and self.hasJac != 1:
            raise ValueError("useJac must be 0 if hasJac is not 1")
        if useMass == 1 and self.hasMass != 1:
            raise ValueError("useMass must be 0 if hasMass is not 1")

        self.hinit = hinit
        self.hmax = hmax
        self.safety = safety
        self.jacRecompute = jacRecompute
        self.newtonStop = newtonStop
        self.stepChangeLB = stepChangeLB
        self.stepChangeUB = stepChangeUB
        self.stepSizeLB = stepSizeLB
        self.stepSizeUB = stepSizeUB
        self.hessenberg = hessenberg
        self.maxNewton = maxNewton
        self.newtonStart = newtonStart
        self.index1dim = index1dim
        self.index2dim = index2dim
        self.index3dim = index3dim
        self.stepSizeStrategy = stepSizeStrategy
        self.DAEstructureM1 = DAEstructureM1
        self.DAEstructureM2 = DAEstructureM2
        self.useJac = useJac
        self.useMass = useMass
        self.verbose = verbose
        self.checkAux = checkAux
        self.calcSpecTimes = calcSpecTimes



class Radau_ODEsystem(ODEsystem):
    """Wrapper for Radau integrator (with no DAE support).

    Uses C functional specifications only."""

    def __init__(self, kw):
        """Use the nobuild key to postpone building of the library, e.g. in
        order to provide additional build options to makeLibSource and
        compileLib methods or to make changes to the C code by hand.
        No build options can be specified otherwise."""

        if 'nobuild' in kw:
            nobuild = kw['nobuild']
            del kw['nobuild']
        else:
            nobuild = False
        ODEsystem.__init__(self, kw)
        self._solver = None
        assert self.funcspec.targetlang == 'c', \
               ('Wrong target language for functional specification. '
                'C needed for this class')
        assert isinstance(self.funcspec, RHSfuncSpec), ('Radau '
                                    'requires RHSfuncSpec type to proceed')
        self._errorcodes = {0: 'Unrecognized error code returned (see stderr output)',
                            -1 : 'input is not consistent',
                            -2 : 'larger nmax is needed',
                            -3 : 'step size becomes too small',
                            -4 : 'the matrix is repeatedly singular (interrupted)'}
        self._paraminfo = {'rtol': 'Relative error tolerance.',
                           'atol': 'Absolute error tolerance.',
                           'safety': 'Safety factor in the step size prediction, default 0.9.',
                           'max_step': 'Maximal step size, default tend-tstart.',
                           'init_step': 'Initial step size, default is a guess computed by the function init_step.',
                           'fac1': 'Parameter for step size selection; the new step size is chosen subject to the restriction  fac1 <= new_step/old_step <= fac2. Default value is 0.2.',
                           'fac2': 'Parameter for step size selection; the new step size is chosen subject to the restriction  fac1 <= new_step/old_step <= fac2. Default value is 0.8.',
                           'stepLB': '',
                           'stepUB': '',
                           'refine': 'Refine output by adding points interpolated using the RK4 polynomial (0, 1 or 2).',
                           'step_strategy': """Switch for step size strategy;
If step_strategy=1  mod. predictive controller (Gustafsson).
If step_strategy=2  classical step size control.
The default value (for step_strategy=0) is step_strategy=1.
the choice step_strategy=1 seems to produce safer results;
for simple problems, the choice step_strategy=2 produces
often slightly faster runs.""",
                           'jac_recompute': """Decides whether the Jacobian should be recomputed;
increase jac_recompute to 0.1 say, when Jacobian evaluations
are costly. for small systems jac_recompute should be smaller
(0.001, say). negative jac_recompute forces the code to
compute the Jacobian after every accepted step.
Default 0.001.""",
                            'newton_start': "",
                            'newton_stop': "",
                            'max_newton': "Maximum number of Newton iterations to take in solving the implicit system at each step (default 7)",
                            'DAEstructureM1': "",
                            'DAEstructureM2': "",
                            'hessenberg': "",
                            'index1dim': "",
                            'index2dim': "",
                            'index3dim': "",
                            'use_special': "Switch for using special times",
                            'specialtimes': "List of special times to use during integration",
                            'check_aux': "Switch",
                            'extraspace': ""
                           }
        # currently the final four of these params are for event handling
        algparams_def = {'init_step': 0,
                        'max_step': 0,
#                        'min_step': 0,
                        'rtol': [1e-9 for i in range(self.dimension)],
                        'atol': [1e-12 for i in range(self.dimension)],
                        'fac1': 1.0,
                        'fac2': 1.2,
                        'stepLB': 0.2,
                        'stepUB': 8.0,
                        'safety': 0.9,
                        'max_pts': 10000,
                        'refine': 0,
                        'maxbisect': [],
                        'maxevtpts': 1000,
                        'eventInt': [],
                        'eventDelay': [],
                        'eventTol': [],
                        'use_special': 0,
                        'specialtimes': [],
                        'check_aux': 1,
                        'extraspace': 100,
                        'verbose': 0,
                        'jac_recompute': 0.001,
                        'step_strategy': 1,
                        'index1dim': -1,
                        'index2dim': 0,
                        'index3dim': 0,
                        'DAEstructureM1': 0,
                        'DAEstructureM2': 0,
                        'hessenberg': 0,
                        'newton_start': 0,
                        'newton_stop': -1,
                        'max_newton': 7,
                        'hasJac': 0,
                        'hasJacP': 0
                        }
        for k, v in algparams_def.iteritems():
            if k not in self._algparams:
                self._algparams[k] = v
        # verify that no additional keys are present in algparams, after
        # defaults are added above
        if len(self._algparams) != len(algparams_def):
            raise ValueError("Invalid keys present in algparams argument: " \
                     + str(remain(self._algparams.keys(),algparams_def.keys())))
        self._outputinfo = {'last_step': 'Predicted step size of the last accepted step (useful for a subsequent call to radau).',
                            'num_steps': 'Number of used steps.',
                            'num_accept': 'Number of accepted steps.',
                            'num_reject': 'Number of rejected steps.',
                            'num_fcns': 'Number of function evaluations.',
                            'num_jacs': 'Number of Jacobian evaluations.',
                            'num_dec': 'Number of LU-decompositions',
                            'num_subs': 'Number of forward-backward substitutions',
                            'errorStatus': 'Error status on completion.'
                            }
        self.outputStats = {}
        # Check for non-constant mass matrix
        if self.haveMass():
            mspec = self.funcspec.auxfns['massMatrix']
            lensig = len(mspec[1])
            body_str = mspec[0][lensig:].replace('\n','')
            qbody = QuantSpec('__body__', body_str, treatMultiRefs=False,
                              ignoreSpecial=['[',']','{','}'])
            self._const_massmat = intersect(['Y_','t'], qbody.usedSymbols) == []
        else:
            self._const_massmat = True
        # File system preparations
        thisplatform = platform.system()
        if thisplatform == 'Windows':
            self._dllext = ".pyd"
        elif thisplatform in ['Linux', 'IRIX', 'Solaris', 'SunOS', 'MacOS', 'Darwin']:
            self._dllext = '.so'
        else:
            print "Shared library extension not tested on this platform."
            print "If this process fails please report the errors to the"
            print "developers."
            self._dllext = '.so'
        self._compilation_tempdir = os.path.join(os.getcwd(),
                                                      "radau5_temp")
        if not os.path.isdir(self._compilation_tempdir):
            try:
                assert not os.path.isfile(self._compilation_tempdir), \
                     "A file already exists with the same name"
                os.mkdir(self._compilation_tempdir)
            except:
                print "Could not create compilation temp directory " + \
                      self._compilation_tempdir
                raise
        self._compilation_sourcedir = os.path.join(_pydstool_path,"integrator")
        self._vf_file = self.name+"_vf.c"
        self._vf_filename_ext = "_"+self._vf_file[:-2]
        self._prepareEventSpecs()
        if not (os.path.isfile(os.path.join(os.getcwd(),
                                "radau5"+self._vf_filename_ext+".py")) and \
                os.path.isfile(os.path.join(os.getcwd(),
                                "_radau5"+self._vf_filename_ext+self._dllext))):
            if not nobuild:
                self.makeLibSource()
                self.compileLib()
            else:
                print "Build the library using the makeLib method, or in "
                print "stages using the makeLibSource and compileLib methods."
        self._inputVarList = []
        self._inputTimeList = []


    def forceLibRefresh(self):
        """forceLibRefresh should be called after event contents are changed,
        or alterations are made to the right-hand side of the ODEs.

        Currently this function does NOT work!"""

        # (try to) free radau module from namespace
        delfiles = True
        self._solver = None
        try:
            del(sys.modules["_radau5"+self._vf_filename_ext])
            del(sys.modules["radau5"+self._vf_filename_ext])
##            del(self._integMod)
        except NameError:
            # modules weren't loaded, so nothing to do
            delfiles = False
        if delfiles:
            gc.collect()
            # still not able to delete these files!!!!! Argh!
##            if os.path.isfile(os.path.join(os.getcwd(),
##                                    "radau5"+self._vf_filename_ext+".py")):
##                os.remove(os.path.join(os.getcwd(),
##                                    "radau5"+self._vf_filename_ext+".py"))
##            if os.path.isfile(os.path.join(os.getcwd(),
##                                 "_radau5"+self._vf_filename_ext+self._dllext)):
##                os.remove(os.path.join(os.getcwd(),
##                                 "_radau5"+self._vf_filename_ext+self._dllext))
        print "Cannot rebuild library without restarting session. Sorry."
        print "Try asking the Python developers to make a working module"
        print "unimport function!"
##        self.makeLibSource()


    def _prepareEventSpecs(self):
        eventActive = []
        eventTerm = []
        eventDir = []
        eventDelay = []
        eventTol = []
        maxbisect = []
        eventInt = []
        # convert event specs (term, active, etc.) into integparam specs
        self._eventNames = self.eventstruct.sortedEventNames()
        for evname in self._eventNames:
            ev = self.eventstruct.events[evname]
            assert isinstance(ev, LowLevelEvent), ("Radau can only "
                                                "accept low level events")
        # if event 'precise' flags set to False then set their tolerances
        # to be > max_step
        maxstep = self._algparams['max_step']
        for evname in self._eventNames:
            ev = self.eventstruct.events[evname]
            eventActive.append(int(ev.activeFlag))
            eventTerm.append(int(ev.termFlag))
            eventDir.append(ev.dircode)
            eventInt.append(ev.eventinterval)
            eventDelay.append(ev.eventdelay)
            if ev.preciseFlag:
                eventTol.append(ev.eventtol)
                maxbisect.append(ev.bisectlimit)
            else:
                eventTol.append(maxstep*1.5)
                maxbisect.append(1)
        self._algparams['eventTol'] = eventTol
        self._algparams['eventDelay'] = eventDelay
        self._algparams['eventInt'] = eventInt
        self._algparams['maxbisect'] = maxbisect
        self._algparams['eventActive'] = eventActive
        self._algparams['eventTerm'] = eventTerm
        self._algparams['eventDir'] = eventDir


    def makeLib(self, libsources=[], libdirs=[], include=[]):
        """makeLib calls makeLibSource and then the compileLib method.
        To postpone compilation of the source to a DLL, call makelibsource()
        separately."""
        self.makeLibSource(include)
        self.compileLib(libsources, libdirs)


    def makeLibSource(self, include=[]):
        """makeLibSource generates the C source for the vector field specification.
        It should be called only once per vector field."""

        # Make vector field (and event) file for compilation
        assert isinstance(include, list), "includes must be in form of a list"
        # codes for library types (default is USERLIB, since compiler will look in standard library d
        STDLIB = 0
        USERLIB = 1
        libinclude = dict([('math.h', STDLIB), ('stdio.h', STDLIB), ('stdlib.h', STDLIB),
                      ('string.h', STDLIB), ('vfield.h', USERLIB), ('events.h', USERLIB)])
        include_str = ''
        for libstr, libtype in libinclude.iteritems():
            if libtype == STDLIB:
                quoteleft = '<'
                quoteright = '>'
            else:
                quoteleft = '"'
                quoteright = '"'
            include_str += "#include " + quoteleft + libstr + quoteright + "\n"
        if include != []:
            assert uniqueList(include), "list of library includes must not contain repeats"
            for libstr in include:
                if libstr in libinclude:
                    # don't repeat libraries
                    print "Warning: library '" + libstr + "' already appears in list"\
                          + " of imported libraries"
                else:
                    include_str += "#include " + '"' + libstr + '"\n'
        allfilestr = "/*  Vector field function and events for Radau integrator.\n" \
            + "  This code was automatically generated by PyDSTool, but may be modified " \
            + "by hand. */\n\n" + include_str + """
extern double *gICs;
extern double **gBds;
extern double globalt0;

static double pi = 3.1415926535897931;

"""
        pardefines = ""
##        parundefines = ""
        vardefines = ""
##        varundefines = ""
        inpdefines = ""
##        inpundefines = ""
        # sorted version of var, par, and input names
        vnames = self._var_ixmap
        pnames = self.funcspec.pars
        inames = self.funcspec.inputs
        pnames.sort()
        inames.sort()
        for i in xrange(self.numpars):
            p = pnames[i]
            # add to defines
            pardefines += self.funcspec._defstr+" "+p+"\tp_["+str(i)+"]\n"
##            # add to undefines
##            parundefines += self.funcspec._undefstr+" "+p+"\n"
        for i in xrange(self.dimension):
            v = vnames[i]
            # add to defines
            vardefines += self.funcspec._defstr+" "+v+"\tY_["+str(i)+"]\n"
##            # add to undefines
##            varundefines += self.funcspec._undefstr+" "+v+"\n"
        for i in xrange(len(self.funcspec.inputs)):
            inp = inames[i]
            # add to defines
            inpdefines += self.funcspec._defstr+" "+inp+"\txv_["+str(i)+"]\n"
##            # add to undefines
##            inpundefines += self.funcspec._undefstr+" "+inp+"\n"            
        allfilestr += "\n/* Variable, parameter, and input definitions: */ \n" \
                      + pardefines + vardefines + inpdefines + "\n"
        # preprocess event code
        allevs = ""
        if self._eventNames == []:
            numevs = 0
        else:
            numevs = len(self._eventNames)
        for evname in self._eventNames:
            ev = self.eventstruct.events[evname]
            evfullfn = ""
            assert isinstance(ev, LowLevelEvent), ("Radau can only "
                                                "accept low level events")
            evsig = ev._LLreturnstr + " " + ev.name + ev._LLargstr
            assert ev._LLfuncstr.index(';') > 1, ("Event function code "
                    "error: Have you included a ';' character at the end of"
                                            "your 'return' statement?")
            fbody = ev._LLfuncstr
            # check fbody for calls to user-defined aux fns
            # and add hidden p argument
            if self.funcspec.auxfns:
                fbody_parsed = addArgToCalls(fbody,
                                        self.funcspec.auxfns.keys(),
                                        "p_, wk_, xv_")
                if 'initcond' in self.funcspec.auxfns:
                    # convert 'initcond(x)' to 'initcond("x")' for
                    # compatibility with C syntax
                    fbody_parsed = wrapArgInCall(fbody_parsed,
                                        'initcond', '"')
            else:
                fbody_parsed = fbody
            evbody = " {\n" + fbody_parsed + "\n}\n\n"
            allevs += evsig + evbody
            allfilestr += evsig + ";\n"
        # add signature for auxiliary functions
        if self.funcspec.auxfns:
            allfilestr += "\n"
            for finfo in self.funcspec.auxfns.values():
                allfilestr += finfo[1] + ";\n"
        assignEvBody = ""
        for evix in range(numevs):
            assignEvBody += "events[%d] = &%s;\n"%(evix,self._eventNames[evix])
        allfilestr += "\nint N_EVENTS = " + str(numevs) + ";\nvoid assignEvents(" \
              + "EvFunType *events){\n " + assignEvBody  \
              + "\n}\n\nvoid auxvars(unsigned, unsigned, double, double*, double*, " \
              + "double*, unsigned, double*, unsigned, double*);\n" \
              + """void jacobian(unsigned, unsigned, double, double*, double*, double**, unsigned, double*, unsigned, double*);
void jacobianParam(unsigned, unsigned, double, double*, double*, double**, unsigned, double*, unsigned, double*);
"""
        if self.funcspec.auxvars == []:
            allfilestr += "int N_AUXVARS = 0;\n\n\n"
        else:
            allfilestr += "int N_AUXVARS = " + str(len(self.funcspec.auxvars)) \
                       + ";\n\n\n"
        if self.funcspec.inputs == []:
            allfilestr += "int N_EXTINPUTS = 0;\n\n\n"
        else:
            allfilestr += "int N_EXTINPUTS = " + str(len(self.funcspec.inputs)) \
                       + ";\n\n\n"
        allfilestr += self.funcspec.spec[0] + "\n\n"
        if self.funcspec.auxfns:
            for fname, finfo in self.funcspec.auxfns.iteritems():
                fbody = finfo[0]
                # subs _p into auxfn-to-auxfn calls (but not to the signature)
                fbody_parsed = addArgToCalls(fbody,
                                        self.funcspec.auxfns.keys(),
                                        "p_, wk_, xv_", notFirst=fname)
                if 'initcond' in self.funcspec.auxfns:
                    # convert 'initcond(x)' to 'initcond("x")' for
                    # compatibility with C syntax, but don't affect the
                    # function signature!
                    fbody_parsed = wrapArgInCall(fbody_parsed,
                                        'initcond', '"', notFirst=True)
                allfilestr += "\n" + fbody_parsed + "\n\n"
        # add auxiliary variables (shell of the function always present)
        # add event functions
        allfilestr += self.funcspec.auxspec[0] + allevs
        # if jacobians or mass matrix not present, fill in dummy
        if not self.haveMass():
            allfilestr += """
void massMatrix(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}
"""
        if not self.haveJacobian():
            allfilestr += """
void jacobian(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}
"""
        if not self.haveJacobian_pars():
            allfilestr += """
void jacobianParam(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}
""" #+ "\n/* Variable and parameter substitutions undefined:*/\n" + parundefines + varundefines + "\n"
        # write out C file
        vffile = os.path.join(self._compilation_tempdir, self._vf_file)
        try:
            file = open(vffile, 'w')
            file.write(allfilestr)
            file.close()
        except IOError, e:
            print "Error opening file "+self._vf_file+" for writing"
            raise IOError, e


    def compileLib(self, libsources=[], libdirs=[]):
        """compileLib generates a python extension DLL with integrator and vector
        field compiled and linked.

        libsources list allows additional library sources to be linked.
        libdirs list allows additional directories to be searched for
          precompiled libraries."""

        if os.path.isfile(os.path.join(os.getcwd(),
                                "_radau5"+self._vf_filename_ext+self._dllext)):
            # then DLL file already exists and we can't overwrite it at this
            # time
            proceed = False
            print "\n"
            print "-----------------------------------------------------------"
            print "Present limitation of Python: Cannot rebuild library"
            print "without exiting Python and deleting the shared library"
            print "   " + str(os.path.join(os.getcwd(),
                                "_radau5"+self._vf_filename_ext+self._dllext))
            print "by hand! If you made any changes to the system you should"
            print "not proceed with running the integrator until you quit"
            print "and rebuild."
            print "-----------------------------------------------------------"
            print "\n"
        else:
            proceed = True
        if not proceed:
            print "Did not compile shared library."
            return
        if self._solver is not None:
            self.forceLibRefresh()
        vffile = os.path.join(self._compilation_tempdir, self._vf_file)
        try:
            ifacefile_orig = open(os.path.join(self._compilation_sourcedir,
                                               "radau5.i"), 'r')
            ifacefile_copy = open(os.path.join(self._compilation_tempdir,
                                       "radau5_"+self._vf_file[:-2]+".i"), 'w')
            firstline = ifacefile_orig.readline()
            ifacefile_copy.write('%module radau5_'+self._vf_file[:-2]+'\n')
            iffilestr = ifacefile_orig.read()
            ifacefile_copy.write(iffilestr)
            ifacefile_orig.close()
            ifacefile_copy.close()
        except IOError:
            print "radau5.i copying error in radau compilation directory"
            raise
        swigfile = os.path.join(self._compilation_tempdir,
                                "radau5"+self._vf_filename_ext+".i")
        radwrapfile = os.path.join(self._compilation_sourcedir, "radau5mod.c")
        if self._const_massmat:
            radfile = os.path.join(self._compilation_sourcedir, "radau5.f")
        else:
            # version for non-constant mass matrices
            radfile = os.path.join(self._compilation_sourcedir, "radau5v.f")
        dc_lapackfile = os.path.join(self._compilation_sourcedir, "dc_lapack.f")
        lapackfile = os.path.join(self._compilation_sourcedir, "lapack.f")
        lapackcfile = os.path.join(self._compilation_sourcedir, "lapackc.f")
        integfile = os.path.join(self._compilation_sourcedir, "integration.c")
        interfacefile = os.path.join(self._compilation_sourcedir, "interface.c")
        eventfile = os.path.join(self._compilation_sourcedir, "eventFinding.c")
        memfile = os.path.join(self._compilation_sourcedir, "memory.c")
        # The following if statement attempts to avoid recompiling the SWIG wrapper
        # if the files mentioned already exist, because in principle the SWIG interface
        # only needs compiling once. But this step doesn't seem to work yet.
        # Instead, it seems that SWIG always gets recompiled with everything else
        # (at least on Win32). Maybe the list of files is incorrect...
        if not (all([os.path.isfile(os.path.join(self._compilation_tempdir,
                   sf)) for sf in ['radau5'+self._vf_filename_ext+'_wrap.o',
                                   'lib_radau5'+self._vf_filename_ext+'.a',
                                   'radau5'+self._vf_filename_ext+'.py',
                                   '_radau5'+self._vf_filename_ext+'.def']])):
            modfilelist = [swigfile]
        else:
            modfilelist = []
##        modfilelist = [swigfile]
        modfilelist.extend([radwrapfile, vffile, integfile, eventfile,
                           interfacefile, memfile])
        modfilelist.extend(libsources)
        fortfilelist=[radfile, lapackcfile, dc_lapackfile, lapackfile]
        # would like temp dir to be '/' but scipy_distutils cannot build the
        # radau library in that configuration (unlike Dopri's identical use of
        # regular python distutils, which cannot deal with Fortran sources).
        script_args = ['-q', 'build', '--build-lib='+os.getcwd(),
                 '-t'+self._compilation_tempdir,
                 '--build-base='+self._compilation_tempdir,
                 '--build-platlib='+os.getcwd(), '--build-purelib='+os.getcwd()]
        if self._compiler != '':
            script_args.append('-c'+str(self._compiler))
        # include directories for libraries
        narraydir = os.path.join(get_python_inc(plat_specific=1), "numarray")
        incdirs = [narraydir, os.getcwd(), self._compilation_tempdir,
                   self._compilation_sourcedir]
        incdirs.extend(libdirs)
        radlibdirs = [self._compilation_sourcedir, self._compilation_tempdir,
                       os.getcwd()]
        # Use distutils to perform the compilation of the selected files
#        scipy_distutils.log.set_verbosity(4)  # doesn't work!
        rout.start()    # redirect stdout
        try:
            distobject = setup(name = "Radau5 integrator",
                  author = "PyDSTool (automatically generated)",
                  script_args = script_args,
                  ext_modules = [Extension("_radau5"+self._vf_filename_ext,
                                 sources=modfilelist,
                                 include_dirs=incdirs,
                                 library_dirs=radlibdirs+['./'],
                                 libraries=['radau5'],
                                 extra_compile_args=['-w', '-D__RADAU__'],
                                 extra_link_args=['-w'])],
#                  library_dirs=radlibdirs,
                  libraries=[('radau5',{'sources': fortfilelist,
                              'library_dirs': radlibdirs+['./']})])
        except:
            print "\nError occurred in generating Radau system..."
            print sys.exc_info()[0], sys.exc_info()[1]
            raise RuntimeError
        rout.stop()    # restore stdout
        # Attempt to unload module through a shutdown() function being
        # added to the SWIG module file. But it didn't work!
##        try:
##            modfilepy = open(os.path.join(self._compilation_tempdir,
##                                    "radau5"+self._vf_filename_ext+".py"), 'a')
##            extfilename = "_radau5"+self._vf_filename_ext
##            modfilepy.write("""# The following addition made by PyDSTool:
##def shutdown():
##    import sys
##    del sys.modules['""" + extfilename + """']
##            """)
####    del """ + extfilename + """
####    del new_doubleArray
####    del delete_doubleArray
####    del doubleArray_getitem
####    del doubleArray_setitem
####    del new_intArray
####    del delete_intArray
####    del intArray_getitem
####    del intArray_setitem
####    del Integrate
##            modfilepy.close()
##        except IOError:
##            print "radau5.py modifying error in radau temp compilation " \
##                  + "directory"
##            raise
        try:
            # move library files into the user's CWD
            if swigfile in modfilelist or not \
                   os.path.isfile(os.path.join(self._compilation_tempdir,
                                "radau5"+self._vf_filename_ext+".py")):
                # temporary hack to fix scipy_distutils bug
                shutil.move(os.path.join(os.getcwd(),
                                     self._compilation_tempdir,'src',
                                     "radau5"+self._vf_filename_ext+".py"),
                                 os.path.join(os.getcwd(),
                             "radau5"+self._vf_filename_ext+".py"))
        except:
            print "\nError occurred in generating Radau system"
            print "(while moving library extension modules to CWD)"
            print sys.exc_info()[0], sys.exc_info()[1]
            raise RuntimeError

    # doesn't work!
#    def _ensureLoaded(self, modname):
##        if modname in sys.modules:
##            _integMod = reload(sys.modules[modname])
##        else:
#        try:
#            _integMod = __import__(modname, globals())
#        except:
#            print "Error in importing compiled vector field and integrator."
#            print "Did you compile the RHS C code?"
#            raise
#        # Initialize integrator
#        assert 'Integrate' in dir(_integMod), \
#               "radau library does not contain Integrate()"
#        return _integMod


    def compute(self, trajname, dirn='f', ics=None):
        continue_integ = ODEsystem.prepDirection(self, dirn)
#        _integMod = __import__("radau5"+self._vf_filename_ext, globals())
        # validate spec if there exists a prior trajectory computation
        if self.defined:
            self.validateSpec()
            self.validateICs()
            self.clearWarnings()
            self.clearErrors()
        if ics is not None:
            self.set(ics=ics)
        if self.inputs and self._extInputsChanged:
            # self._extInputsChanged if globa t0 changed so that can
            # adjust times given to the integrator (it is blind to global t0
            # when accesses input variable times)
            self._inputVarList = []
            self._inputTimeList = []
            # inputVarList is a list of Variables or Pointsets
            for inp in sortedDictValues(self.inputs):
                if isinstance(inp, Variable):
                    pts = inp.getDataPoints()
                    if pts is None:
                        raise TypeError("Can only pass external input Variable objects if based on"
                                        " an underlying mesh")
                    else:
                        tvals = copy(pts[inp.indepvarname])
                        tvals -= self.globalt0
                    self._inputVarList.append(pts[inp.coordname].tolist())
                    self._inputTimeList.append(tvals.tolist())
                elif isinstance(inp, Pointset):
                    tvals = copy(inp.indepvararray)
                    tvals -= self.globalt0
                    self._inputVarList.append(inp[inp.coordname].tolist())
                    self._inputTimeList.append(tvals.tolist())
                else:
                    raise TypeError, "Invalid type of input"
        if isinstance(self._algparams['rtol'], list):
            if len(self._algparams['rtol']) != self.dimension:
                raise ValueError, 'rtol list must have same length as phase dimension'
        else:
            rtol = self._algparams['rtol']
            self._algparams['rtol'] = [rtol for i in range(self.dimension)]
        if isinstance(self._algparams['atol'], list):
            if len(self._algparams['atol']) != self.dimension:
                raise ValueError, 'atol list must have same length as phase dimension'
        else:
            atol = self._algparams['atol']
            self._algparams['atol'] = [atol for i in range(self.dimension)]
        anames = self.funcspec.auxvars
        # Check i.c.'s are well defined (finite)
        self.checkInitialConditions()
        self.setEventICs(self.initialconditions, self.globalt0)
        # update active & terminal events in case changed since last run
        self._algparams['eventActive'] = [int(self.eventstruct.events[evname].activeFlag) for evname in self._eventNames]
        self._algparams['eventTerm'] = [int(self.eventstruct.events[evname].termFlag) for evname in self._eventNames]
        self._algparams['eventDir'] = [int(self.eventstruct.events[evname].dircode) for evname in self._eventNames]
        # Main integration
        t0 = self.indepvariable.depdomain.get(0)
        t1 = self.indepvariable.depdomain.get(1)
        plist = sortedDictValues(self.pars)
        if self._solver is None:
#            _integMod = self._ensureLoaded("radau5"+self._vf_filename_ext)
            self._solver = radau("radau5"+self._vf_filename_ext,
                                 rhs=self.name, phaseDim=self.dimension,
                                 paramDim=len(plist), nAux=len(anames),
                                 nEvents=len(self._eventNames),
                                 nExtInputs=len(self.inputs),
                                 hasJac=self.haveJacobian(),
                                 hasJacP=self.haveJacobian_pars(),
                                 hasMass=self.haveMass(),
                                 extraSpace=self._algparams['extraspace'])
            try:
                genDB.register(self)
            except PyDSTool_KeyError:
                errstr = "Generator " + self.name + ": this vector field's " +\
                         "DLL is already in use"
                raise RuntimeError, errstr
        if self._dircode == 1:
            tbegin = t0
            tend = t1
        elif self._dircode == -1:
            # radau does reverse time integration simply by switching t0 and t1
            tbegin = t1
            tend = t0
        if len(self._algparams['specialtimes'])>0:
            use_special = self._algparams['use_special']
        else:
            use_special = 0
        bounds = [[],[]]  # lower, then upper
        for v in self.funcspec.vars:
            bds = self._xdomain[v]
            try:
                bounds[0].append(bds[0])
                bounds[1].append(bds[1])
            except TypeError:
                print v, type(bds), bds
                print self._xdomain
                raise
        for p in self.funcspec.pars:
            bds = self._pdomain[p]
            try:
                bounds[0].append(bds[0])
                bounds[1].append(bds[1])
            except TypeError:
                print type(bds), bds
                raise
        if continue_integ:
            x0 = self._solver.lastPoint
            # overwrite t0 from self.indepvariable.domain, but use its t1
            tbegin = self._solver.lastTime
            self._algparams['init_step'] = self._solver.lastStep
            if abs(t1-tbegin) < self._algparams['init_step']:
                raise ValueError("Integration end point too close to initial "
                                 "point")
#            if self.inputs and self._extInputsChanged:
#                self._extInputsChanged = False
#                self._solver.setContParams(tend, plist,
#                                           use_special,
#                                           self._algparams['verbose'],
#                                           True, deepcopy(self._inputVarList),
#                                           deeppcopy(self._inputTimeList))
        else:
            if self._solver.numRuns > 0:
                self._solver.clearAll()
            x0 = sortedDictValues(self.initialconditions, self.funcspec.vars)
            self._solver.setInteg(maxpts=self._algparams['max_pts'],
                rtol=self._algparams['rtol'], atol=self._algparams['atol'])
            self._solver.setRunParams(ic=x0, params=plist,
                                  t0=tbegin, tend=tend, gt0=self.globalt0,
                                  refine=self._algparams['refine'],
                                  specTimes=self._algparams['specialtimes'],
                                  bounds=bounds)
            if self.inputs:
                self._extInputsChanged = False
                self._solver.setExtInputs(True,
                                      deepcopy(self._inputVarList),
                                      deepcopy(self._inputTimeList))
        # hinit only set if not continue_integ
        if len(anames)>0:
            check_aux = self._algparams['check_aux']
        else:
            check_aux = 0
        if self._algparams['max_step'] == 0:
            max_step = abs(tend-tbegin)
        else:
            max_step = self._algparams['max_step']
        if continue_integ:
            alltData, X, A, Stats, H, Err, Evtimes, \
                Evpoints = self._solver.Continue(tend, plist,
                                  use_special, self._algparams['verbose'],
                                  self._extInputsChanged,
                                  deepcopy(self._inputVarList),
                                  deepcopy(self._inputTimeList), bounds)
            self._extInputsChanged = False
        else:
            self._solver.setEvents(eventActive=self._algparams['eventActive'],
                            eventTerm=self._algparams['eventTerm'],
                            eventDir=self._algparams['eventDir'],
                            eventDelay=self._algparams['eventDelay'],
                            eventInt=self._algparams['eventInt'],
                            eventTol=self._algparams['eventTol'],
                            maxevtpts=self._algparams['maxevtpts'],
                            maxbisect=self._algparams['maxbisect'])
            alltData, X, A, Stats, H, Err, Evtimes, \
                Evpoints = self._solver.Run(self._algparams['init_step'],
                                    max_step,
                                    check_aux,
                                    use_special,
                                    self._algparams['verbose'],
                                    self._algparams['safety'],
                                    self._algparams['jac_recompute'],
                                    self._algparams['newton_stop'],
                                    self._algparams['fac1'],
                                    self._algparams['fac2'],
                                    self._algparams['stepLB'],
                                    self._algparams['stepUB'],
                                    self._algparams['hessenberg'],
                                    self._algparams['max_newton'],
                                    self._algparams['newton_start'],
                                    self._algparams['index1dim'],
                                    self._algparams['index2dim'],
                                    self._algparams['index3dim'],
                                    self._algparams['step_strategy'],
                                    self._algparams['DAEstructureM1'],
                                    self._algparams['DAEstructureM2'],
                                    self.haveJacobian(),
                                    self.haveMass())
        self.outputStats = {'last_step': H,
                            'last_time': self._solver.lastTime,
                            'last_point': self._solver.lastPoint,
                            'num_fcns': Stats[0],
                            'num_jacs': Stats[1],
                            'num_steps': Stats[2],
                            'num_accept': Stats[3],
                            'num_reject': Stats[4],
                            'num_dec': Stats[5],
                            'num_subs': Stats[6],
                            'errorStatus': Err
                            }
        if self._dircode == -1:
            # reverse the numarray object (no reverse method!)
            alltData = alltData[::-1]
            X = X[:,::-1]
            if anames != []:
                A = A[:,::-1]
        xnames = self._var_ixmap
        try:
            allxDataDict = dict(zip(xnames,X))
        except IndexError:
            print "Integration returned variable values of unexpected dimensions."
            print "  Did you change the system and not refresh the C library" \
                  + " using the forcelibrefresh() method?"
            raise
        # storage of all auxiliary variable data
        anames = self.funcspec.auxvars
        try:
            if anames != []:
                allaDataDict = dict(zip(anames,A))
        except IndexError:
            print "Integration returned auxiliary values of unexpected dimensions."
            print "  Did you change the system and not refresh the C library" \
                  + " using the forcelibrefresh() method?"
            raise
        # Package up computed trajectory in Variable variables
        # Add external inputs warnings to self.warnings, if any
##        for f in inputVarList:
##            for winfo in f.warnings:
##                self.warnings.append((W_NONTERMSTATEBD,
##                                     (winfo[0], f.name, winfo[1],
##                                      f.depdomain)))
        eventslist = self.eventstruct.query(['lowlevel', 'active'])
        termevents = self.eventstruct.query(['term'], eventslist)
        if self._eventNames != []:
            # build self.warnings because events happened --
            # and keep a record of which times terminal events happened because
            # Model.py's event handling procedure assumes multiple events
            # happening at one time are listed in one warning
            termevtimes = {}
            nontermevtimes = {}
            try:
                for evix in range(len(self._eventNames)):
                    if Evpoints[evix] is None:
                        continue
                    numevs = len(Evtimes[evix])
                    if self._algparams['eventTerm'][evix]:
                        if numevs > 1:
                            print "Event info:", Evpoints, Evtimes
                        assert numevs <= 1, ("Internal error: more than one "
                                         "terminal event of same type found")
                        # For safety, we should assert that this event
                        # also appears in termevents, but we don't
                        evname = self._eventNames[evix]
                        if Evtimes[evix][0] in termevtimes.keys():
                            # append event name to this warning
                            warning_ix = termevtimes[Evtimes[evix][0]]
                            self.warnings[warning_ix][1][1].append(evname)
                        else:
                            # make new termevtime entry for the new warning
                            termevtimes[Evtimes[evix][0]] = len(self.warnings)
                            self.warnings.append((W_TERMEVENT,
                                             (Evtimes[evix][0],
                                             [self._eventNames[evix]])))
                    else:
                        for ev in range(numevs):
                            if Evtimes[evix][ev] in nontermevtimes.keys():
                                # append event name to this warning
                                warning_ix = nontermevtimes[Evtimes[evix][ev]]
                                self.warnings[warning_ix][1][1].append(evname)
                            else:
                                # make new nontermevtime entry for the new warning
                                nontermevtimes[Evtimes[evix][ev]] = \
                                                            len(self.warnings)
                                self.warnings.append((W_NONTERMEVENT,
                                                 (Evtimes[evix][ev],
                                                  [self._eventNames[evix]])))
            except IndexError:
                print "Events returned from integrator are the wrong size."
                print "  Did you change the system and not refresh the C " \
                      + "library using the forcelibrefresh() method?"
                raise
        termcount = 0
        for (w,i) in self.warnings:
            if w == W_TERMEVENT or w == W_TERMSTATEBD:
                if termcount > 0:
                    raise ValueError, \
                            "Internal error: more than one terminal event found"
                termcount += 1
        # Create variables (self.variables contains no actual data)
        variables = copyVarDict(self.variables)
        # build event pointset information (reset previous trajectory's)
        self.trajevents = {}
        for evix in range(len(self._eventNames)):
            evname = self._eventNames[evix]
            if Evpoints[evix] is None:
                self.trajevents[evname] = None
            else:
                self.trajevents[evname] = Pointset({'coordnames': xnames,
                                               'indepvarname': 't',
                                               'coordarray': Evpoints[evix],
                                               'indepvararray': Evtimes[evix]})
        if int(Err) > 0:
            for x in xnames:
                if len(alltData) > 1:
                    variables[x] = Variable(interp1d(alltData,allxDataDict[x]),
                                             't', x, x)
                else:
                    raise PyDSTool_ValueError, "Fewer than 2 data points computed"
            for a in anames:
                if len(alltData) > 1:
                    variables[a] = Variable(interp1d(alltData,allaDataDict[a]),
                                             't', a, a)
                else:
                    raise PyDSTool_ValueError, "Fewer than 2 data points computed"
            # final checks
            self.validateSpec()
            self.defined = True
            return Trajectory(trajname, variables.values(),
                              self.globalt0, self.checklevel)
        else:
            print 'Integration failed'
            try:
                self.errors.append((E_COMPUTFAIL, (self._solver.lastTime,
                                              self._errorcodes[int(Err)])))
            except TypeError:
                # errcode messed up from Radau
                print "Error information: ", Err
                self.errors.append((E_COMPUTFAIL, (self._solver.lastTime,
                                              self._errorcodes[0])))
            if self._solver.verbose:
                info(self.outputStats, "Output statistics")
            self.defined = False
            if self.outputStats['num_steps'] > self._algparams['max_pts']:
                print "max_pts algorithmic parameter too small: current value is %i"%self._algparams['max_pts']
#                avstep = (self._algparams['init_step']+self.outputStats['last_step'])/2.
                if self.outputStats['last_time']-tbegin > 0:
                    ms = int(round(self._algparams['max_pts']/(self.outputStats['last_time']-tbegin)*(tend-tbegin)))
                else:
                    ms = Inf
                print "(recommended value for this trajectory segment is estimated to be %s)"%str(ms)
            raise PyDSTool_ExistError, "No trajectory created"


    computeTraj = compute


    def Rhs(self, t, xdict, pdict, idict={}):
        x = sortedDictValues(xdict)
        p = sortedDictValues(pdict)
        if self._solver is None:
            x0=sortedDictValues(self.initialconditions)
            self._solver = radau("radau5"+self._vf_filename_ext,
                                 rhs=self.name, phaseDim=self.dimension,
                                 paramDim=len(p), nAux=len(self.funcspec.auxvars),
                                 nEvents=len(self._eventNames),
                                 hasJac=self.haveJacobian(),
                                 hasJacP=self.haveJacobian_pars(),
                                 hasMass=self.haveMass(),
                                 extraSpace=self._algparams['extraspace'])
            # tend value unimportant here
            self._solver.setRunParams(ic=x0, params=p,
                                  t0=0, tend=1, gt0=self.globalt0,
                                  refine=0, specTimes=[])
        return self._solver.Rhs(t, x, p)[0]


    def Jacobian(self, t, xdict, pdict, idict={}):
        x = sortedDictValues(xdict)
        p = sortedDictValues(pdict)
        if self._solver is None:
            x0=sortedDictValues(self.initialconditions)
            self._solver = radau("radau5"+self._vf_filename_ext,
                                 rhs=self.name, phaseDim=self.dimension,
                                 paramDim=len(p), nAux=len(self.funcspec.auxvars),
                                 nEvents=len(self._eventNames),
                                 hasJac=self.haveJacobian(),
                                 hasJacP=self.haveJacobian_pars(),
                                 hasMass=self.haveMass(),
                                 extraSpace=self._algparams['extraspace'])
            # tend value unimportant here
            self._solver.setRunParams(ic=x0, params=p,
                                  t0=0, tend=1, gt0=self.globalt0,
                                  refine=0, specTimes=[])
        return self._solver.Jacobian(t, x, p)[0]


    def JacobianP(self, t, xdict, pdict, idict={}):
        x = sortedDictValues(xdict)
        p = sortedDictValues(pdict)
        if self._solver is None:
            x0=sortedDictValues(self.initialconditions)
            self._solver = radau("radau5"+self._vf_filename_ext,
                                 rhs=self.name, phaseDim=self.dimension,
                                 paramDim=len(p), nAux=len(self.funcspec.auxvars),
                                 nEvents=len(self._eventNames),
                                 hasJac=self.haveJacobian(),
                                 hasJacP=self.haveJacobian_pars(),
                                 hasMass=self.haveMass(),
                                 extraSpace=self._algparams['extraspace'])
            # tend value unimportant here
            self._solver.setRunParams(ic=x0, params=p,
                                  t0=0, tend=1, gt0=self.globalt0,
                                  refine=0, specTimes=[])
        return self._solver.JacobianP(t, x, p)[0]


    def AuxVars(self, t, xdict, pdict, idict={}):
        x = sortedDictValues(xdict)
        p = sortedDictValues(pdict)
        if self._solver is None:
            x0=sortedDictValues(self.initialconditions)
            self._solver = radau("radau5"+self._vf_filename_ext,
                                 rhs=self.name, phaseDim=self.dimension,
                                 paramDim=len(p), nAux=len(self.funcspec.auxvars),
                                 nEvents=len(self._eventNames),
                                 hasJac=self.haveJacobian(),
                                 hasJacP=self.haveJacobian_pars(),
                                 hasMass=self.haveMass(),
                                 extraSpace=self._algparams['extraspace'])
            # tend value unimportant here
            self._solver.setRunParams(ic=x0, params=p,
                                  t0=0, tend=1, gt0=self.globalt0,
                                  refine=0, specTimes=[])
        return self._solver.AuxFunc(t, x, p)[0]


    def MassMatrix(self, t, xdict, pdict, idict={}):
        x = sortedDictValues(xdict)
        p = sortedDictValues(pdict)
        if self._solver is None:
            x0 = sortedDictValues(self.initialconditions)
            self._solver = radau("radau5"+self._vf_filename_ext,
                                 rhs=self.name, phaseDim=self.dimension,
                                 paramDim=len(p), nAux=len(self.funcspec.auxvars),
                                 nEvents=len(self._eventNames),
                                 hasJac=self.haveJacobian(),
                                 hasJacP=self.haveJacobian_pars(),
                                 hasMass=self.haveMass(),
                                 extraSpace=self._algparams['extraspace'])
            # tend value doesn't matter
            self._solver.setRunParams(ic=x0, params=p,
                                  t0=0, tend=1, gt0=self.globalt0,
                                  refine=0, specTimes=[])
        return self._solver.MassMatrix(t, x, p)[0]


    def __del__(self):
        genDB.unregister(self)
        ODEsystem.__del__(self)



# Register this Generator with the database

symbolMapDict = {'abs': 'fabs', 'sign': 'sgn'}
# in future, provide appropriate mappings for libraries math,
# random, etc. (for now it's left to FuncSpec)
theGenSpecHelper.add(Radau_ODEsystem, symbolMapDict, 'c')
