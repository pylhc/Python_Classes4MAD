# Basic integrator interface class
# Erik Sherwood, September 2006

from errors import PyDSTool_InitError as InitError
from scipy import isinf, Inf

# Need to add checks for when tolerances that should be nonneg are
# less than 0
class integrator:
    def __init__(self, rhs='default_name', phaseDim=0, paramDim=0, nAux=0,
                 nEvents=0, nExtInputs=0, hasJac=0, hasJacP=0, hasMass=0,
                 extraSpace=0, defaultBound=1e8):

        # Internal state variables
        self.initBasic = False
        self.initIntegrate = False
        self.initEvents = False
        self.initExtInputs = False
        self.setParams = False
        self.numRuns = 0
        self.numContinues = 0
        self.canContine = False

        # Run output
        self.times = []
        self.points = []
        self.auxPoints = []
        self.eventTimes = []
        self.eventPoints = []
        self.errors = []
        self.stats = []
        self.step = []

        self.lastPoint = []
        self.lastTime = []

        # Run parameters
        self.ic = []; self.params = []
        self.t0 = []; self.tend = []; self.gt0 = []
        self.refine = []
        self.upperBounds = []
        self.lowerBounds = []
        self.defaultBound = abs(float(defaultBound))

        # Event variables
        self.maxevtpts = []
        self.eventActive = []; self.eventDir = []; self.eventTerm = []
        self.eventInt = []; self.eventDelay=[]; self.eventTol=[]
        self.maxbisect = []

        # Integ variables
        self.maxpts = []; self.rtol = []; self.atol = []

        # Specific to integrator
        self.hinit = []; self.hmax = []; self.verbose = []
        self.specTimes = []
        self.direction = 0

        self.checkAux = 0
        self.calcSpecTimes = 0

        if self.checkBasic(rhs=rhs, phaseDim=phaseDim, paramDim=paramDim,
                           nAux=nAux, nEvents=nEvents, nExtInputs=nExtInputs,
                           hasJac=hasJac, hasJacP=hasJacP, hasMass=hasMass,
                           extraSpace=extraSpace):

            self.rhs = rhs
            self.phaseDim = phaseDim
            self.paramDim = paramDim
            self.extraSpace = extraSpace
            self.nAux = nAux
            self.nEvents = nEvents
            self.nExtInputs = nExtInputs

            if hasJac:
                self.hasJac = 1
            else:
                self.hasJac = 0

            if hasJacP:
                self.hasJacP = 1
            else:
                self.hasJacP = 0

            if hasMass:
                self.hasMass = 1
            else:
                self.hasMass = 0
        else:
            raise InitError, 'Integrator initialization failed!'

        # At this point, we expect the child class to set the self._integ field
        # and then call the initBasic method of the shared library.
        # The child class sets self.initBasic to True


    def __del__(self):
        try:
            self._integMod.CleanUp()
        except:
            pass


    def checkBasic(self, rhs, phaseDim, paramDim, nAux, nEvents, nExtInputs,
                   hasJac, hasJacP, hasMass, extraSpace):
        retval = True
        # Check that inputs to this function are correct
        if type(rhs) is not str:
            raise TypeError("rhs must be a string")

        if type(phaseDim) is not int:
            raise TypeError("phaseDim must be an integer")
        if phaseDim <= 0:
            raise ValueError("phaseDim must be positive")

        if type(paramDim) is not int:
            raise TypeError("paramDim must be an integer")
        if phaseDim < 0:
            raise ValueError("paramDim must be non-negative")

        if type(hasJac) is not int and type(hasJac) is not bool:
            raise TypeError("hasJac must be an integer or bool")
        if type(hasJac) is int and hasJac not in [0,1]:
            raise ValueError("integer hasJac must be 0 or 1");

        if type(hasJacP) is not int and type(hasJac) is not bool:
            raise TypeError("hasJacP must be an integer or bool")
        if type(hasJacP) is int and hasJacP not in [0,1]:
            raise ValueError("integer hasJacP must be 0 or 1");

        if type(hasMass) is not int and type(hasMass) is not bool:
            raise TypeError("hasMass must be an integer or bool")
        if type(hasMass) is int and hasMass not in [0,1]:
            raise ValueError("integer hasMass must be 0 or 1");

        if type(nAux) is not int:
            raise TypeError("nAux must be an integer")
        if nAux < 0:
            raise ValueError("nAux must be non-negative")

        if type(nEvents) is not int:
            raise TypeError("nEvents must be an integer")
        if nEvents < 0:
            raise ValueError("nEvents must be non-negative")

        if type(nExtInputs) is not int:
            raise TypeError("nExtInputs must be an integer")
        if nExtInputs < 0:
            raise ValueError("nExtInputs must be non-negative")

        if type(extraSpace) is not int:
            raise TypeError("extraSpace must be an integer")
        if extraSpace < 0:
            raise ValueError("extraSpace must be non-negative")

        return retval


    def setEvents(self, maxevtpts=1000, eventActive=[], eventDir=[], eventTerm=[],
                  eventInt=0.005, eventDelay=1e-4, eventTol=1e-6, maxbisect=100,
                  eventNearCoef=1000):

        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before setting events. (initBasic)'

        if self.initEvents == True:
            raise InitError, 'You must clear events before setting them. Use clearEvents()'

        # Currently we will not raise an error, but instead ignore setting events
        # if nEvents is zero. Just set to some default values and pass to the
        # shared library

        if self.nEvents <= 0:
            maxevtpts = 0
            eventActive = []
            eventDir = []
            eventTerm = []
            eventInt = 0.005
            eventDelay = 1e-4
            eventTol = 1e-6
            maxbisect = 100
            eventNearCoef=1000

        if self.checkEvents(maxevtpts, eventActive, eventDir, eventTerm,
                            eventInt, eventDelay, eventTol, maxbisect, eventNearCoef):
            self.maxevtpts = maxevtpts

            if type(eventActive) is int:
                evec = [eventActive for x in range(self.nEvents)]
                self.eventActive = evec
            else:
                self.eventActive = eventActive

            if type(eventDir) is int:
                evec = [eventDir for x in range(self.nEvents)]
                self.eventDir = evec
            else:
                self.eventDir = eventDir

            if type(eventTerm) is int:
                evec = [eventTerm for x in range(self.nEvents)]
                self.eventTerm = evec
            else:
                self.eventTerm = eventTerm

            if type(maxbisect) is int:
                evec = [maxbisect for x in range(self.nEvents)]
                self.maxbisect = evec
            else:
                self.maxbisect = maxbisect

            if type(eventInt) is float or type(eventInt) is int:
                evec = [float(eventInt) for x in range(self.nEvents)]
                self.eventInt = evec
            else:
                self.eventInt = eventInt

            if type(eventDelay) is float or type(eventDelay) is int:
                evec = [float(eventDelay) for x in range(self.nEvents)]
                self.eventDelay = evec
            else:
                self.eventDelay = eventDelay

            if type(eventTol) is float or type(eventTol) is int:
                evec = [float(eventTol) for x in range(self.nEvents)]
                self.eventTol = evec
            else:
                self.eventTol = eventTol

            self.eventNearCoef = eventNearCoef

            retval = self._integMod.InitEvents(self.maxevtpts, self.eventActive, self.eventDir,
                                               self.eventTerm, self.eventInt, self.eventDelay,
                                               self.eventTol, self.maxbisect, self.eventNearCoef)

            if retval[0] != 1:
                raise InitError, 'InitEvents call failed!'

            self.canContinue = False
            self.initEvents = True


    def clearEvents(self):
        if self.initBasic != True:
            raise ClearError, 'You must initialize the integrator before clearing events.'
        if self.initEvents == True:
            retval = self._integMod.ClearEvents()
            if retval[0] != 1:
                raise ClearError, 'ClearEvents call failed!'
            self.canContinue = False
            self.initEvents = False
        else:
            pass


    def checkEvents(self, maxevtpts, eventActive, eventDir, eventTerm,
                            eventInt, eventDelay, eventTol, maxbisect, eventNearCoef):
        retval = True

        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before checking events. (initBasic)'

        if type(maxevtpts) is not int:
            raise TypeError("maxevtpts must be integer.")
        if maxevtpts < 0:
            raise ValueError("maxevtpts must be non-negative.")

        if type(eventActive) not in [int, list]:
            raise TypeError("eventActive must be an integer or list.")
        if type(eventActive) is int:
            if eventActive not in [0, 1]:
                raise ValueError("integer eventActive must be 0 or 1.")
            else:
                pass
        if type(eventActive) is list:
            if len(eventActive) != self.nEvents:
                raise InitError, "list eventActive length must equal nEvents"
            for x in eventActive:
                if type(x) is not int:
                    raise TypeError("list eventActive entries must be integer")
                if x not in [0, 1]:
                    raise ValueError("list eventActive entries must be -1, 0 or 1")

        if type(eventDir) not in [int, list]:
            raise TypeError("eventDir must be an integer or list.")
        if type(eventDir) is int:
            if eventDir not in [-1, 0, 1]:
                raise ValueError("integer eventDir must be -1, 0 or 1.")
            else:
                pass
        if type(eventDir) is list:
            if len(eventDir) != self.nEvents:
                raise InitError, "list eventDir length must equal nEvents"
            for x in eventDir:
                if type(x) is not int:
                    raise TypeError("list eventDir entries must be integer")
                if x not in [-1, 0, 1]:
                    raise ValueError("list eventDir entries must be 0 or 1")

        if type(eventTerm) not in [int, list]:
            raise TypeError("eventTerm must be an integer or list.")
        if type(eventTerm) is int:
            if eventTerm not in [0, 1]:
                raise ValueError("integer eventTerm must be 0 or 1.")
            else:
                pass
        if type(eventTerm) is list:
            if len(eventTerm) != self.nEvents:
                raise InitError, "list eventTerm length must equal nEvents"
            for x in eventTerm:
                if type(x) is not int:
                    raise TypeError("list eventTerm entries must be integer")
                if x not in [0, 1]:
                    raise ValueError("list eventTerm entries must be 0 or 1")

        if type(eventInt) not in [int, float, list]:
            raise TypeError("eventInt must be an integer, float, or list.")
        if type(eventInt) in [int, float]:
            if eventInt < 0:
                raise ValueError("integer,float eventInt must be non-negative.")
            else:
                pass
        if type(eventInt) is list:
            if len(eventInt) != self.nEvents:
                raise InitError, "list eventInt length must equal nEvents"
            for x in eventInt:
                if type(x) not in [int, float]:
                    raise TypeError("list eventInt entries must be integer or float")
                if x < 0:
                    raise ValueError("list eventInt entries must be non-negative")

        if type(eventDelay) not in [int, float, list]:
            raise TypeError("eventDelay must be an integer, float, or list.")
        if type(eventDelay) in [int, float]:
            if eventDelay < 0:
                raise ValueError("integer,float eventDelay must be non-negative.")
            else:
                pass
        if type(eventDelay) is list:
            if len(eventDelay) != self.nEvents:
                raise InitError, "list eventDelay length must equal nEvents"
            for x in eventDelay:
                if type(x) not in [int, float]:
                    raise TypeError("list eventDelay entries must be integer or float")
                if x < 0:
                    raise ValueError("list eventDelay entries must be non-negative")

        if type(eventTol) not in [int, float, list]:
            raise TypeError("eventTol must be an integer, float, or list.")
        if type(eventTol) in [int, float]:
            if eventTol <= 0:
                raise ValueError("integer, float eventTol must be positive.")
            else:
                pass
        if type(eventTol) is list:
            if len(eventTol) != self.nEvents:
                raise InitError, "list eventTol length must equal nEvents"
            for x in eventTol:
                if type(x) not in [int, float]:
                    raise TypeError("list eventTol entries must be integer or float")
                if x <= 0:
                    raise ValueError("list eventTol entries must be positive")

        if type(maxbisect) not in [int, list]:
            raise TypeError("maxbisect must be an integer or list.")
        if type(maxbisect) is int:
            if maxbisect < 0:
                raise ValueError("integer maxbisect must be non-negative.")
            else:
                pass
        if type(maxbisect) is list:
            if len(maxbisect) != self.nEvents:
                raise InitError, "list maxbisect length must equal nEvents"
            for x in maxbisect:
                if type(x) is not int:
                    raise TypeError("list maxbisect entries must be integer")
                if x < 0:
                    raise ValueError("list maxbisect entries must be non-negative")

        if type(eventNearCoef) not in [int, float]:
            raise TypeError("eventNearCoef must be an integer or float")
        if type(eventNearCoef) in [int, float]:
            if eventNearCoef <= 0:
                raise ValueError("eventNearCoef must be positive.")

        return retval


    def setInteg(self, maxpts=100000, rtol=1e-9, atol=1e-12):
        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before setting integ. (initBasic)'

        if self.initIntegrate:
            raise InitError, 'You must clear integ before setting it. Use clearInteg()'

        if self.checkInteg(maxpts, rtol, atol):
            self.maxpts = maxpts

            if type(rtol) in [int, float]:
                tvec = [rtol for x in range(self.phaseDim)]
                self.rtol = tvec
            else:
                self.rtol = rtol

            if type(atol) in [int, float]:
                tvec = [atol for x in range(self.phaseDim)]
                self.atol = tvec
            else:
                self.atol = atol

            retval = self._integMod.InitInteg(self.maxpts, self.rtol, self.atol)

            if retval[0] != 1:
                raise InitError, 'InitInteg call failed!'

            self.canContinue = False
            self.initIntegrate = True


    def clearInteg(self):
        if not self.initBasic:
            raise ClearError, 'You must initialize the integrator before clearing events.'
        if self.initIntegrate:
            retval = self._integMod.ClearInteg()
            if retval[0] != 1:
                raise ClearError, 'ClearInteg call failed!'
            self.canContinue = False
            self.initIntegrate = False
        else:
            pass


    def checkInteg(self, maxpts, rtol, atol):
        retval = True

        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before checking integ. (initBasic)'

        if type(maxpts) is not int:
            raise TypeError("maxpts must be integer.")
        if maxpts <= 0:
            raise ValueError("maxpts must be positive")

        if type(rtol) not in [int, float, list]:
            raise TypeError("rtol must be an integer, float, or list.")
        if type(rtol) in [int, float]:
            if rtol <= 0:
                raise ValueError("integer, float rtol must be positive.")
            else:
                pass
        if type(rtol) is list:
            if len(rtol) != self.phaseDim:
                raise InitError, "list rtol length must equal phaseDim"
            for x in rtol:
                if type(x) not in [int, float]:
                    raise TypeError("list rtol entries must be integer or float")
                if x <= 0:
                    raise ValueError("list rtol entries must be positive")

        if type(atol) not in [int, float, list]:
            raise TypeError("atol must be an integer, float, or list.")
        if type(atol) in [int, float]:
            if atol <= 0:
                raise ValueError("integer, float atol must be positive.")
            else:
                pass
        if type(atol) is list:
            if len(atol) != self.phaseDim:
                raise InitError, "list atol length must equal phaseDim"
            for x in atol:
                if type(x) not in [int, float]:
                    raise TypeError("list atol entries must be integer or float")
                if x <= 0:
                    raise ValueError("list atol entries must be positive")

        return retval


    def setExtInputs(self, doCheck, extInputVals, extInputTimes):
        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before setting external Inputs. (initBasic)'

        if self.initExtInputs:
            raise InitError, 'You must clear extInputs before setting it. Use clearInteg()'

        if self.nExtInputs > 0:
            if doCheck:
                check = self.checkExtInputs(extInputVals, extInputTimes)
            else:
                check = True
            if check:
                self.extInputLens = []
                for x in range(self.nExtInputs):
                    self.extInputLens.append(len(extInputTimes[x]))

                IVals = extInputVals[0]
                ITimes = extInputTimes[0]
                for x in range(self.nExtInputs - 1):
                    IVals += extInputVals[x+1]
                    ITimes += extInputTimes[x+1]

                self.extInputVals = extInputVals
                self.extInputTimes = extInputTimes

            retval = self._integMod.InitExtInputs(self.nExtInputs,
                           self.extInputLens, IVals, ITimes)

            if retval[0] != 1:
                raise InitError, 'InitExtInputs call failed!'

            self.canContinue = False
            self.initExtInputs = True


    def clearExtInputs(self):
        if self.initBasic != True:
            raise ClearError, 'You must initialize the integrator before clearing external inputs.'
        if self.nExtInputs <= 0:
            self.extInputLens = []
            self.extInputVals = []
            self.extInputTimes = []
            self.initExtInputs = False
        elif self.initExtInputs:
            retval = self._integMod.ClearExtInputs()
            if retval[0] != 1:
                raise ClearError, 'ClearExtInputs call failed!'
            self.extInputLens = []
            self.extInputVals = []
            self.extInputTimes = []
            self.canContinue = False
            self.initExtInputs = False
        else:
            pass


    def checkExtInputs(self, extInputVals, extInputTimes):
        retval = True

        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before checking external inputs. (initBasic)'

        if type(extInputVals) is not list:
            raise TypeError("extInputVals must be list.")
        if len(extInputVals) != self.nExtInputs:
            raise ValueError("length of extInputVals must match nExtInputs")
        for x in range(self.nExtInputs):
            for y in range(len(extInputVals[x])):
                if type(extInputVals[x][y]) not in [int, float]:
                    raise TypeError("extInputVals entries must be int, float")

        for x in range(self.nExtInputs):
            for y in range(len(extInputTimes[x])):
                if type(extInputVals[x][y]) not in [int, float]:
                    raise TypeError("extInputTimes entries must be int, float")
            if len(extInputTimes[x]) > 1:
                # Check orientation
                orientation = extInputTimes[x][-1] - extInputTimes[x][0]
                if orientation == 0:
                    raise ValueError("Each extInputTime must be distinct; first and last times are identical with len > 1")
                if orientation < 0:
                    for y in range(len(extInputTimes[x])-1):
                        if extInputTimes[x][y] <= extInputTimes[x][y+1]:
                            raise ValueError("extInputTimes must be ordered consistently")
                if orientation > 0:
                    for y in range(len(extInputTimes[x])-1):
                        if extInputTimes[x][y] >= extInputTimes[x][y+1]:
                            print x, y, extInputTimes[x][y],extInputTimes[x][y+1]
                            raise ValueError("extInputTimes must be ordered consistently")

        return retval


    def setRunParams(self, ic=[], params=[], t0=[], tend=[], gt0=[], refine=0,
                     specTimes=[], bounds=[]):
        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before setting params. (initBasic)'

        #if self.initParams == True:
        #    raise InitError, 'You must clear params before setting them. Use clearParams()'

        if self.checkRunParams(ic, params, t0, tend, gt0, refine, specTimes,
                               bounds):
            self.ic = ic
            self.params = params
            self.t0 = float(t0)
            self.tend = float(tend)
            self.gt0 = float(gt0)
            self.refine = int(refine)
            self.specTimes = specTimes

            if self.t0 < self.tend:
                self.direction = 1
            else:
                self.direction = -1

            # Set bounds
            if bounds != []:
                self.upperBounds = bounds[1]
                self.lowerBounds = bounds[0]
                for i in range(self.phaseDim + self.paramDim):
                    if isinf(self.upperBounds[i]) and self.upperBounds[i] > 0:
                        self.upperBounds[i] = abs(float(self.defaultBound))
                    elif isinf(self.upperBounds[i]) and self.upperBounds[i] < 0:
                        self.upperBounds[i] = -abs(float(self.defaultBound))

                    if isinf(self.lowerBounds[i]) and self.lowerBounds[i] > 0:
                        self.lowerBounds[i] = abs(float(self.defaultBound))
                    elif isinf(self.lowerBounds[i]) and self.lowerBounds[i] < 0:
                        self.lowerBounds[i] = -abs(float(self.defaultBound))
            else:
                self.upperBounds = [abs(float(self.defaultBound)) for x in range(self.phaseDim + self.paramDim)]
                self.lowerBounds = [-abs(float(self.defaultBound)) for x in range(self.phaseDim + self.paramDim)]

            retval = self._integMod.SetRunParameters(self.ic, self.params,
                                 self.gt0, self.t0, self.tend, self.refine,
                                 len(self.specTimes), self.specTimes,
                                 self.upperBounds, self.lowerBounds)

            if retval[0] != 1:
                raise InitError, 'SetRunParameters call failed!'

            self.canContinue = False
            self.setParams = True


    def clearRunParams(self):
        if not self.initBasic:
            raise ClearError, 'You must initialize the integrator before clearing params.'
        if self.setParams:
            retval = self._integMod.ClearParams()
            if retval[0] != 1:
                raise ClearError, 'ClearParams call failed!'
            self.canContinue = False
            self.setParams = False
        else:
            pass


    def checkRunParams(self, ic, params, t0, tend, gt0, refine, specTimes,
                       bounds):
        retval = True

        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before checking run params. (initBasic)'

        if type(ic) is not list:
            raise TypeError("ic must be list")
        if len(ic) != self.phaseDim:
            raise ValueError, 'ic must have length equal to phaseDim'
        for x in ic:
            if type(x) not in [int, float]:
                raise TypeError("ic entries must be float, int")

        if type(params) is not list:
            raise TypeError("params must be list")
        if len(params) != self.paramDim:
            raise ValueError("params must have length equal to phaseDim")
        if len(params) > 0:
            for x in params:
                if type(x) not in [int, float]:
                    raise TypeError("params entries must be float, int")

        if type(refine) is not int:
            raise TypeError("refine must be int")
        if refine < 0:
            raise ValueError("refine must be non-negative")

        if type(t0) not in [int, float]:
            raise TypeError("t0 must be float, int")
        if type(tend) not in [int, float]:
            raise TypeError("tend must be float, int")
        if t0 == tend:
            raise ValueError("t0 must differ from tend")
        if t0 < tend:
            direction = 1
        else:
            direction = -1

        if type(specTimes) is not list:
            raise TypeError("specTimes must be list")
        if len(specTimes) > 0:
            if type(specTimes[0]) not in [int, float]:
                raise TypeError("specTimes entries must be int, float")
            if direction == 1:
                if specTimes[0] < t0 or specTimes[0] > tend:
                    raise ValueError("specTimes entries must be within [t0,tend]")
            else:
                if specTimes[0] > t0 or specTimes[0] < tend:
                    raise ValueError("specTimes entries must be within [tend,t0]")

        if len(specTimes) > 1:
            for x in range(len(specTimes)-1):
                if type(specTimes[x]) not in [int, float]:
                    raise TypeError("specTimes entries must be int, float")
            if direction == 1:
                if specTimes[x] < t0 or specTimes[x] > tend:
                    raise ValueError("specTimes entries must be within [t0,tend]")
                if specTimes[x] > specTimes[x+1]:
                    raise ValueError("specTimes must be non-decreasing")
            else:
                if specTimes[x] > t0 or specTimes[x] < tend:
                    raise ValueError("specTimes entries must be within [tend,t0]")
                if specTimes[x] < specTimes[x+1]:
                    raise ValueError("specTimes must be non-increasing")


        # Check the parameter, phase space bounds
        if type(bounds) is not list:
            raise TypeError("bounds must be list")
        if bounds != []:
            if len(bounds) != 2:
                raise TypeError("non-empty bounds must be a 2-list")
            for i in range(2):
                if type(bounds[i]) is not list:
                    raise TypeError("non-empty bounds must be a 2-list")
                if len(bounds[i]) != self.phaseDim + self.paramDim:
                    raise ValueError("bounds have incorrect size")
                for j in range(self.phaseDim + self.paramDim):
                    if type(bounds[i][j]) not in [int, float]:
                        raise TypeError("bounds entries must be int, float")

        return True


    def setContParams(self, tend, params, calcSpecTimes, verbose,
                      extInputChanged, extInputVals, extInputTimes, bounds):
        if self.direction > 0:
            if tend < self.tend:
                raise ValueError("new tend must be > old tend")
        if self.direction < 0:
            if tend > self.tend:
                raise ValueError("new tend must be < old tend")

        if type(params) is not list:
            raise TypeError("params must be list")
        if params != []:
            if len(params) != self.paramDim:
                raise ValueError("params must have length equal to phaseDim")
            if len(params) > 0:
                for x in params:
                    if type(x) not in [int, float]:
                        raise TypeError("params entries must be float, int")

        if type(calcSpecTimes) is not int:
            raise TypeError("calcSpecTimes must be 0 or 1")
        if calcSpecTimes not in [0,1]:
            raise TypeError("calcSpecTimes must be 0 or 1")
        if calcSpecTimes == 1 and len(self.specTimes) <= 0:
            raise ValueError("calcSpecTimes cannot be 1 if specTimes is empty")

        if type(verbose) is not int:
            raise TypeError("verbose must be 0 or 1")
        if verbose not in [0,1]:
            raise TypeError("verbose must be 0 or 1")

        if extInputChanged:
            if self.nExtInputs <= 0:
                raise ValueError("Cannot change extInputs if nExtInputs is 0")
            if self.initExtInputs:
                if self.checkExtInputs(extInputVals, extInputTimes):
                    self.clearExtInputs()
                    self.setExtInputs(False, extInputVals, extInputTimes)
            else:
                self.setExtInputs(True, extInputVals, extInputTimes)

        # Check and set the parameter, phase space bounds
        if type(bounds) is not list:
            raise TypeError("bounds must be list")
        if bounds != []:
            if len(bounds) != 2:
                raise TypeError("non-empty bounds must be a 2-list")
            for i in range(2):
                if type(bounds[i]) is not list:
                    raise TypeError("non-empty bounds must be a 2-list")
                if len(bounds[i]) != self.phaseDim + self.paramDim:
                    raise ValueError("bounds have incorrect size")
                for j in range(self.phaseDim + self.paramDim):
                    if type(bounds[i][j]) not in [int, float]:
                        raise TypeError("entries in bounds must be int, float")

            self.upperBounds = bounds[0]
            self.lowerBounds = bounds[1]
            for i in range(self.phaseDim + self.paramDim):
                if self.upperBounds[i] == Inf:
                    self.upperBounds[i] = abs(float(self.defaultBound))
                elif self.upperBounds[i] == -Inf:
                    self.upperBounds[i] = -abs(float(self.defaultBound))

                if self.lowerBounds[i] == Inf:
                    self.lowerBounds[i] = abs(float(self.defaultBound))
                elif self.lowerBounds[i] == -Inf:
                    self.lowerBounds[i] = -abs(float(self.defaultBound))

        # in case params blank, leave alone
        if params != []:
            self.params = params
        self.tend = tend

        retval = self._integMod.SetContParameters(self.tend, self.params,
                                         self.upperBounds, self.lowerBounds)
        if retval[0] != 1:
            raise InitError, 'Call to SetContParams failed!'

    def clearAll(self):
        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before clearing'
        self.clearRunParams()
        self.clearEvents()
        self.clearInteg()
        self.clearExtInputs()

    def Run(*args):
        if self.__class__==integrator:
            raise NotImplementedError, "Call Run on a concrete subclass"

    def Continue(*args):
        if self.__class__==integrator:
            raise NotImplementedError, "Call Continue on a concrete subclass"

    def Reset(self):
        if not self.initBasic:
            raise InitError, 'You must initialize the integrator before resetting'
        self._integMod.Reset()

        # What to do now?

    def Rhs(self, time, x, p):
        if self.initBasic:
            if self.nExtInputs > 0 and not self.initExtInputs:
                return None
            else:
                return self._integMod.Vfield(time, x, p)
        return None

    def Jacobian(self, time, x, p):
        if self.initBasic:
            if self.nExtInputs > 0 and not self.initExtInputs:
                return None
            if self.hasJac:
                return self._integMod.Jacobian(time, x, p)
        return None

    def JacobianP(self, time, x, p):
        if self.initBasic:
            if self.nExtInputs > 0 and not self.initExtInputs:
                return None
            if self.hasJacP:
                return self._integMod.JacobianP(time, x, p)
        return None

    def MassMatrix(self, time, x, p):
        if self.initBasic:
            if self.nExtInputs > 0 and not self.initExtInputs:
                return None
            if self.hasMass:
                return self._integMod.MassMatrix(time, x, p)
        return None

    def AuxFunc(self, time, x, p):
        if self.initBasic:
            if self.nExtInputs > 0 and not self.initExtInputs:
                return None
            if self.nAux > 0:
                return self._integMod.AuxFunc(time, x, p)
        return None

