#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      01 Sep 2006
#$Revision: 1.1 $
"""Model Constructor classes.

   Instantiate abstract Model Specifications into concrete simulation code.

   Robert Clewley, September 2005.


Overview of steps that ModelConstructor takes:

1. Build Generators for each v.f.
 * Take flattened spec and select a compatible Gen type and language
   * include check that compatibleGen selection actually exists,
      and is compatible with specType/domain, and with cts/discrete
      time variables)
   * presently, each v.f. must have all the same Gen type
 * Map things like x^y and x**y power specs -> pow(x,y)
    and abs() function -> fabs() in C.

2. Given the specific target lang/Gen with a ModelSpec v.f.
 * Associate additional events, code inserts, reused vars
 * Choose which domains map to domain verification events (certainly none that
    are infinite -- although semi-infinite can be checked in one direction)
 * All non-infinite domain ends would have an event created for them, but
    default is for them to be inactive.

 The 'mspec' argument to the GeneratorConstructor class must be a complete
 ModelSpec (up to introduction of global references including definitions of
 external inputs).

 The resulting model output from getModel contains the ModelSpec mspec in
 order to use its structure in resolving information about the relationship
 between variables.
"""

# ----------------------------------------------------------------------------
# VERSION HISTORY
#
# Version 1.1
#   Moved bounds event detection to Generator/baseclasses.py
#   Fixed bug where ModelConstructor's getModel() failed to pass on inputs
#     information in initialization of GeneratorConstructor
#   Fixed bug in makeGenInfoEntry to include names of terminal events in
#     'special reasons' list of reasons a simulationg could end. This fixes
#     use of 'embed' on systems with terminal events.
#
# Version 1.0.8
#   Fixed bug in EventMapping constructor.
#
# Version 1.0.7
#   Changed event mapping constructor to allow a dictionary of activation and
#    inactivations of events.
#
# Version 1.0.6
#   Fixed embed() to use generator's initial condition if none given.
#   Changed 'xdatadict' Generator __init__ key to 'ics'.
#
# Version 1.0.5
#   Made embed() copy the Generator argument.
#
# Version 1.0.4
#   Moved makeEvMapping() function from Toolbox.
#
# Version 1.0.3
#   Addition of EvMapping class.
#   Minor bug fixes.
#
# Version 1.0.2
#   Made embed() use the single Generator's name for the model by default.
#   Added usage of new simplify() method for QuantSpecs()
#
# Version 1.0.1
#   Change processReused() to make u_subsMap before getting into loops
#     so that all auxvars are always substituted in any auxvar-auxvar
#     inter-dependency.
#
# ----------------------------------------------------------------------------

# -----------------------------------------------------------------------------

# PyDSTool imports
from errors import *
from common import *
from utils import info, remain, intersect
import Model, Generator, ModelSpec, Symbolic, Events
from Generator import ctsGen, discGen, theGenSpecHelper
from parseUtils import symbolMapClass, NAMESEP, isNumericToken

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue, any, all
from numarray import array, arange, zeros, ones, Int, Float, Complex, Int32, \
     Float64, Complex64, concatenate, swapaxes, take
import sys, types, copy

# Exports
__all__ = ['GeneratorConstructor', 'ModelConstructor', 'makeGenInfo',
           'makeGenInfoEntry', 'embed', 'EvMapping', 'makeEvMapping']

# -----------------------------------------------------------------------------

mathNameMap = dict(zip(Symbolic.allmathnames_symbolic,Symbolic.allmathnames))

class GeneratorConstructor(object):
    def __init__(self, mspec=None, userevents=[], unravelInfo=True,
                 inputs={}, checklevel=2, activateAllBounds=False,
                 activatedBounds={},
                 targetGen="", algparams={}, indepvar=('t',[-Inf,Inf]),
                 parvalues={}, icvalues={}, reuseTerms={},
                 options={}, abseps=None, eventtol=1e-8):

        self.mspec = mspec
        # user events are additional to the intrinsic constraint events
        # that are made automatically from the variables' bounds information
        self.userevents = copy.copy(userevents)
        self.unravelInfo = unravelInfo   # presently just a Boolean
        if type(targetGen)==str:
            self.targetGen = targetGen
        else:
            raise TypeError, "targetGen argument must be a string"
        self.algparams = copy.copy(algparams)
        self.indepvarname = indepvar[0]
        self.indepvardomain = indepvar[1]
        self.inputs = copy.copy(inputs)
        self.parvalues = copy.copy(parvalues)
        self.icvalues = copy.copy(icvalues)
        self.checklevel = checklevel
        self.forcedAuxVars = []
        self.optDict = copy.copy(options)
        self.reuseTerms = copy.copy(reuseTerms)
        self.vfcodeinsert_start = ""
        self.vfcodeinsert_end = ""
        self.activatedBounds = copy.copy(activatedBounds)
        self.activateAllBounds = activateAllBounds  # overrides activatedBounds
        self.abseps = abseps
        self.eventtol = eventtol


    def setForcedAuxVars(self, vlist):
        self.forcedAuxVars = vlist

    def setReuseTerms(self, rdict):
        self.reuseTerms = rdict

    def setVfCodeInsertStart(self, codestr):
        self.vfcodeinsert_start = codestr

    def setVfCodeInsertEnd(self, codestr):
        self.vfcodeinsert_end = codestr

    def setOptions(self, optDict):
        # e.g. for 'nobuild' option for C-based ODE integrators
        self.optDict = optDict

    def addEvents(self, evtarg):
        if isinstance(evtarg, list):
            self.userevents.extend(evtarg)
        elif isinstance(evtarg, Event):
            self.userevents.append(evtarg)
        else:
            raise TypeError("Invalid event or event list")

    def activateBounds(self, varname=None, which_bounds='all'):
        """which_bounds argument is either 'lo', 'hi', or a pair ('lo', 'hi').
        Calling with no arguments activates all bounds."""
        if varname is None and which_bounds=='all':
            self.activateAllBounds = True
        else:
            entry = [False,False]
            if 'hi' in which_bounds:
                entry[1] = True
            if 'lo' in which_bounds:
                entry[0] = True
            self.activatedBounds[varname] = entry


    def getGenerator(self):
        """Build and return a Generator instance from an abstract
        specification."""
        ### Instantiate (flatten) target model structured specification
        ## Flatten ModelSpec self.mspec using indepvarname global and inputs
        # and using connectivity bindings (latter not yet implemented)
        globalRefs = [self.indepvarname] + self.inputs.keys()
        try:
            flatspec = self.mspec.flattenSpec(self.unravelInfo, globalRefs,
                                              ignoreInputs=True)
        except:
            print "Problem flattening Model Spec '%s'"%self.mspec.name
            print "Global refs: ", globalRefs
            raise
        FScompatibleNames = flatspec['FScompatibleNames']
        FScompatibleNamesInv = flatspec['FScompatibleNamesInv']
        ## Check target Generator info
        if self.targetGen in self.mspec.compatibleGens \
                   and self.targetGen in theGenSpecHelper:
            gsh = theGenSpecHelper(self.targetGen)
            if gsh.lang not in self.mspec.targetLangs:
                raise ValueError("Incompatible target language between supplied"
                                 " ModelSpec and target Generator")
        else:
            print "ModelSpec's compatible Generators:", \
                  ", ".join(self.mspec.compatibleGens)
            print "ModelConstructor target Generator:", self.targetGen
            raise ValueError("Target Generator mismatch during generator "
                             "construction")
        self.targetLang = gsh.lang
        ## Make Generator initialization argument dictionary
        a = args()
        if self.abseps is not None:
            a.abseps = self.abseps
        a.pars = {}
        parnames = flatspec['pars'].keys()
        for p, valstr in flatspec['pars'].iteritems():
            if valstr == '':
                if FScompatibleNamesInv(p) not in self.parvalues:
                    raise ValueError("Parameter %s is missing a value"%FScompatibleNamesInv(p))
            else:
                try:
                    a.pars[p] = float(valstr)
                except ValueError:
                    raise ValueError("Invalid parameter value set in ModelSpec"
                                     " for '%s'"%p)
        # override any par vals set in ModelSpec with those explicitly set
        # here
        for p, val in self.parvalues.iteritems():
            try:
                pr = FScompatibleNames(p)
            except KeyError:
                raise NameError("Parameter '%s' missing from ModelSpec"%p)
            if pr not in flatspec['pars']:
                raise NameError("Parameter '%s' missing from ModelSpec"%p)
            a.pars[pr] = val
        if self.icvalues != {}:
            a.ics = {}
            for v, val in self.icvalues.iteritems():
                try:
                    vr = FScompatibleNames(v)
                except KeyError:
                    raise NameError("Variable '%s' missing from ModelSpec"%v)
                if vr not in flatspec['vars']:
                    raise NameError("Variable '%s' missing from ModelSpec"%v)
                a.ics[vr] = val
        a.tdomain = self.indepvardomain
        # a.ttype = 'float' or 'int' ?
        a.inputs = self.inputs
        a.name = self.mspec.name
        xdomain = {}
        pdomain = {}
        for k, d in flatspec['domains'].iteritems():
            # e.g. d == (float, Continuous, [-Inf, Inf])
            assert d[1] == gsh.domain, "Domain mismatch with target Generator"
            if k in flatspec['vars']:
                assert len(d[2]) == 2, "Domain spec must be a valid interval"
                xdomain[k] = d[2]
            elif k in flatspec['pars']:
                assert len(d[2]) == 2, "Domain spec must be a valid interval"
                pdomain[k] = d[2]
        a.xdomain = xdomain
        a.pdomain = pdomain
        exp_vars = [v for (v,t) in flatspec['spectypes'].items() \
                         if t == 'ExpFuncSpec']
        rhs_vars = [v for (v,t) in flatspec['spectypes'].items() \
                         if t == 'RHSfuncSpec']
        imp_vars = [v for (v,t) in flatspec['spectypes'].items() \
                         if t == 'ImpFuncSpec']
        if gsh.specType == 'RHSfuncSpec':
            assert imp_vars == [], "Cannot use implicitly defined variables"
            assert self.forcedAuxVars == [], "Cannot force auxiliary variables"
            varnames = rhs_vars
            auxvarnames = exp_vars
        elif gsh.specType == 'ExpFuncSpec':
            assert imp_vars == [], "Cannot use implicitly defined variables"
            assert rhs_vars == [], "Cannot use RHS-type variables"
            varnames = exp_vars
            invalid_auxvars = remain(self.forcedAuxVars, varnames)
            if invalid_auxvars == []:
                # then all forced aux varnames were legitimate
                # so remove them from varnames and put them in auxvarnames
                varnames = remain(varnames, self.forcedAuxVars)
                auxvarnames = self.forcedAuxVars
            else:
                print "Invalid auxiliary variable names:"
                print invalid_auxvars
                raise ValueError("Forced auxiliary variable names were invalid")
        elif gsh.specType == 'ImpFuncSpec':
            assert rhs_vars == [], "Cannot use RHS-type variables"
            varnames = imp_vars
            auxvarnames = exp_vars
        # search for explicit variable interdependencies and resolve by
        # creating 'reuseterms' declarations, substituting in the cross-ref'd
        # definitions
        # e.g. state variables v and w, and explicit aux vars are given by:
        #     x = 1+v
        #     y = f(x) + w
        # Here, y illegally depends on x, so define a 'reused' temporary
        # definition, and re-write in terms of that:
        #     temp = 1+v
        #     x = temp
        #     y = f(temp) + w
        # e.g. state variables v and w, aux var x:
        #     v' = 1-v -f(x)
        # Here, v illegally uses an auxiliary variable on the RHS, so make
        # a 'reused' substitution as before
        #
        # first pass to find which substitutions are needed
        reuseTerms, subsExpr = processReused(varnames+auxvarnames, auxvarnames,
                                       flatspec, self.mspec._registry,
                                       FScompatibleNames, FScompatibleNamesInv)
        clash_reused = intersect(reuseTerms.keys(), self.reuseTerms.keys())
        if clash_reused != []:
            print "Clashing terms:", clash_reused
            raise ValueError("User-supplied reused terms clash with auto-"
                             "generated terms")
        # second pass, this time to actually make the substitutions
        for v in subsExpr:
            flatspec['vars'][v] = subsExpr[v]
        reuseTerms.update(self.reuseTerms)
        a.reuseterms = reuseTerms
        a.varspecs = dict(zip(varnames+auxvarnames, [flatspec['vars'][v] \
                                      for v in varnames+auxvarnames]))
        a.auxvars = auxvarnames
        try:
            a.fnspecs = flatspec['auxfns']
        except KeyError:
            # no aux fns defined!
            pass
        a.checklevel = self.checklevel
        a.algparams = self.algparams
        if self.vfcodeinsert_start != "":
            a.vfcodeinsert_start = self.vfcodeinsert_start
        if self.vfcodeinsert_end != "":
            a.vfcodeinsert_end = self.vfcodeinsert_end
        ## Events
        events = []
        # make events from bound constraints (activated accordingly)
        # (parameter bounds only useful for continuation with PyCont)
        domnames = varnames+parnames
        for xname in domnames:
            hier_name_lo = FScompatibleNamesInv(xname)+"_domlo"
            FScompatibleNames[hier_name_lo] = xname+"_domlo"
            FScompatibleNamesInv[xname+"_domlo"] = hier_name_lo
            hier_name_hi = FScompatibleNamesInv(xname)+"_domi"
            FScompatibleNames[hier_name_hi] = xname+"_domhi"
            FScompatibleNamesInv[xname+"_domhi"] = hier_name_hi
        if self.activateAllBounds:
            a.activatedbounds = {}.fromkeys(domnames,(True,True))
        else:
            a.activatedbounds = self.activatedBounds
        a.enforcebounds = True
        # add events from user events
        for e in self.userevents:
            if e not in events:
                events.append(e)
            else:
                raise ValueError, "Repeated event definition!"
#        events.extend(self.userevents)
        a.events = events
        # Add any additional special options (e.g. 'nobuild' directive)
        for k,v in self.optDict.iteritems():
            if hasattr(a, k):
                raise KeyError("'%s' already exists as a Generator argument"%k)
            a.k = v
        # keep a copy of the arguments in self for users to see what was done
        self.conargs = a
        self.FScompatibleNames = FScompatibleNames
        self.FScompatibleNamesInv = FScompatibleNamesInv
        ## Make Generator
        try:
            return eval('Generator.'+self.targetGen + "(a)")
        except:
            print "Problem initializing target Generator '%s'"%self.targetGen
            raise


    def __del__(self):
        del self.userevents
        if hasattr(self, 'conargs'):
            try:
                del self.conargs.events
            except AttributeError:
                pass


# -----------------------------------------------------------------------------


class ModelConstructor(object):
    def __init__(self, name, userevents={}, unravelInfo=True,
                 inputs={}, checklevel=2, activateAllBounds=False,
                 generatorspecs={}, indepvar=('t',[-Inf,Inf]),
                 parvalues={}, icvalues={}, reuseTerms={},
                 abseps=None, eventtol=1e-8):
        self.name = name
        self.forcedIntVars = []
        self._generators = generatorspecs
        self._events = copy.copy(userevents)
        self.indepvar = indepvar
        self.eventmaps = {}
        self.reuseTerms = copy.copy(reuseTerms)
        self.inputs = copy.copy(inputs)
        self.parvalues = copy.copy(parvalues)
        self.icvalues = copy.copy(icvalues)
        self.checklevel = checklevel
        self.unravelInfo = unravelInfo   # presently just a Boolean
        self.activatedBounds = {}
        self.activateAllBounds = activateAllBounds  # overrides activatedBounds
        self.abseps = abseps
        self.eventtol = eventtol

    def __repr__(self):
        return "ModelConstructor %s"%self.name

    def addGenInfo(self, genSpec, genTarg, genAlgPars={}, unravelInfo={},
                   genOpts={}):
        self._generators[genSpec.name] = {'modelspec': genSpec,
                                          'target': genTarg,
                                          'algparams': copy.copy(genAlgPars),
                                          'unravelinfo': copy.copy(unravelInfo),
                                          'options': copy.copy(genOpts)
                                          }        

    def getModel(self):
        """Build and return (hybrid) model made up of declared Generators and
        the mappings between events used to change vector fields in a hybrid
        system."""
        # 1. build constituent generators
        # 2. combine all generators' FScompatibleNames symbol maps
        FScompatibleNames = {}
        FScompatibleNamesInv = {}
        allGenNames = []
        genObjs = {}
        for gname, geninfodict in self._generators.iteritems():
            genSpec = geninfodict['modelspec']
            allGenNames.append(genSpec.name)
            assert gname == genSpec.name, "Generator name mismatch in geninfodict"
            genTarg = geninfodict['target']
            if 'algparams' in geninfodict:
                genAlgPars = geninfodict['algparams']
            else:
                genAlgPars = {}
            if self.inputs is not None:
                genInputs = self.inputs
            else:
                genInputs = {}
            if 'unravelinfo' in geninfodict:
                genUnravelInfo = geninfodict['unravelinfo']
            else:
                genUnravelInfo = True
            if 'options' in geninfodict:
                genOpts = geninfodict['options']
            else:
                genOpts = {}
            try:
                genEvents = self._events[gname]
            except KeyError:
                genEvents = []
            # extract par values and ic values relevant to this generator
            genPars = {}
            for p, val in self.parvalues.iteritems():
                # don't bother to check that p is a valid param name
                # for this generator -- that will be checked by
                # GeneratorConstructor
                if p in genSpec._registry:
                    genPars[p] = val
            genICs = {}
            for v, val in self.icvalues.iteritems():
                # don't bother to check that v is a valid variable name
                # for this generator -- that will be checked by
                # GeneratorConstructor
                if v in genSpec._registry:
                    genICs[v] = val
            genCon = GeneratorConstructor(genSpec, checklevel=self.checklevel,
                                          userevents=genEvents,
                                          targetGen=genTarg,
                                          algparams=genAlgPars,
                                          indepvar=self.indepvar,
                                          parvalues=genPars,
                                          inputs=genInputs,
                                          icvalues=genICs,
                                          options=genOpts,
                                          unravelInfo=genUnravelInfo,
                                          reuseTerms=self.reuseTerms,
                                          abseps=self.abseps,
                                          eventtol=self.eventtol,
                                          activateAllBounds=self.activateAllBounds,
                                          activatedBounds=self.activatedBounds)
            genObjs[gname] = genCon.getGenerator()
            # assume that there won't be any name clashes (there shouldn't be)
            FScompatibleNames.update(genCon.FScompatibleNames.lookupDict)
            FScompatibleNamesInv.update(genCon.FScompatibleNamesInv.lookupDict)
        # 3. build event mappings
        genInfoEntries = {}
        for hostGen, genObj in genObjs.iteritems():
            allGenTermEvents = genObj.eventstruct.getTermEvents()
            allGenTermEvNames = [e[0] for e in allGenTermEvents]
            if hostGen in self.eventmaps:
                genMaps = self.eventmaps[hostGen]
                genMapNames = []
                for gmtuple in genMaps:
                    genMapNames.append(gmtuple[0])
                for evname in remain(allGenTermEvNames, genMapNames):
                    genMaps.append((evname, 'terminate'))
                genInfoEntries[hostGen] = makeGenInfoEntry(genObj,
                                                           allGenNames,
                                                           genMaps)
            else:
                # default for a generator without an event mapping is to
                # terminate when its time is up.
                genMaps = [('time', 'terminate')]
                for evname in allGenTermEvNames:
                    genMaps.append((evname, 'terminate'))
                if not isfinite(genObj.indepvariable.depdomain.get(1)):
                    print "Warning: Generator %s has no termination event"%genObj.name
                    print "because it has an non-finite end computation time..."
                genInfoEntries[hostGen] = makeGenInfoEntry(genObj,
                                                           allGenNames,
                                                           genMaps)
        genInfoDict = makeGenInfo(genInfoEntries.values())
        # 4. build model
        mod_args = {'name': self.name,
               'genInfo': genInfoDict,
               'mspecdict': copy.copy(self._generators),
               'FScompatibleNames': symbolMapClass(FScompatibleNames),
               'FScompatibleNamesInv': symbolMapClass(FScompatibleNamesInv)}
        model = Model.Model(mod_args)
        model.forceIntVars(self.forcedIntVars)
        del genObjs
        del genInfoEntries
        del genInfoDict
        return model

    def addEvents(self, hostGen, evTarg):
        if hostGen not in self._events:
            self._events[hostGen] = []
        if isinstance(evTarg, list):
            self._events[hostGen].extend(evTarg)
        elif isinstance(evTarg, Events.Event):
            self._events[hostGen].append(evTarg)
        else:
            raise TypeError("Invalid event or event list")

    def setReuseTerms(self, rdict):
        self.reuseTerms = rdict

    def activateBounds(self, varname=None, which_bounds='all'):
        """which_bounds argument is either 'lo', 'hi', or a pair ('lo', 'hi').
        Calling with no arguments activates all bounds."""
        if varname is None and which_bounds=='all':
            self.activateAllBounds = True
        else:
            entry = [0,0]
            if 'hi' in which_bounds:
                entry[1] = 1
            if 'lo' in which_bounds:
                entry[0] = 1
            self.activatedBounds[varname] = entry

    def setInternalVars(self, arg):
        if isinstance(arg, list):
            self.forcedIntVars = arg
        elif isinstance(arg, str):
            self.forcedIntVars = [arg]
        # !! Should check whether these are valid variable names of model

    def mapEvent(self, hostGen, eventname, target, eventmapping=None):
        # have to have declared all generators first
        allGenNames = []
        for gname, geninfodict in self._generators.iteritems():
            allGenNames.append(geninfodict['modelspec'].name)
        if target not in allGenNames and target != 'terminate':
            raise ValueError("Unknown target Generator %s"%target)
        if hostGen not in self._generators:
            raise ValueError("Unknown host Generator %s"%hostGen)
        try:
            genEvs = self._events[hostGen]
        except KeyError:
            genEvs = []
        evNames = [ev.name for ev in genEvs]
        if eventname not in evNames and eventname != 'time':
            raise ValueError("Unknown event '%s' for host Generator"
                             " '%s'"%(eventname, hostGen))
        if hostGen in self.eventmaps:
            self.eventmaps[hostGen].append((eventname, (target, eventmapping)))
        else:
            self.eventmaps[hostGen] = [(eventname, (target, eventmapping))]


# ---------------------------------------------------------------------------

## Utility functions
def embed(gen, icdict=None, name=None):
    """Only use this function for building non-hybrid models with single
    Generators. Otherwise, use the ModelConstructor class.

    NB The supplied Generator is *copied* into the model."""
    assert compareBaseClass(gen, Generator.Generator), ("gen argument "
                                        "must be a Generator object")
    if name is None:
        name = gen.name + "_model"
    genInfoEntry = makeGenInfoEntry(copy.copy(gen), [gen.name])
    modelArgs = {'name': name,
                 'genInfo': makeGenInfo([genInfoEntry])}
    if icdict is None:
        modelArgs['ics'] = gen.initialconditions
    else:
        modelArgs['ics'] = icdict
    return Model.Model(modelArgs)


def makeGenInfo(arg):
    if len(arg) == 1 and isinstance(arg, dict):
        genList = [arg]
    else:
        genList = arg
    allGenNames = []
    returnDict = {}
    for infodict in genList:
        assert len(infodict) == 1, \
                   "Incorrect length of Generator info dictionary"
        genName = infodict.keys()[0]
        if genName not in allGenNames:
            allGenNames.append(genName)
            returnDict.update(infodict)
        else:
            raise ValueError("clashing Generator names in Generator info "
                             "dictionaries")
        assert isinstance(infodict.values()[0], tuple), \
                             'Expected tuple for info dict'
        assert len(infodict.values()[0]) == 2, \
                             'Incorrect size of Generator info tuple'
    return returnDict


class EvMapping(object):
    """Event mapping class, for use by makeGenInfoEntry and, when
    instantiated, the Model class.

    defStrings (list of statements) overrides assignDict if supplied at
    initialization, to permit full flexibility in the contents of the
    event mapping function."""

    def __init__(self, assignDict={}, defString="",
                 activeDict={}):
        self.assignDict = copy.copy(assignDict)
        self.defString = defString
        self.activeDict = copy.copy(activeDict)
        self.makeCallFn()


    def makeCallFn(self):
        indent = "  "
        fnString = """def evmapping(self, xdict, pdict, idict, estruct, t):"""
        if self.defString == "" and self.assignDict == {} and self.activeDict == {}:
            # default is the "identity mapping" (do nothing)
            fnString += indent + "pass\n"
        elif len(self.defString) >= 13 and self.defString[:13] == "def evmapping":
            # already defined, probably rebuilding after save/load object
            fnString = self.defString
        else:
            if len(self.assignDict) > 0:
                for lhs, rhs in self.assignDict.iteritems():
                    if not(type(lhs)==type(rhs)==str):
                        raise TypeError("Assignment dictionary for event "
                                        "mapping must consist of strings for "
                                        "both keys and values")
                fnString += "\n" + indent + ("\n"+indent).join(["%s = %s"%(l,r) \
                                    for l, r in self.assignDict.items()])
            if len(self.defString) > 0:
                fnString += "\n" + indent + ("\n"+indent).join(self.defString.split("\n"))
            if len(self.activeDict) > 0:
                for evname, state in self.activeDict.iteritems():
                    if not(type(evname)==str and type(state)==bool):
                        raise TypeError, "Invalid types given for setting active events"
                fnString += "\n" + indent + \
                         ("\n"+indent).join(["estruct.setActiveFlag('%s',%s)"%(evname,str(state)) \
                                    for evname, state in self.activeDict.items()])
            self.defString = fnString
        try:
            exec fnString
        except:
            print 'Invalid function definition for event mapping:'
            print fnString
            raise
        setattr(self, 'evmapping', types.MethodType(locals()['evmapping'],
                                                      self, self.__class__))

    def __getstate__(self):
        d = copy.copy(self.__dict__)
        try:
            del d['evmapping']
        except KeyError:
            print "'evmapping' local function not in self.__dict__"
        return d

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.makeCallFn()


def makeEvMapping(mappingDict, varnames, parnames):
    evMapDict = {}
    namemap = {}
    for varname in varnames:
        namemap[varname] = "xdict['"+varname+"']"
    for parname in parnames:
        namemap[parname] = "pdict['"+parname+"']"
    for k, v in mappingDict.iteritems():
        v_dummyQ = Symbolic.QuantSpec('dummy', v)
        v_dummyQ.mapNames(namemap)
        evMapDict["xdict['%s']"%k] = v_dummyQ()
    return EvMapping(evMapDict)


def makeGenInfoEntry(generator, allgen_names, swmap_list=[]):
    """Create an entry for the genInfo attribute of a Model."""

    assert isinstance(allgen_names, list), \
                             "'allgen_names' argument must be a list"
    special_reasons = ['time'] + generator.variables.keys()
    assert generator.name not in special_reasons + ['terminate'], \
         "Cannot use variable names or internal names 'time' and 'terminate' as generator names"
    try:
        allEndReasonNames = map(lambda (n,o):n,
                                generator.eventstruct.getTermEvents()) \
                            + special_reasons
    except AttributeError:
        # no events associated with the generator
        allEndReasonNames = special_reasons
    assert generator.name in allgen_names, ('Generator`s name not in list of all '
                                         'available Generator names!')
    assert alltrue([name not in allEndReasonNames for name in allgen_names]), \
           'Generator names overlapped with event or variable names'
    alltarg_names = allgen_names + ['terminate']
    # if no event map function specified, assume the identity fn
    seenReasons = []
    swmap_pairs = []
    if swmap_list != []:
        for mapentry in swmap_list:
            # check the entries of swmap_list and turn into a
            # (reason, infopair) pair, adding a default event map function
            # to some entries
            reason = mapentry[0]
            mapping_info = mapentry[1]
            if len(mapentry) > 2:
                raise ValueError("mapping entry must be (reason, infopair) tuple")
            if isinstance(mapping_info, tuple):
                genTargetName = mapping_info[0]
                numargs = len(mapping_info)
            elif isinstance(mapping_info, str):
                genTargetName = mapentry[1]
                numargs = 1
            else:
                raise TypeError("Invalid event mapping entry")
            if numargs == 2:
                epmap = mapping_info[1]
                assert isinstance(epmap, EvMapping), "Must supply EvMapping class"
                swmap_pairs.append((reason, mapping_info))
            elif numargs == 1:
                # use default identity mapping fn for event
                # and make this entry into a three-tuple
                swmap_pairs.append((reason, (genTargetName, EvMapping())))
            else:
                raise ValueError("Expected 2 or 3 arguments to Generator "
                                 "switch map entry")
            assert reason not in seenReasons, ('reason cannot appear more than'
                                               ' once in map domain')
            seenReasons.append(reason)
            assert reason in allEndReasonNames, ("name '"+reason+"' in map "
                                                "domain is missing")
            assert genTargetName in alltarg_names, ("name '"+genTargetName+"' in "
                                                "map range is missing")
    else:
        # There had better be only a single Generator in allgen_names,
        # otherwise we need a map
        assert len(allgen_names) == 1, ("There must be an event mapping "
                  "specified when there is more than one Generator in the Model")
    unseen_sr = remain(allEndReasonNames, seenReasons)
    if unseen_sr != []:
        # then there are 'end reasons' that do not have switch rules,
        # so give them defaults (terminate)
        for r in unseen_sr:
            swmap_pairs.append((r, ('terminate', EvMapping())))
    if len(swmap_pairs) != len(allEndReasonNames):
        info(dict(swmap_pairs))
        print "(%i in total),versus:"%len(swmap_pairs)
        print allEndReasonNames, "(%i in total)"%len(allEndReasonNames)
        raise ValueError('Incorrect number of map pairs given in argument')
    return {generator.name: (generator, dict(swmap_pairs))}


def processReused(sourcenames, auxvarnames, flatspec, registry,
                  FScompatibleNames, FScompatibleNamesInv):
    """Find and process reused terms in abstract specification. To avoid
    RHS specs depending on auxiliary variables, temp variables will be declared
    in FuncSpec.py and used in both the RHS and auxiliary variables in the
    target language specification."""
    reuseTerms={}
    subsExpr={}
    num_reused=0
    # auxvarnames are those that sourcename definitions cannot use
    # build auxiliary token map to get rid of auxvar - auxvar inter-
    # dependencies
    u_subsMap = {}
    for auxtok in auxvarnames:
        tokobj = registry[FScompatibleNamesInv(auxtok)].obj
        addtokbraces = tokobj.spec.isCompound()
#        u_new_reusedname = "__"+auxtok+str(num_reused)
        FScompat_spec = "".join(FScompatibleNames(tokobj.spec[:]))
        u_subsMap[auxtok] = "("*addtokbraces + \
                       FScompat_spec + ")"*addtokbraces
#        u_subsMap[auxtok] = "".join(FScompatibleNames( \
#                                     tokobj.spec[:]))
    # some of these u_subsMap targets may contain auxiliary variables
    # themselves, so we must purge them now in repeated passes to u_subsMap.
    # put in a trap for infinite loop of inter-dependencies!
    loopCount = 0
    loopMaxDepth = 15
    purgeDone = {}.fromkeys(auxvarnames, False)
    while not all(purgeDone.values()) and loopCount < loopMaxDepth:
        loopCount += 1
#        print "Loop count: ", loopCount
        tempMap = {}
        for auxtok, sx in u_subsMap.iteritems():
#            print "** ", auxtok
            if purgeDone[auxtok]:
#                print "  Continue 1"
                continue
            dummyQ = Symbolic.QuantSpec('dummy', sx)
            if not any([auxname in dummyQ \
                            for auxname in auxvarnames]):
                # no auxvar names appear in the subs expr, so this is cleared
                purgeDone[auxtok] = True
#                print "  Continue 2"
                continue
            dummyQ.mapNames(u_subsMap)
            tempMap[auxtok] = dummyQ()
        # update name map with any new substitutions
#        if tempMap != {}:
#            info(tempMap)
        u_subsMap.update(tempMap)
    if not purgeDone and len(auxvarnames)>0:
        # then must have maxed out
        print "Declared auxilary variables:", auxvarnames
        raise RuntimeError("You probably have an infinite loop of auxiliary "
                       "variable inter-dependencies: recursion depth of "
                       "more than %i encountered during model build"%loopCount)
    for v in sourcenames:
        if v not in flatspec['vars']:
            # v could be a parameter, a function name, or a constant (in a
            # recursive call), so ignore
            continue
        subsMap = {}
        dummyQ = Symbolic.QuantSpec('dummy', flatspec['vars'][v])
        for u in dummyQ.usedSymbols:
            if u in auxvarnames:
                new_reusedname = "__"+u
                if new_reusedname in reuseTerms.values():
                    # simple way to avoid name clashes
                    new_reusedname += '_'+str(num_reused)
                    num_reused += 1
                spec_text = flatspec['vars'][u]
                testQ = Symbolic.QuantSpec('dummy', spec_text)
                testQ.mapNames(mathNameMap)
                # add test for unary minus otherwise no braces around
                # testQ will lead to both signs disappearing on reuseTerm
                # substitution, leaving two symbols adjoined without any
                # operator!
                addbraces = testQ.isCompound() or testQ()[0] == '-'
                # no subs expression for auxvar that points to a constant
                noSubsExpr = not addbraces and \
                      (FScompatibleNamesInv(spec_text) in registry \
                                          or isNumericToken(spec_text))
                # make substitutions for any aux vars appearing in
                # spec_text (testQ)
                testQ.mapNames(u_subsMap)
                # update addbraces after mapping
                addbraces = testQ.isCompound() or testQ()[0] == '-'
                testQ.simplify()
                spec_text_new = "("*addbraces + testQ() + ")"*addbraces
#                spec_text_new = testQ()
                if not noSubsExpr:
                    if u in subsExpr:
                        # putting braces around auxtok in u_subsMap means
                        # that some of the expressions won't have the same
                        # bracketing as spec_text_new, so don't bother with
                        # this check
                        pass
#                        if subsExpr[u] != spec_text_new:
#                            print subsExpr[u]
#                            print spec_text_new
#                            raise RuntimeError("Different subs expr for %s in subsExpr"%u)
                    else:
                        subsExpr[u] = spec_text_new
                    if testQ()[0] == '-':
                        reuse_term = spec_text_new
                    else:
                        reuse_term = testQ()
                    if reuse_term not in reuseTerms:
                        reuseTerms[reuse_term] = new_reusedname
                if u in subsMap:
                    raise RuntimeError("%s already in subsMap!"%u)
                else:
                    subsMap[u] = spec_text_new
        # use QuantSpec's inbuilt tokenized version of exp_var definition
        # to make substitutions using the name mapping subsMap
        dummyQ.mapNames(subsMap)
        dummyQ.simplify()
        # uses addvbraces is use addbraces above, otherwise get clash
##        addvbraces = dummyQ.isCompound()
##        subsExpr[v] = "("*addvbraces + dummyQ() + ")"*addvbraces
        dummyQ.mapNames(mathNameMap)
        subsExpr[v] = dummyQ()
    return reuseTerms, subsExpr
