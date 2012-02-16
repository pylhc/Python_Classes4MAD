#! /sw/bin/python2.3
#Author:    Drew LaMar
#Date:      12 September 2006
#$Revision: 0.3.1 $
""" ContClass stores continuation curves for a specified model.

    Drew LaMar, March 2006
"""

# ----------------------------------------------------------------------------
# Version History:
#   Version 0.2.2, April 2006
#       Processes initial direction (initdirec)
#       Got rid of CurvesDict (redundant)
#       Added computeEigen to compute eigen information for all curves
#       Added exportGeomview method
#
#   Version 0.3.0, June 2006
#       Removed toggleLabels method
#       Added loadAutoMod, makeAutoLib, makeAutoLibSource, compileAutoLib,
#           forcedAutoLibRefresh methods.
#
#   Version 0.3.1, September 2006
#       BUG FIX: Added extra compile argument '-std=c99' to fix NAN issue on linux.
#           (also added '#include <math.h>' in auto_c.h)
#
# ----------------------------------------------------------------------------

from Continuation import EquilibriumCurve, FoldCurve, HopfCurveOne, HopfCurveTwo, \
                         FixedPointCurve, LimitCycleCurve
from misc import *
from Plotting import pargs, initializeDisplay

from PyDSTool.Model import Model, findInitialGenerator
from PyDSTool.Generator import Generator
from PyDSTool.ModelConstructor import embed
from PyDSTool import Point, Pointset
from PyDSTool.common import pickle, Utility, args
import PyDSTool.Redirector as redirc
try:
    from PyDSTool.matplotlib_import import *
except ImportError:
    from PyDSTool.matplotlib_unavailable import *
    print "Warning: matplotlib failed to import properly and so is not"
    print "  providing a graphing interface"

from scipy import Inf, NaN, isfinite, r_, c_, sign, mod, mat, \
    subtract, divide, matrixmultiply, transpose, eye, real, imag, \
    all
from scipy import array as scipy_array
from numarray import array, Float, Complex, Int, Float64, Complex64, Int32, \
    zeros, divide, subtract
from numarray import transpose as num_transpose
from numarray.numarraycore import NumArray
    
from PyDSTool.parseUtils import addArgToCalls, wrapArgInCall
from distutils.core import setup, Extension
from distutils.sysconfig import get_python_inc

import os, platform, shutil, sys
from PyDSTool import __path__ as _pydstool_path

#####
_pydstool_path = _pydstool_path[0]

_classes = ['ContClass']

_constants = ['curve_list', 'curve_args_list', 'auto_list']

__all__ = _classes + _constants
#####

rout = redirc.Redirector(redirc.STDOUT)
rerr = redirc.Redirector(redirc.STDERR)

curve_list = {'EP-C': EquilibriumCurve, 'LP-C': FoldCurve,
              'H-C1': HopfCurveOne, 'H-C2': HopfCurveTwo, 
              'FP-C': FixedPointCurve, 'LC-C': LimitCycleCurve}
              
curve_args_list = ['verbosity']

auto_list = ['LC-C']

class ContClass(Utility):
    """Stores continuation curves for a specified model."""
    def __init__(self, model):
        if isinstance(model, Generator):
            self.model = embed(model)
        else:
            self.model = model
        self.gensys = findInitialGenerator(self.model.genInfo,
                                0, self.model.icdict, self.model.pars,
                    self.model.intvars)[0]
        self._autoMod = None
        self.curves = {}
        self.plot = pargs()

    def __getitem__(self, name):
        if not self.curves.has_key(name):
            raise KeyError('No curve named ' + name)
        else:
            return self.curves[name]

    def __copy__(self):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)

    def __deepcopy__(self, memo=None, _nil=[]):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)

    def newCurve(self, initargs):
        """Create new curve with arguments specified in the dictionary initargs."""
        curvetype = initargs['type']
        curvetype = curvetype.upper()
        
        if curvetype not in curve_list:
            raise TypeError(curvetype + ' not an allowable curve type')

        # Check name
        if initargs['name'] in self.curves.keys():
            raise AttributeError('Ambiguous name field: ' + initargs['name'] + ' already exists')
        
        # Check parameters
        if self.model.pars == {}:
            raise ValueError('No parameters defined for this system!')

        # Process initial point
        if 'initpoint' not in initargs:
            # Default to initial conditions for model
            if self.model.icdict == {}:
                raise ValueError('No initial point defined for this system!')
            initargs['initpoint'] = self.model.icdict.copy()
            #for p in initargs['freepars']:
            #    initargs['initpoint'][p] = self.model.pars[p]
        else:
            if isinstance(initargs['initpoint'], dict):
                initargs['initpoint'] = initargs['initpoint'].copy()
                #for p in initargs['freepars']:
                #    if p not in initargs['initpoint'].keys():
                #        initargs['initpoint'][p] = self.model.pars[p]
            elif isinstance(initargs['initpoint'], str):
                curvename, pointname = initargs['initpoint'].split(':')
                pointtype = pointname.strip('0123456789')
                if not self.curves.has_key(curvename):
                    raise KeyError('No curve of name ' + curvename + ' exists.')
                else:
                    point = self.curves[curvename].getSpecialPoint(pointtype, pointname)
                    if point is None:
                        raise KeyError('No point of name ' + pointname + ' exists.')
                    else:
                        initargs['initpoint'] = point
                    
            # Separate from if-else above since 'str' clause returns type Point
            if isinstance(initargs['initpoint'], Point):
                # Check to see if point contains a cycle.  If it does, assume
                #   we are starting at a cycle and save it in initcycle
                for v in initargs['initpoint'].labels.itervalues():
                    if v.has_key('cycle'):
                        initargs['initcycle'] = v   # Dictionary w/ cycle, name, and tangent information

                # Save initial point information
                initargs['initpoint'] = initargs['initpoint'].todict()
                #for p in initargs['freepars']:
                #    if p not in initargs['initpoint'].keys():
                #        initargs['initpoint'][p] = self.model.pars[p]
                        
        # Process cycle
        if 'initcycle' in initargs:
            if isinstance(initargs['initcycle'], NumArray):
                c0 = {}
                c0['data'] = args(V = {'udotps': None, 'rldot': None})
                c0['cycle'] = Pointset({'coordnames': self.gensys.funcspec.vars,
                                        'coordarray': num_transpose(initargs['initcycle'][:,1:]).copy(),
                                        'indepvarname': 't',
                                        'indepvararray': num_transpose(initargs['initcycle'][:,0]).copy()
                                       })
                initargs['initcycle'] = c0
            elif isinstance(initargs['initcycle'], Pointset):
                c0 = {}
                c0['data'] = args(V = {'udotps': None, 'rldot': None})
                c0['cycle'] = initargs['initcycle']
                initargs['initcycle'] = c0

        # Load auto module if required
        automod = None
        if curvetype in auto_list:
            if self._autoMod is None:
                self.loadAutoMod()
            automod = self._autoMod
            
        self.curves[initargs['name']] = curve_list[curvetype](self.model, self.gensys, automod, self.plot, initargs)

    def exportGeomview(self, coords=None, filename="geom.dat"):
        if  coords is not None and len(coords) == 3:
            GeomviewOutput = "(progn (geometry " + self.model.name + " { LIST {: axes_" + self.model.name + "}"
            for cname, curve in self.curves.iteritems():
                GeomviewOutput += " {: " + cname + "}"
            GeomviewOutput += "}))\n\n"
            
            # Get axes limits
            alim = [[Inf,-Inf],[Inf,-Inf],[Inf,-Inf]]
            for cname, curve in self.curves.iteritems():
                for n in range(len(coords)):
                    alim[n][0] = min(alim[n][0], min(curve.sol[coords[n]].toarray()))
                    alim[n][1] = max(alim[n][1], max(curve.sol[coords[n]].toarray()))
            
            GeomviewOutput += "(progn (hdefine geometry axes_" + self.model.name + " { appearance { linewidth 2 } SKEL 4 3 " + \
                "0 0 0 1 0 0 0 1 0 0 0 1 " + \
                "2 0 1 1 0 0 1 2 0 2 0 1 0 1 2 0 3 0 0 1 1})\n\n"
                
            for cname, curve in self.curves.iteritems():
                GeomviewOutput += "(hdefine geometry " + cname + " { LIST {: curve_" + cname + "} {: specpts_" + cname + "}})\n\n"
                
                GeomviewOutput += "(hdefine geometry curve_" + cname + " { appearance { linewidth 2 } SKEL " + \
                    repr(len(curve.sol)) + " " + repr(len(curve.sol)-1)
                for n in range(len(curve.sol)):
                    GeomviewOutput += " " + repr((curve.sol[n][coords[0]]-alim[0][0])/(alim[0][1]-alim[0][0])) + \
                                           " " + repr((curve.sol[n][coords[1]]-alim[1][0])/(alim[1][1]-alim[1][0])) + \
                                           " " + repr((curve.sol[n][coords[2]]-alim[2][0])/(alim[2][1]-alim[2][0]))
                for n in range(len(curve.sol)-1):
                    GeomviewOutput += " 2 " + repr(n) + " " + repr(n+1) + " 0 0 0 1"
                
                GeomviewOutput += "})\n\n"
                
            GeomviewOutput += ")\n"
            
            f = open(filename, "w")
            f.write(GeomviewOutput)
            f.close()
        else:
            raise Warning("Coordinates not specified or not of correct dimension.")
        
    def display(self, coords=None, curves=None, figure=None, axes=None, stability=False, domain=False, **plot_args):
        """Plot all curves in coordinates specified by coords.
        
           Inputs:
                    
               coords -- tuple of coordinates (None defaults to the first free
                   parameter and the first state variable)
        """
        if coords is not None and len(coords) == 3:
            self.exportGeomview(coords)
            return
            
        if curves is None:
            curves = self.curves.keys()

        plot_curves = []
        for curve in curves:
            if self.curves.has_key(curve):
                plot_curves.append(curve)
            else:
                print "Warning: Curve " + curve + " does not exist."
                
        if len(plot_curves) > 0:
            initializeDisplay(self.plot, figure=figure, axes=axes)
            
        for curve in plot_curves:
            self.curves[curve].display(coords, figure=figure, axes=axes, stability=stability, domain=domain, init_display=False, **plot_args)
        
    def computeEigen(self):
        for curve in self.curves.itervalues():
            curve.computeEigen()

    def info(self):
        print self.__repr__()
        print "  Variables : %s"%', '.join(self.model.allvars)
        print "  Parameters: %s\n"%', '.join(self.model.pars.keys())
        print "Containing curves: "
        for c in self.curves:
            print "  " + c + " (type " + self.curves[c].curvetype + ")"
            
    def update(self, args):
        """Update parameters for all curves."""
        for c in args.keys():
            if c not in curve_args_list:
                args.pop(c)
                
        for v in self.curves.values():
            v.update(args)
            
    def loadAutoMod(self):
        thisplatform = platform.system()
        if thisplatform == 'Windows':
            self._dllext = ".pyd"
        elif thisplatform in ['Linux', 'IRIX', 'Solaris', 'SunOS', 'Darwin']:
            self._dllext = '.so'
        else:
            print "Shared library extension not tested on this platform."
            print "If this process fails please report the errors to the"
            print "developers."
            self._dllext = '.so'
        
        self._compilation_tempdir = os.path.join(os.getcwd(),
                                                      "auto_temp")
        if not os.path.isdir(self._compilation_tempdir):
            try:
                assert not os.path.isfile(self._compilation_tempdir), \
                     "A file already exists with the same name"
                os.mkdir(self._compilation_tempdir)
            except:
                print "Could not create compilation temp directory " + \
                      self._compilation_tempdir
                raise
        self._compilation_sourcedir = os.path.join(_pydstool_path,"PyCont/auto/module")
        self._vf_file = self.gensys.name+"_vf.c"
        self._vf_filename_ext = "_"+self._vf_file[:-2]
        if not (os.path.isfile(os.path.join(os.getcwd(),
                                "auto"+self._vf_filename_ext+".py")) and \
                os.path.isfile(os.path.join(os.getcwd(),
                                "_auto"+self._vf_filename_ext+self._dllext))):
            if True:
                self.funcspec = self.gensys.funcspec.recreate('c')
                self.makeAutoLibSource()
                self.compileAutoLib()
            else:
                print "Build the library using the makeLib method, or in "
                print "stages using the makeLibSource and compileLib methods."

        try:
            self._autoMod = __import__("auto"+self._vf_filename_ext, globals())
        except:
            print "Error loading auto module."
            raise
        
    def forceAutoLibRefresh(self):
        """forceAutoLibRefresh should be called after event contents are changed,
        or alterations are made to the right-hand side of the ODEs.

        Currently this function does NOT work!"""

        # (try to) free auto module from namespace
        delfiles = True
        try:
            del(sys.modules["_auto"+self._vf_filename_ext])
            del(sys.modules["auto"+self._vf_filename_ext])
        except NameError:
            # modules weren't loaded, so nothing to do
            delfiles = False
        if delfiles:
            gc.collect()
            # still not able to delete these files!!!!! Argh!
        print "Cannot rebuild library without restarting session. Sorry."
        print "Try asking the Python developers to make a working module"
        print "unimport function!"

    def makeAutoLib(self, libsources=[], libdirs=[], include=[]):
        """makeAutoLib calls makeAutoLibSource and then the compileAutoLib method.
        To postpone compilation of the source to a DLL, call makeAutoLibSource()
        separately."""
        self.makeAutoLibSource(include)
        self.compileAutoLib(libsources, libdirs)

    def makeAutoLibSource(self, include=[]):
        """makeAutoLibSource generates the C source for the vector field specification.
        It should be called only once per vector field."""

        # Make vector field (and event) file for compilation
        assert isinstance(include, list), "includes must be in form of a list"
        # codes for library types (default is USERLIB, since compiler will look in standard library d
        STDLIB = 0
        USERLIB = 1
        libinclude = dict([('math.h', STDLIB), ('stdio.h', STDLIB), ('stdlib.h', STDLIB),
                      ('string.h', STDLIB), ('autovfield.h', USERLIB),
                      ('auto_c.h', USERLIB)])
        include_str = '#include "auto_f2c.h"\n' # This must come first
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
        
        # f2c auto conventions (dirty trick!)
        #define_str = "\n#define double    doublereal\n#define int   integer\n\n"
        define_str = ""
        
        allfilestr = "/*  Vector field and other functions for Auto continuer.\n " \
            + "  This code was automatically generated by PyDSTool, but may be modified " \
            + "by hand. */\n\n" + include_str + define_str + """
double *gICs;
double **gBds;
double globalt0;

"""
        pardefines = ""
##        parundefines = ""
        vardefines = ""
##        varundefines = ""
        inpdefines = ""
##        inpundefines = ""
        # sorted version of var, par, and input names
        vnames = self.gensys._var_ixmap
        pnames = self.funcspec.pars
        inames = self.funcspec.inputs
        pnames.sort()
        inames.sort()
        for i in xrange(self.gensys.numpars):
            p = pnames[i]
            # add to defines (WATCH OUT FOR PERIOD _T!!!)
            if (i < 10):
                pardefines += self.funcspec._defstr+" "+p+"\tp_["+str(i)+"]\n"
            elif (i >= 10):
                pardefines += self.funcspec._defstr+" "+p+"\tp_["+str(i+40)+"]\n"
        # add period _T
        pardefines += self.funcspec._defstr+" _T\tp_[10]\n"
        for i in xrange(self.gensys.dimension):
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
        # add signature for auxiliary functions
        if self.funcspec.auxfns:
            allfilestr += "\n"
            for finfo in self.funcspec.auxfns.values():
                allfilestr += finfo[1] + ";\n"
        allfilestr += "\nvoid auxvars(unsigned, unsigned, double, double*, double*, " \
              + "double*, unsigned, double*, unsigned, double*);\n" \
              + """void jacobian(unsigned, unsigned, double, double*, double*, double**, unsigned, double*, unsigned, double*);
void jacobianParam(unsigned, unsigned, double, double*, double*, double**, unsigned, double*, unsigned, double*);
"""
        if self.funcspec.auxvars == []:
            allfilestr += "int N_AUXVARS = 0;\n\n\n"
        else:
            allfilestr += "int N_AUXVARS = " + str(len(self.funcspec.auxvars)) \
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
        allfilestr += self.funcspec.auxspec[0]
        # if jacobians or mass matrix not present, fill in dummy
        if not self.gensys.haveJacobian():
            allfilestr += """
void jacobian(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}
"""
        if not self.gensys.haveJacobian_pars():
            allfilestr += """
void jacobianParam(unsigned n_, unsigned np_, double t, double *Y_, double *p_, double **f_, unsigned wkn_, double *wk_, unsigned xvn_, double *xv_) {
}
""" #+ "\n/* Variable and parameter substitutions undefined:*/\n" + parundefines + varundefines + "\n" 
        # write out C file
        vffile = os.path.join(self._compilation_tempdir, self._vf_file)
        try:
            file = open(vffile, 'w')
            #allfilestr = allfilestr.replace("double ","doublereal ")
            #allfilestr = allfilestr.replace("double*","doublereal *")
            #allfilestr = allfilestr.replace("double**","doublereal **")
            file.write(allfilestr)
            file.close()
        except IOError, e:
            print "Error opening file "+self._vf_file+" for writing"
            raise IOError, e

    def compileAutoLib(self, libsources=[], libdirs=[]):
        """compileAutoLib generates a python extension DLL with continuer and vector
        field compiled and linked.

        libsources list allows additional library sources to be linked.
        libdirs list allows additional directories to be searched for
          precompiled libraries."""

        if os.path.isfile(os.path.join(os.getcwd(),
                                "_auto"+self._vf_filename_ext+self._dllext)):
            # then DLL file already exists and we can't overwrite it at this
            # time
            proceed = False
            print "\n"
            print "-----------------------------------------------------------"
            print "Present limitation of Python: Cannot rebuild library"
            print "without exiting Python and deleting the shared library"
            print "   " + str(os.path.join(os.getcwd(),
                                "_auto"+self._vf_filename_ext+self._dllext))
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
        if self._autoMod is not None:
            self.forceAutoLibRefresh()
        vffile = os.path.join(self._compilation_tempdir, self._vf_file)
        try:
            ifacefile_orig = open(os.path.join(self._compilation_sourcedir,
                                               "automod.i"), 'r')
            ifacefile_copy = open(os.path.join(self._compilation_tempdir,
                                       "auto_"+self._vf_file[:-2]+".i"), 'w')
            firstline = ifacefile_orig.readline()
            ifacefile_copy.write('%module auto_'+self._vf_file[:-2]+'\n')
            iffilestr = ifacefile_orig.read()
            ifacefile_copy.write(iffilestr)
            ifacefile_orig.close()
            ifacefile_copy.close()
        except IOError:
            print "automod.i copying error in auto compilation directory"
            raise
            
        swigfile = os.path.join(self._compilation_tempdir,
                                "auto"+self._vf_filename_ext+".i")
        automodfile = os.path.join(self._compilation_sourcedir, "automod.c")
        interfacefile = os.path.join(self._compilation_sourcedir, "interface.c")
        
        # source files
        if not (all([os.path.isfile(os.path.join(self._compilation_tempdir,
                           sf)) for sf in ['auto'+self._vf_filename_ext+'_wrap.o',
                                           'auto'+self._vf_filename_ext+'.py',
                                           '_auto'+self._vf_filename_ext+'.def']])):
            modfilelist = [swigfile]
        else:
            modfilelist = []
        # FOR DIST (ADD)
        modfilelist.extend([os.path.join(self._compilation_sourcedir, "../src/"+x) \
            for x in ['auto.c','autlib1.c','autlib2.c','autlib3.c','autlib4.c','autlib5.c', \
            'eispack.c','conpar.c','setubv.c','reduce.c','dmatrix.c','fcon.c','libf2c/cabs.c','libf2c/d_lg10.c', \
            'libf2c/i_nint.c','libf2c/pow_di.c','libf2c/r_lg10.c','libf2c/z_exp.c','libf2c/d_imag.c', \
            'libf2c/d_sign.c','libf2c/i_dnnt.c','libf2c/pow_dd.c','libf2c/pow_ii.c','libf2c/z_abs.c', \
            'libf2c/z_log.c']])
        
        modfilelist.extend([automodfile, interfacefile, vffile])
        # FOR DIST (SUBTRACT)
        #modfilelist.extend(libsources)
        
        # script args
        script_args = ['-q', 'build', '--build-lib='+os.getcwd(), # '-t/',
                 '-t'+self._compilation_tempdir,
                 '--build-base='+self._compilation_sourcedir]
        if self.gensys._compiler != '':
            script_args.append('-c'+str(self.gensys._compiler))
            
        # include directories for libraries
        narraydir = os.path.join(get_python_inc(plat_specific=1), "numarray")
        incdirs = [narraydir, os.getcwd(), os.path.join(self._compilation_sourcedir,"include"),
            self._compilation_tempdir, os.path.join(_pydstool_path,"PyCont/auto/src/include")]
        incdirs.extend(libdirs)
        
        # libraries
        # FOR DIST (SUBTRACT)
        #libdirs.append(os.path.join(_pydstool_path, "PyCont/auto/lib"))
        #libsources.append('auto2000')
        
        # Use distutils to perform the compilation of the selected files
        rout.start()    # redirect stdout
        try:
            distobject = setup(name = "Auto 2000 continuer",
                  author = "PyDSTool (automatically generated)",
                  script_args = script_args,
                  ext_modules = [Extension("_auto"+self._vf_filename_ext,
                                 sources=modfilelist,
                                 include_dirs=incdirs,
                                 extra_compile_args=['-w', '-D__PYTHON__', '-std=c99'],
                                 library_dirs=libdirs,
                                 libraries=libsources)])
        except:
            print "\nError occurred in generating Auto system..."
            print sys.exc_info()[0], sys.exc_info()[1]
            raise RuntimeError
        rout.stop()    # restore stdout
        try:
            # move library files into the user's CWD
            if swigfile in modfilelist or not \
               os.path.isfile(os.path.join(self._compilation_tempdir,
                                "auto"+self._vf_filename_ext+".py")):
                shutil.move(os.path.join(self._compilation_tempdir,
                                 "auto"+self._vf_filename_ext+".py"),
                            os.path.join(os.getcwd(),
                                 "auto"+self._vf_filename_ext+".py"))
        except:
            print "\nError occurred in generating Auto system"
            print "(while moving library extension modules to CWD)"
            print sys.exc_info()[0], sys.exc_info()[1]
            raise RuntimeError

    def __repr__(self):
        return 'ContClass of model %s'%self.model.name

    __str__ = __repr__
