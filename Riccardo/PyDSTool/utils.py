#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      19 Sep 2006
#$Revision: 1.0.8 $
"""
    User utilities.

    Robert Clewley, October 2005.
"""

# ----------------------------------------------------------------------------
# VERSION HISTORY
#
# Version 1.0.8, Sep 2006
#   Added function find to find elements of a 1D array.
#   Added stdout redirection for calls to minpack functions.
#
# Version 1.0.7, Feb 2006
#   Fixed bug when including independent variable in exportPointset call.
#   Added compareList function to check that equality of elements (like sets).
#
# Version 1.0.6, Jan 2006
#   Minor bug fixes.
#   Updated exportPointset to accept empty list of coord names to default to
#     all coord names, and to make use of scipy.io.write_array function.
#
# Version 1.0.5, Nov 2005
#   Moved convertPowers to ModelSpec.py. (now in parseUtils.py -- Mar 06)
#   Added findClosestArray.
#   Added use of cPickle for non-Win32 platforms, to improve efficiency of
#     saveObjects, loadObjects
#   Added makeMfileFunction for partial Matlab support.
#
# Version 1.0.4, Sep 2005
#   Added info pretty printer utility.
#   convertPowers now accepts x**p as well as x^p using additional argument.
#
# Version 1.0.3, Aug 2005
#   Added convertPowers() function
#   Changed name of importDataFile to importPointset and exportDataFiles
#     to exportPointset, and added functionality to importer.
#
# Version 1.0.2, July 2005
#   Added functionality to importDataFile function to support
#     internally-provided time array.
#
# Version 1.0.1, June 2005
#   Added saveObjects and loadObjects functions.
#   Using a fixed version of Python pickle to get around bug in loading
#      of Inf and Nan.
#
# ----------------------------------------------------------------------------

from errors import *
from common import *
import Redirector as redirc
from parseUtils import joinStrs

from scipy import Inf, NaN, isfinite, less, greater, sometrue, alltrue, \
     searchsorted, take, argsort
from scipy import integer as Numeric_Int
from numarray import array, Float, Int, Complex, Float64, Int32, \
        Complex64, Complex32, Float32, swapaxes, asarray, zeros, transpose
from scipy.optimize import minpack
from scipy.io import write_array
import time, sys, os
import copy

# --------------------------------------------------------------------

# EXPORTS

_classes = []

_functions = ['intersect', 'remain', 'orderEventData', 'makeDataDict',
              'importPointset', 'exportPointset', 'makeImplicitFunc',
              'saveObjects', 'loadObjects', 'info', 'compareList',
              'findClosestArray', 'makeMfileFunction', 'find']

_mappings = ['_implicitSolveMethods', '_1DimplicitSolveMethods']

__all__ = _classes + _functions + _mappings


rout = redirc.Redirector(redirc.STDOUT)
rerr = redirc.Redirector(redirc.STDERR)

## ------------------------------------------------------------------

## Utility functions

def makeMfileFunction(name, argname, defs):
    """defs is a dictionary of left-hand side -> right-hand side definitions"""
    # writeout file <name>.m
    mfile = open(name+".m", 'w')
    mfile.write("function %s = %s(%s)\n"%(name,name,argname))
    for k, v in defs.iteritems():
        if k != name:
            mfile.write("%s = %s;\n"%(k,v))
    # now the final definition of tau_recip or inf
    mfile.write("%s = %s;\n"%(name,defs[name]))
    mfile.write("return\n")
    mfile.close()


def info(x, specName="Contents", offset=1, recurseDepth=1,
                 recurseDepthLimit=2):
    """Pretty printer for showing argument lists and dictionary
    specifications."""

    if hasattr(x, 'iteritems'):
        print specName + ":"
        x_keys = sortedDictKeys(x)
        for k in x_keys:
            v = x[k]
            kstr = object2str(k)
            basestr = " "*(offset-1) + kstr
            if isinstance(v, dict):
                info(v, basestr, offset+4, recurseDepth+1,
                             recurseDepthLimit)
            else:
                vStrList = object2str(v).split(', ')
                outStrList = [basestr+": "]
                for i in range(len(vStrList)):
                    if len(vStrList[i] + outStrList[-1]) < 78:
                        outStrList[-1] += ", "*(i>0) + vStrList[i]
                    else:
                        if i>0:
                            if i != len(vStrList):
                                # add trailing comma to previous line
                                outStrList[-1] += ","
                            # start on new line
                            outStrList.append(" "*(len(kstr)+3) + vStrList[i])
                        else:
                            # too long for line and string has no commas
                            # could work harder here, but for now, just include
                            # the long line
                            outStrList[-1] += vStrList[i]
                if recurseDepth==1 and len(outStrList)>1:
                    # print an extra space between topmost level entries
                    # provided those entries occupy more than one line.
                    print "\n"
                for s in outStrList:
                    print s
    elif hasattr(x, '__dict__') and recurseDepth <= recurseDepthLimit:
        info(x.__dict__, specName, offset, recurseDepth,
                     recurseDepthLimit)
    else:
        print x


_implicitSolveMethods = ['newton', 'bisect', 'steffe', 'fsolve']
_1DimplicitSolveMethods = ['newton', 'bisect', 'steffe']


def makeImplicitFunc(f, x0, fprime=None, extrafargs=(), xtolval=1e-8,
                        maxnumiter=100, solmethod='newton', standalone=True):
    """Builds an implicit function representation of an N-dimensional curve
    specified by (N-1) equations. Thus argument f is a function of 1 variable.
    In the case of the 'fsolve' method, f may have dimension up to N-1.

    Available solution methods are: newton, bisect, steffensen, fsolve.
    All methods utilize SciPy's Minpack wrappers to Fortran codes.

    Steffenson uses Aitken's Delta-squared convergence acceleration.
    fsolve uses Minpack's hybrd and hybrj algorithms.

    Standalone option (True by default) returns regular function. If False,
    an additional argument is added, so as to be compatible as a method
    definition."""
    
    if solmethod == 'bisect':
        assert type(x0) in [Array, NArray, list, tuple], \
               "Invalid type '"+str(type(x0))+"' for x0 = "+str(x0)
        assert len(x0) == 2
    elif solmethod == 'fsolve':
        assert type(x0) in [Array, NArray, list, tuple, int, float], \
               "Invalid type '"+str(type(x0))+"' for x0 = "+str(x0)
    else:
        assert type(x0) in [int, float], \
               "Invalid type '"+str(type(x0))+"' for x0 = "+str(x0)

    # define the functions that could be used
    # scipy signatures use y instead of t, but this naming is consistent
    # with that in the Generator module
    try:
        if standalone:
            def newton_fn(t):
                rout.start()
##                rerr.start()
                res = minpack.newton(f, x0, args=(t,)+extrafargs, tol=xtolval,
                                     maxiter=maxnumiter, fprime=fprime)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                # return {'result': res, 'warnings': warns}
                return res

            def bisect_fn(t):
                rout.start()
##                rerr.start()
                res = minpack.bisection(f, x0[0], x0[1], args=(t,)+extrafargs,
                                      xtol=xtolval, maxiter=maxnumiter)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                return res

            def steffe_fn(t):
                rout.start()
##                rerr.start()
                res = minpack.fixed_point(f, x0, args=(t,)+extrafargs,
                                           xtol=xtolval, maxiter=maxnumiter)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                return res

            def fsolve_fn(t):
                rout.start()
##                rerr.start()
                res = minpack.fsolve(f, x0, args=(t,)+extrafargs,
                                      xtol=xtolval, maxfev=maxnumiter,
                                      fprime=fprime)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                return res
        else:
            def newton_fn(s, t):
                rout.start()
##                rerr.start()
                res = minpack.newton(f, x0, args=(t,)+extrafargs, tol=xtolval,
                                     maxiter=maxnumiter, fprime=fprime)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                return res

            def bisect_fn(s, t):
                rout.start()
##                rerr.start()
                res = minpack.bisection(f, x0[0], x0[1], args=(t,)+extrafargs,
                                      xtol=xtolval, maxiter=maxnumiter)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                return res

            def steffe_fn(s, t):
                rout.start()
##                rerr.start()
                res = minpack.fixed_point(f, x0, args=(t,)+extrafargs,
                                           xtol=xtolval, maxiter=maxnumiter)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                return res

            def fsolve_fn(s, t):
                rout.start()
##                rerr.start()
                res = minpack.fsolve(f, x0, args=(t,)+extrafargs,
                                      xtol=xtolval, maxfev=maxnumiter,
                                      fprime=fprime)
                rout.stop()
##                warns = rout.stop()
##                rerr.stop()
                return res

    except TypeError, e:
        if solmethod == 'bisect':
            infostr = " (did you specify a pair for x0?)"
        else:
            infostr = ""
        raise TypeError("Could not create function" +infostr + ": "+str(e))

    if solmethod == 'newton':
        return newton_fn
    elif solmethod == 'bisect':
        if fprime is not None:
            print "Warning: fprime argument unused for bisection method"
        return bisect_fn
    elif solmethod == 'steffe':
        if fprime is not None:
            print "Warning: fprime argument unused for aitken method"
        return steffe_fn
    elif solmethod == 'fsolve':
        return fsolve_fn
    else:
        raise ValueError, "Unrecognized type of implicit function solver"


def findClosestArray(input_array, target_array, tol):
    """

    Find the set of elements in input_array that are closest to
    elements in target_array.  Record the indices of the elements in
    target_array that are within tolerance, tol, of their closest
    match. Also record the indices of the elements in target_array
    that are outside tolerance, tol, of their match.

    For example, given an array of observations with irregular
    observation times along with an array of times of interest, this
    routine can be used to find those observations that are closest to
    the times of interest that are within a given time tolerance.

    NOTE: input_array must be sorted! The array, target_array, does not have to be sorted.

    Inputs:
      input_array:  a sorted Float64 numarray
      target_array: a Float64 numarray
      tol:          a tolerance

    Returns:
      closest_indices:  the array of indices of elements in input_array that are closest to elements in target_array
      accept_indices:  the indices of elements in target_array that have a match in input_array within tolerance
      reject_indices:  the indices of elements in target_array that do not have a match in input_array within tolerance
    
    Author: Gerry Wiener, 2004
    Version 1.0
    """

    input_array_len = len(input_array)
    closest_indices = searchsorted(input_array, target_array) # determine the locations of target_array in input_array
#    acc_rej_indices = [-1] * len(target_array)
    curr_tol = [tol] * len(target_array)

    est_tol = 0.0
    for i in xrange(len(target_array)):
        best_off = 0          # used to adjust closest_indices[i] for best approximating element in input_array

        if closest_indices[i] >= input_array_len:
            # the value target_array[i] is >= all elements in input_array so check whether it is within tolerance of the last element
            closest_indices[i] = input_array_len - 1
            est_tol = target_array[i] - input_array[closest_indices[i]]
            if est_tol < curr_tol[i]:
                curr_tol[i] = est_tol
#                acc_rej_indices[i] = i
        elif target_array[i] == input_array[closest_indices[i]]:
            # target_array[i] is in input_array
            est_tol = 0.0
            curr_tol[i] = 0.0
#            acc_rej_indices[i] = i
        elif closest_indices[i] == 0:
            # target_array[i] is <= all elements in input_array
            est_tol = input_array[0] - target_array[i]
            if est_tol < curr_tol[i]:
                curr_tol[i] = est_tol
#                acc_rej_indices[i] = i
        else:
            # target_array[i] is between input_array[closest_indices[i]-1] and input_array[closest_indices[i]]
            # and closest_indices[i] must be > 0
            top_tol = input_array[closest_indices[i]] - target_array[i]
            bot_tol = target_array[i] - input_array[closest_indices[i]-1]
            if bot_tol <= top_tol:
                est_tol = bot_tol
                best_off = -1           # this is the only place where best_off != 0
            else:
                est_tol = top_tol

            if est_tol < curr_tol[i]:
                curr_tol[i] = est_tol
#                acc_rej_indices[i] = i

        if est_tol <= tol:
            closest_indices[i] += best_off

#    accept_indices = compress(greater(acc_rej_indices, -1),
#                                       acc_rej_indices)
#    reject_indices = compress(equal(acc_rej_indices, -1),
#                                       arange(len(acc_rej_indices)))
    return closest_indices #, accept_indices, reject_indices)


def find(x, v, next_largest=1, indices=None):
    """Returns the index into the 1D array x corresponding to the
    element of x that is either equal to v or the nearest to
    v. x is assumed to contain unique elements.

    if v is outside the range of values in x then the index of the
    smallest or largest element of x is returned.

    If next_largest == 1 then the nearest element taken is the next
    largest, otherwise if next_largest == 0 then the next smallest
    is taken.

    The optional argument indices speeds up multiple calls to this
    function if you pre-calculate indices=argsort(x).
    """
    if indices is None:
        indices=argsort(x)
    xs=take(x, indices)
    assert next_largest in [0,1], "next_largest must be 0 or 1"
    eqmask=(xs==v).tolist()
    try:
        ix = eqmask.index(1)
    except ValueError:
        if next_largest:
            mask=(xs<v).tolist()
        else:
            mask=(xs>v).tolist()
        try:
            ix=min([max([0,mask.index(1-next_largest)+next_largest-1]),len(mask)-1])
        except ValueError:
            ix = 0+next_largest-1
    return indices[ix]


def orderEventData(edict, evnames=None, nonames=False):
    """Time-order event data dictionary items.

    Returns time-ordered list of (eventname, time) tuples.
    If 'evnames' argument included, this restricts output to only the named
    events.
    The 'nonames' flag forces routine to return only the event times, with
    no associated event names.
    """
    
    if evnames is None:
        evnames = edict.keys()
    else:
        assert remain(evnames, edict.keys()) == [], "Invalid event names passed"
    # put times as first tuple entry of etuplelist
    if nonames:
        alltlist = []
        for (evname,tlist) in edict.items():
            if evname in evnames:
                alltlist.extend(tlist)
        alltlist.sort()
        return alltlist
    else:
        etuplelist = []
        for (evname,tlist) in edict.items():
            if evname in evnames:
                etuplelist.extend([(t,evname) for t in tlist])
        # sort by times
        etuplelist.sort()
        # swap back to get event names as first tuple entry
        return [(evname,t) for (t,evname) in etuplelist]


def saveObjects(objlist, filename, force=False):
    """Store PyDSTool objects to file."""
        
    # passing protocol = -1 to pickle means it uses highest available
    # protocol (e.g. binary format)
    if not force:
        if os.path.isfile(filename):
            raise ValueError("File '" + filename + "' already exists")
    pklfile = file(filename, 'wb')
    # Win32 only: in call to pickle.dump ...
    # DO NOT use binary option (or HIGHESTPROTOCOL) because
    # IEE754 special values are not UNpickled correctly in Win32
    # (you'll see no exception raised). This is a known bug
    # and fixedpickle.py is a work-around. (June 2005)
    if os.name == 'nt':
        opt = None
    else:
        opt = 0
    if not isinstance(objlist, list):
        objlist=[objlist]
    for obj in objlist:
        try:
            pickle.dump(obj, pklfile, opt)
        except:
            if hasattr(obj, 'name'):
                print "Failed to save '%s'"%obj.name
            else:
                print "Failed to save object '%s'"%str(obj)
            raise
    pklfile.close()



def loadObjects(filename, namelist=[]):
    """Retrieve PyDSTool objects from file. Returns list of objects.

    Optional namelist argument selects objects to return by name,
    provided that the objects have name fields (otherwise they are ignored)."""

    # Since names are not intended to be unique in PyDSTool, the while
    # loop always goes to the end of the file, and pulls out *all*
    # occurrences of the names.
    if not os.path.isfile(filename):
        raise ValueError("File '" + filename + "' not found")
    if not isinstance(namelist, list):
        if isinstance(namelist, str):
            namelist = [copy.copy(namelist)]
        else:
            raise TypeError("namelist must be list of strings or singleton string")
    if not uniqueList(namelist):
        raise ValueError("Names must only appear once in namelist argument")
    pklfile = file(filename, 'rb')
    if namelist == []:
        getall = True
    else:
        getall = False
    objlist = []
    notDone = True
    while notDone:
        try:
            if getall:
                objlist.append(pickle.load(pklfile))
            else:
                tempobj = pickle.load(pklfile)
                if hasattr(tempobj, 'name'):
                    if tempobj.name in namelist:
                        objlist.append(tempobj)
        except EOFError:
            notDone = False
        except:
            print "Error in un-pickling:"
            print "Was the object created with an old version of PyDSTool?"
            pklfile.close()
            raise
    pklfile.close()
    if objlist == []:
        if getall:
            print "No objects found in file"
        else:
            print "No named objects found in file"
    return objlist


# convert pointset (e.g. sampled from a trajectory) to a
# set of ascii whitespace- (or user-defined character-) separated data
# files. option to list each variable's data in rows ('across') or in
# columns ('down'). existing files of the same names will be overwritten,
# unless the 'append' boolean option is set.
# NB if argument 'ext' present without a leading dot, one will be added
#
# infodict should consist of: keys = filenames, values = tuples of pointset
#   variable names to export
def exportPointset(thepointset, infodict, separator='   ', linesep='\n',
                   precision=12, suppress_small=0, varvaldir='col',
                   ext='', append=False):

    assert varvaldir in ['col', 'row'], \
           "invalid variable value write direction"
    # in order to avoid import cycles, cannot explicitly check that
    # thepointset is of type Pointset, because Points.py imports this file
    # (utils.py), so check an attribute instead.
    try:
        thepointset.coordnames
    except AttributeError:
        raise TypeError, "Must pass Pointset to this function: use arrayToPointset first!"
    infodict_usedkeys = []
    for key, info in infodict.iteritems():
        if isinstance(info, str):
            infodict_usedkeys += [info]
        elif info == []:
            infodict[key] = copy.copy(thepointset.coordnames)
            infodict_usedkeys.extend(thepointset.coordnames)
        else:
            infodict_usedkeys += list(info)
    allnames = copy.copy(thepointset.coordnames)
    if thepointset._parameterized:
        allnames.append(thepointset.indepvarname)
    remlist = remain(infodict_usedkeys, allnames+range(len(allnames)))
    if remlist != []:
        print "Coords not found in pointset:", remlist
        raise ValueError, \
              "invalid keys in infodict - some not present in thepointset"
    assert isinstance(ext, str), "'ext' extension argument must be a string"
    if ext != '':
        if ext[0] != '.':
            ext = '.'+ext
    if append:
        assert varvaldir == 'col', ("append mode not supported for row"
                                     "format of data ordering")
        modestr = 'a'
    else:
        modestr = 'w'
    totlen = len(thepointset)
    if totlen == 0:
        raise ValueError, ("Pointset is empty")
    for fname, tup in infodict.iteritems():
        try:
            f = open(fname+ext, modestr)
        except IOError:
            print "There was a problem opening file "+fname+ext
            raise
        try:
            if isinstance(tup, str):
                try:
                    varray = thepointset[tup]
                except TypeError:
                    raise ValueError, "Invalid specification of coordinates"
            elif isinstance(tup, int):
                try:
                    varray = thepointset[:,tup].toarray()
                except TypeError:
                    raise ValueError, "Invalid specification of coordinates"
            elif type(tup) in [list, tuple]:
                if alltrue([type(ti)==str for ti in tup]):
                    thetup=list(tup)
                    if thepointset.indepvarname in tup:
                        tix = thetup.index(thepointset.indepvarname)
                        thetup.remove(thepointset.indepvarname)
                    try:
                        vlist = thepointset[thetup].toarray().tolist()
                    except TypeError:
                        raise ValueError, "Invalid specification of coordinates"
                    if len(thetup)==1:
                        vlist = [vlist]
                    if thepointset.indepvarname in tup:
                        vlist.insert(tix, thepointset.indepvararray.tolist())
                    varray = array(vlist)
                elif alltrue([type(ti)==int for ti in tup]):
                    try:
                        varray = thepointset[:,tup].toarray()
                    except TypeError:
                        raise ValueError, "Invalid specification of coordinates"
                else:
                    raise ValueError, "Invalid specification of coordinates"
            else:
                f.close()
                raise TypeError, \
                   "infodict values must be singletons or tuples/lists of strings or integers"
        except IOError:
            f.close()
            print "Problem writing to file"+fname+ext
            raise
        except KeyError:
            f.close()
            raise KeyError, ("Keys in infodict not found in pointset")
        if varvaldir == 'row':
            write_array(f, varray, separator, linesep,
                        precision, suppress_small, keep_open=0)
        else:
            write_array(f, transpose(varray), separator, linesep,
                        precision, suppress_small, keep_open=0)


# import text files containing data points into numeric arrays (not pointsets)
def importPointset(xFileName, t=None, indices=None, sep=" "):
    if indices is None:
        indices = []
    xFile = open(xFileName, 'r')
    xFileStrList = xFile.readlines()
    numpoints = len(xFileStrList)
    assert numpoints > 1, "Only 1 data point found in variables datafile"
    x_dummy_all = xFileStrList[0].rstrip("\n")
    x_dummy_vallist = filter(lambda s: s != '', x_dummy_all.split())
    tVals = zeros(numpoints, Float32)
    if t is None:
        get_t = 0
    elif isinstance(t, str):
        tFileName = t
        tFile = open(tFileName, 'r')
        tFileStrList = tFile.readlines()
        if len(tFileStrList) != numpoints:
            raise ValueError, ("Length of data and time files must be equal"
                           " -- are there any blank lines in the files?")
        get_t = 1
    elif type(t) in [list, type(array([1]))]:
        if len(t) != numpoints:
            raise ValueError, ("Length of data file and t array must be "
                       "equal -- are there any blank lines in the files?")
        tVals = t
        get_t = 0
    elif isinstance(t, int):
        # t represents column index to find time data in data file
        if t >= len(x_dummy_vallist) or t < 0:
            raise ValueError, "t index out of range"
        get_t = 2
    if indices == []:
        if get_t == 2:
            dim = len(x_dummy_vallist)-1
            indices = remain(range(0,dim+1),[t])
        else:
            dim = len(x_dummy_vallist)
            indices = range(0,dim)
    else:
        dim = len(indices)
        if get_t == 2:
            if t in indices:
                raise ValueError, ("You specified column "+str(t)+" as time "
                    "data, but you have specified it as a data column in "
                    "indices argument")
    xVals = zeros([numpoints, dim], Float32)
    for i in xrange(numpoints):
        vLine = xFileStrList[i].rstrip("\n")
        if vLine == '':
            continue
        vLineVals = filter(lambda s: s != '', vLine.split(sep))
        if get_t == 1:
            # Additional left strip of space char in case sep is different
            tLine = tFileStrList[i].rstrip("\n").lstrip(sep).lstrip(" ")
            if len(tLine.split(sep)) != 1:
                raise ValueError, ("Only one t value expected per line of"
                                   " datafile")
            if tLine == '':
                continue
            tVals[i] = float(tLine)
        elif get_t == 2:
            tVals[i] = float(vLineVals[t])
        try:
            xLineVals = [vLineVals[ix] for ix in indices]
        except IndexError:
            print "Valid indices were: 0-", len(vLineVals)-1
            raise
        if len(xLineVals) != dim:
            raise ValueError, ("Exactly "+str(dim)+" values expected per "
                               "line of datafile")
        xVals[i] = array([float(xstr) for xstr in xLineVals], Float32)
    xFile.close()
    if get_t == 1:
        tFile.close()
    xVals.transpose()
    if t is None:
        return xVals
    else:
        return {'t': tVals, 'vararray': xVals}


def intersect(a, b):
    """Find intersection of two lists, sequences, etc.
    Includes repetitions if they occur in the input lists."""
    return filter(lambda e : e in b, a)

def remain(a, b):
    """Find remainder of two lists, sequences, etc., after intersection.
    Includes repetitions if they occur in the input lists."""
    return filter(lambda e : e not in b, a)

def compareList(a, b):
    """Compare elements of lists, ignoring order (like sets)."""
    return len(intersect(a,b))==len(a)==len(b)


def makeDataDict(fieldnames, fieldvalues):
    """Zip arrays of field names and values into a dictionary.
    For instance, to use in Generator initialization arguments."""
    if type(fieldvalues) in [Array, NArray]:
        return dict(zip(fieldnames, [a.tolist() for a in fieldvalues]))
    else:
        return dict(zip(fieldnames, fieldvalues))
