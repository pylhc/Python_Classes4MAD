#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      14 June 2006
#$Revision: 1.0.8 $
"""
    Internal utilities.

    Robert Clewley, September 2005.
"""

# ----------------------------------------------------------------------------
# Version History
#
# v 1.0.8
#  Renamed sortDictByValue to sortedDictLists, to be consistent with similar
#    functions. Added sortedDictItems.
#
# v 1.0.7
#  Fixed bug in args class copy method.
#
# v 1.0.6
#  Added interp0d class for piecewise-constant interpolation.
#  Added NArray constant for numeric/scipy array class identification.
#
# v 1.0.5
#  Added Struct class.
#  Added filteredDict.
#  Modification to args class to provide __repr__ and info methods, and to
#    make dict(args()) work.
#  Fixed problems with diff.
#  Added arraymax, simplifyMatrixRepr utilities.
#
# v 1.0.4
#  Added getSuperClasses.
#  Added getSeqUnique.
#  Added insertInOrder.
#  Added DefaultDict.
#  Added sortDictByValue.
#  Renamed function 'deriv' to 'diff'.
#
# v 1.0.3
#  Updated compareBaseClass().
#  Added compareClassAndBases() to deal with a list of base classes and the
#    check that the first argument's type is actually the second argument.
#  Added className function and object2str conversion util for pretty printer.
#  Added args container class.
#  Removed parsing utilities to separate file.
#  Updated makeUniqueFn to provide option to add idstr to unique name.
#
# v 1.0.2
#  Adapted and improved function 'deriv' to accept scalar functions
#
# v 1.0.1
#  Added numerical derivative function
#
# ----------------------------------------------------------------------------

from errors import *

import sys, types
from scipy import Inf, NaN, atleast_1d, clip, less, greater, logical_or, \
        searchsorted, isfinite, shape, mat, sign, any, all, sometrue, alltrue
from scipy import integer as Numeric_Int
from scipy import array as Numeric_array
##from scipy import zeros as Numeric_zeros
##from scipy import ones as Numeric_ones
from numarray import array, Float, Int, Complex, Float64, Int32, Complex64, \
        Complex32, Float32, swapaxes, asarray, zeros, ones, \
        take, less_equal, putmask
import time
from copy import copy, deepcopy
import os
if os.name == 'nt':
    # slow object copying for you guys
    import fixedpickle as pickle
else:
    import cPickle as pickle

# ----------------------------------------------------------------------------
### EXPORTS

_classes = ['Verbose', 'interpclass', 'interp0d', 'interp1d', 'Utility',
            'args', 'Array', 'NArray', 'DefaultDict', 'Struct', 'pickle']

_mappings = ['_pytypename', '_pynametype', '_namefromtypecode',
             '_typefromtypecode', '_internaltypes', '_typefrompytype',
             '_pytypefromtype', '_nummaxmin', 'LargestInt32']

_functions = ['uniqueList', 'makeArrayIxMap', 'diff', 'className',
              'compareBaseClass', 'compareClassAndBases', 'timestamp',
              'makeUniqueFn', 'copyVarDict', 'concatStrDict',
              'invertMap', 'makeSeqUnique', 'insertInOrder',
              'sortedDictKeys', 'sortedDictValues', 'sortedDictItems',
              'sortedDictLists',
              'listid', 'idfn', 'noneFn', 'isincreasing', 'extent',
              'linearInterp', 'object2str', 'getSuperClasses',
              'filteredDict', 'arraymax', 'simplifyMatrixRepr']

_constants = ['Continuous', 'Discrete', 'targetLangs', 'SeqTypes']

__all__ = _functions + _mappings + _classes + _constants

# ----------------------------------------------------------------------------

# global reference for supported target languages
targetLangs = ['c', 'python', 'matlab'] #, 'xpp', 'dstool']

# useful reference to the numarray type
Array = type(array([1.0]))
NArray = type(Numeric_array([1.0]))

class Struct(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)
        
    def __repr__(self):
        attributes = [attr for attr in dir(self) if attr[0] != '_']
        return 'Struct(' + ', '.join(attributes) + ')'


class DefaultDict(dict):
    """Dictionary with a default value for unknown keys.
    
    Written by Peter Norvig."""
    def __init__(self, default):
        self.default = default

    def __getitem__(self, key):
        if key in self: return self.get(key)
        return self.setdefault(key, deepcopy(self.default))



class args(object):
    """Mapping object class for building arguments for class initialization
    calls. Treat as a dictionary.
    """

    def __init__(self, **kw):
        self.__dict__ = kw

    def __repr__(self):
        if len(self.__dict__) > 0:
            res = "Arguments:"
            for k, v in self.__dict__.iteritems():
                res += "\n %s : %s"%(k,str(v))
            return res
        else:
            return "No arguments defined"

    def info(self):
        print repr(self)

    __str__ = __repr__

    def values(self):
        return self.__dict__.values()

    def keys(self):
        return self.__dict__.keys()

    def items(self):
        return self.__dict__.items()

    def itervalues(self):
        return self.__dict__.itervalues()

    def iterkeys(self):
        return self.__dict__.iterkeys()

    def iteritems(self):
        return self.__dict__.iteritems()

    def __getitem__(self, k):
        return self.__dict__[k]

    def __setitem__(self, k, v):
        self.__dict__.__setitem__(k, v)

    def update(self, d):
        self.__dict__.update(d)

    def copy(self):
        return copy(self)

    def clear(self):
        self.__dict__.clear()

    def get(self, k, d=None):
        return self.__dict__.get(k, d)

    def has_key(self, k):
        return self.__dict__.has_key(k)

    def pop(self, k, d=None):
        return self.__dict__.pop(k, d)

    def popitem(self):
        raise NotImplementedError

    def __contains__(self, v):
        return self.__dict__.__contains__(v)

    def fromkeys(self, S, v=None):
        raise NotImplementedError

    def setdefault(self, d):
        raise NotImplementedError

    def __delitem__(self, k):
        del self.__dict__[k]

    def __cmp__(self, other):
        return self.__dict__ == other

    def __eq__(self, other):
        return self.__dict__ == other

    def __ne__(self, other):
        return self.__dict__ != other

    def __gt__(self, other):
        return self.__dict__ > other

    def __ge__(self, other):
        return self.__dict__ >= other

    def __lt__(self, other):
        return self.__dict__ < other

    def __le__(self, other):
        return self.__dict__ <= other

    def __len__(self):
        return len(self.__dict__)

    def __iter__(self):
        return iter(self.__dict__)

##    def copy(self):
##        return args(**self.__dict__)
##
##    def __copy__(self):
##        return args(**self.__dict__)
##
##    __deepcopy__ = __copy__


## Internally used mappings

LargestInt32 = 2147483647

_pytypename = {Float: 'float', Float64: 'float',
               Int: 'int', Int32: 'int',
               Complex: 'complex', Complex64: 'complex'}

_pynametype = {'float': Float, 'int': Int, 'complex': Complex}

_nummaxmin = {Float: [-Inf, Inf],
             Int: [-LargestInt32-1, LargestInt32],
             Complex: [-Inf-Inf*1.0j, Inf+Inf*1.0j],
             Float64: [-Inf, Inf],
             Int32: [-LargestInt32-1, LargestInt32],
             Complex64: [-Inf-Inf*1.0j, Inf+Inf*1.0j]
             }

_typefrompytype = {float: Float, int: Int, complex: Complex}
_pytypefromtype = {Float: float, Int: int, Complex: complex}

_namefromtypecode = {'i' : 'Int32',
                 'l' : 'Int32',   # should be long (?)
                 'f' : 'Float32',
                 'd' : 'Float64',
                 'F' : 'Complex32',
                 'D' : 'Complex64',
                 }
_typefromtypecode = {'i' : Int32,
                 'l' : Int32,   # should be long (?)
                 'f' : Float32,
                 'd' : Float64,
                 'F' : Complex32,
                 'D' : Complex64,
                 }

_typecodefromtype = {Int32: 'l',
                     Complex64: 'D',
                     Float32: 'f',
                     Complex32: 'F',
                     Float64: 'd'
                     }

# equivalents of Python internal types (as an unordered list)
_internaltypes = [Float, Complex, Int, Float64, Complex64, Int32]


SeqTypes = [list, tuple, type(array([1]))]


## ------------------------------------------------------------------

## Internally used functions

def filteredDict(d, keys, neg=False):
    out_d = {}
    if neg:
	out_keys = remain(d.keys(), keys)
    else:
	out_keys = keys
    for k in out_keys:
	try:
	    out_d[k] = d[k]
	except KeyError:
	    pass
    return out_d


def concatStrDict(d, order=[]):
    """Concatenates all entries of a dictionary (assumed to be
    lists of strings), in optionally specified order."""
    retstr = ''
    if d != {}:
        if order == []:
            order = d.keys()
        for key in order:
            itemlist = d[key]
            for strlist in itemlist:
                retstr += ''.join(strlist)
    return retstr


def copyVarDict(vardict):
    """Copy dictionary of Variable objects."""
    return dict(zip(sortedDictKeys(vardict), [copy(v) for v in \
                                              sortedDictValues(vardict)]))


def insertInOrder(sourcelist, inslist, return_ixs=False):
    """Insert elements of inslist into sourcelist, sorting these
    lists in case they are not already in increasing order.

    The function will not create duplicate entries in the list, and will
    change neither the first or last entries of the list.

    If sourcelist is an array, an array is returned."""
    try:
        sorted_inslist = inslist.tolist()
    except AttributeError:
        sorted_inslist = copy(inslist)
    sorted_inslist.sort()
    try:
        sorted_sourcelist = sourcelist.tolist()
        was_array = True
    except AttributeError:
        sorted_sourcelist = copy(sourcelist)
        was_array = False
    sorted_sourcelist.sort()
    tix = 0
    # optimize by having separate versions of loop
    if return_ixs:
        ins_ixs = []
        for t in sorted_inslist:
            tcond = less_equal(sorted_sourcelist[tix:], t).tolist()
            try:
                tix = tcond.index(0) + tix  # lowest index for elt > t
                if sorted_sourcelist[tix-1] != t and tix >= 0:
                    sorted_sourcelist.insert(tix, t)
                    ins_ixs.append(tix)
            except ValueError:
                # no 0 value in tcond, so t might be equal to the final value
                pass
        if was_array:
            return array(sorted_sourcelist), ins_ixs
        else:
            return sorted_sourcelist, ins_ixs
    else:
        for t in sorted_inslist:
            tcond = less_equal(sorted_sourcelist[tix:], t).tolist()
            try:
                tix = tcond.index(0) + tix  # lowest index for elt > t
                if sorted_sourcelist[tix-1] != t and tix >= 0:
                    sorted_sourcelist.insert(tix, t)
            except ValueError:
                # no 0 value in tcond, so t might be equal to the final value
                pass
        if was_array:
            return array(sorted_sourcelist)
        else:
            return sorted_sourcelist


def arraymax(a1,a2,t=Float):
    """Element-wise comparison of maximum values for two arrays."""
    o=[]
    try:
        for x, y in zip(a1,a2):
            o.append(max(x,y))
    except TypeError:
        print "Problem with type of arguments in arraymax:"
        print "Received a1 =", a1
        print "         a2 =", a2
        raise
    return array(o,t)

def simplifyMatrixRepr(m):
    """Convert Scipy Matrix object to a compact array
    representation or numeric value."""
    ma=array(m)
    l = len(shape(ma))
    if l == 0:
	return m
    elif l>0 and shape(ma)[0] == 1:
        return simplifyMatrixRepr(ma[0])
    elif l>1 and shape(ma)[1] == 1:
        return simplifyMatrixRepr(ma[:,0])
    else:
        return ma


def diff(func, x0, vars=None, axes=None, dir=1, dx=None):
    """Numerical 1st derivative of R^N -> R^M function about x0 by finite
    differences.

    dir=1 uses finite forward difference.
    dir=-1 uses finite backward difference.
    List valued dx rescales finite differencing in each axis separately.
    vars argument specifies which elements of x0 are to be treated as
      variables for the purposes of taking the Jacobian.
    If axes argument is unused or set to be all axes, the Jacobian of the
      function evaluated at x0 with respect to the variables is returned,
      otherwise a sub-matrix of it is returned.
    If dx is not given an appropriate step size is chosen (proportional to
      sqrt(machine precision)).
    """

    if type(x0) in SeqTypes:
        x0type = 'seq'
    elif type(x0) in [int, float]:
        x0type = 'num'
    else:
        # assume Point type (will get errors if not!)
        x0type = 'point'
    if vars is None:
        if x0type == 'seq':
            dim = len(x0)
            vars = range(dim)
        elif x0type == 'num':
            dim = 1
            vars = [0]
        else:
            # Point type
            dim = x0.dimension
            vars = x0.coordnames
    else:
        assert type(vars) in SeqTypes, "vars argument must be a sequence type"
        if x0type in ['seq', 'num']:
            assert all(vars>=0), \
                    "vars argument must hold non-negative integers"
        else:
            assert all([isinstance(vars[i], str) \
                 for i in range(len(vars))]), "vars argument must hold strings"
        dim = len(vars)
    fx0 = func(x0)
    try:
        # ensure fx0 is a vector or at least only a D x 1 matrix
        assert shape(fx0)[1] == 1
    except IndexError:
        # then shape is of form (D,) and that's fine
        pass
    except AssertionError:
        print "fx0 shape is", shape(fx0)
        print fx0
        raise ValueError("function should return an N-vector or N x 1 matrix,"
                 " but it returned a matrix with shape %s" % str(shape(fx0)))
    if type(fx0) in [int, float]:
        dimf = 1
    else:
        try:
            dimf = shape(fx0)[0]
        except IndexError:
            dimf = 1
    if axes is None:
        if x0type in ['seq', 'num']:
            try:
                axes = range(shape(fx0)[0])
            except IndexError:
                # then singleton (scalar) was returned
                axes = [0]
        else:
            axes = fx0.coordnames
    else:
        assert type(axes) in SeqTypes, "axes argument must be a sequence type"
        if x0type in ['seq', 'num']:
            assert all(axes>=0), \
                   "axes argument must hold non-negative integers"
        else:
            assert all([isinstance(axes[i], str) \
                 for i in range(len(axes))]), "axes argument must hold strings"
    if dx is None:
        # Use formula for h = sqrt(eps)*|x|, eps = machine precision,
        # which is ~ 1e-16 for double precision floating point.
        # Avoid abs(x0) being smaller than ones, but extra term
        # + x0==zeros(dim) avoids case when x0 is zeros so that sign(x0) is
        # all zeros making dx = zeros ( --> division by zero !) ... need
        # to coerce the boolean array to Float so that it can be added to
        # sign(x0)!
        if x0type == 'seq':
            dx = 1e-8*arraymax(abs(x0[vars]),ones(dim))*(sign(x0[vars]) \
                                    + array(x0[vars]==zeros(dim),'Float'))
        elif x0type == 'num':
            dx = 1e-8*max(abs(x0),1)*(sign(x0)+int(x0==0))
        else:
            x0a = x0[vars].toarray()
            dx = dict(zip(vars,
                          1e-8*arraymax(abs(x0a),ones(dim))*(sign(x0a) + \
                                       array(x0a==zeros(dim),'Float'))))
    else:
        if isinstance(dx, float):
            assert dx > 0, "dx scaling array must be strictly positive"
            if x0type in ['seq', 'num']:
                dx = ones(dim)*dx
            else:
                dx = dict(zip(vars,ones(dim)*dx))
        else:
            assert len(dx) == len(vars), \
                   "dx scaling array has length mismatch with vars"
            if x0type in ['seq', 'num']:
                assert all(dx>0), \
                       "dx scaling array must be strictly positive"
            else:
                assert all(dx.values()>0), \
                       "dx scaling array must be strictly positive"
    assert dir==1 or dir==-1, "Direction code must be -1 or 1"
    dim_mat = len(axes)
    df = zeros([dim_mat,dim], 'Float')
    if x0type == 'seq':
        for i in range(dim):
            vix = vars[i]
            try:
                # for numarrays (otherwise copy returns a regular 'array'!)
                x0_d = x0.copy()
            except AttributeError:
                x0_d = copy(x0)
            x0_d[vix] += dir * dx[vix]
            fx0_d = func(x0_d)
            if dim_mat > 1:
                fx0_d_v = array([fx0_d[n] for n in axes])
                fx0_v = array([fx0[n] for n in axes])
            else:
                if dimf > 1:
                    fx0_d_v = fx0_d[axes[0]]
                    fx0_v = fx0[axes[0]]
                else:
                    fx0_d_v = fx0_d
                    fx0_v = fx0
            df[:,vix] = dir*(fx0_d_v - fx0_v)/dx[vix]
        return mat(df)
    elif x0type == 'num':
        x0_d = x0 + dir*dx
        fx0_d = func(x0_d)
        df = dir*(fx0_d - fx0)/dx
        return df
    else:
        # Point type
        for i in range(dim):
            vname = vars[i]
            x0_d = copy(x0)
            x0_d[vname] = x0_d(vname) + dir * dx[vname]
            fx0_d = func(x0_d)(axes)
            fx0_v = fx0(axes)
            df[:,i] = dir*(fx0_d - fx0_v)/dx[vname]
        return mat(df)


def linearInterp(y0, ygoal, y1, x0, x1):
    """Linearly interpolate data points."""
    return ( x1 * (ygoal - y0) + x0 * ( y1 - ygoal) ) / (y1 - y0)


def makeUniqueFn(fstr, tdigits=7, idstr=None):
    """Add unique ID to function names.

    Used when functions are executed in global namespace to avoid name
    clashes, and need to be distinguished when DS objects are copied."""
    # check for syntax errors
    try:
        code = compile(fstr, 'test', 'exec')
    except:
        print " Cannot make unique function because of a syntax (or other) error " \
              "in supplied code:\n"
        print fstr
        raise
    bracepos = fstr.index("(")
    if idstr is None:
        idstr_insert = ""
    else:
        idstr_insert = "_" + idstr
    fname = fstr[4:bracepos] + idstr_insert + "_" + timestamp(tdigits)
    fstr_new = "def " + fname + fstr[bracepos:]
    return (fstr_new, fname)


# return a unique timestamp string for the session. useful for ensuring
# unique function identifiers, etc.
def timestamp(tdigits=8):
    c = str(time.clock())
    return c.replace(".", "")[:tdigits+1]


def uniqueList(objlist):
    """check that list contains items only once"""
    if len(objlist):
        return reduce(bool.__and__, [objlist.count(obj) == 1 \
                          for obj in objlist])
    else:
        return True

# makeSeqUnique adapted from code by Raymond Hettinger, 2002
def makeSeqUnique(seq, asarray=False, preserveorder=True):
    """Return a 1D sequence that only contains the unique values in seq."""
    set = {}
    if preserveorder:
        if asarray:
            return array([set.setdefault(e,e) for e in seq if e not in set])
        else:
            return [set.setdefault(e,e) for e in seq if e not in set]


def getSuperClasses(obj, limitClasses=None):
    if limitClasses == None:
        limitClassNames = ['object']
    elif isinstance(limitClasses, list):
        limitClassNames = [className(lc) for lc in limitClasses]
    else:
        # singleton class
        limitClassNames = [className(limitClasses)]
    # ensure "object" safety net is present
    if 'object' not in limitClassNames:
        limitClassNames.append('object')
    search_obj = [obj.__class__]
    sclasses = [className(search_obj[0])]
    # don't start while loop if obj is already of a type in limitClasses
    done = (sclasses[0] in limitClassNames)
    c = 0
    while not done and c < 10:
        c += 1
        search_temp = []
        for so in search_obj:
            search_temp.extend(list(so.__bases__))
        search_obj = search_temp
        for b in search_obj:
            sclass = className(b)
            done = sclass in limitClassNames
            if done:
                break
            else:
                sclasses.append(sclass)
    return sclasses
    

def className(obj, addPrefix=False):
    """Return human-readable string of class name."""
    if isinstance(obj, type):
        class_str = obj.__name__
        if addPrefix:
            classPrefix = "Class "
        else:
            classPrefix = ""
    elif isinstance(obj, types.ModuleType):
        class_str = obj.__name__
        if addPrefix:
            classPrefix = "Module "
        else:
            classPrefix = ""
    else:
        try:
            class_str = obj.__class__.__name__
        except AttributeError:
            class_str = str(type(obj))
        classPrefix = ""
    return classPrefix + class_str


def object2str(x):
    """Convert occurrences of types / classes,
    to pretty-printable strings."""
    try:
        if type(x) in [types.InstanceType, types.TypeType]:
            return className(x, True)
        elif isinstance(x, list):
            # search through any iterable parts (that aren't strings)
            rx = "["
            if len(x)>0:
                for o in x:
                    rx += object2str(o) + ", "
                return rx[:-2]+"]"
            else:
                return rx+"]"
        elif isinstance(x, tuple):
            rx = "("
            if len(x)>0:
                for o in x:
                    rx += object2str(o) + ", "
                return rx[:-2]+")"
            else:
                return rx+")"
        elif isinstance(x, dict):
            rx = "{"
            if len(x)>0:
                for k, o in x.iteritems():
                    rx += object2str(k) + ": " + object2str(o) + ", "
                return rx[:-2]+"}"
            else:
                return rx+"}"
        elif isinstance(x, str):
            # this removes extraneous single quotes around dict keys, for instance
            return x
        else:
            return repr(x)
    except:
        raise TypeError, "object2str cannot format this object type"


# compareBaseClass is slightly adapted from AUTO 2000 code's findBaseClass:
# Booleans replace `0`, `1`, we allow objects to be passed as well as classes,
#  and we compare class name strings not the class types themselves.
#  The class types can show different roots when they originate from
#  different parts of the PyDSTool package -- it might be a bug.
#  e.g. baseClass might be Generator.Generator, but here this type will be
#  <class 'PyDSTool.Generator.baseclasses.Generator'>
#  and input.__class__ will boil down to
#  <class 'Generator.baseclasses.Generator'>
#  even though these classes are identical (constructed from the same class
#  in the same module!)
def compareBaseClass(input, baseClass):
    try:
        test_has_bases = input.__bases__
        # if got this far then we were passed a class
        inputClass = input
    except AttributeError:
        try:
            inputClass = input.__class__
        except AttributeError:
            # try this as a last resort
            return className(input) == className(baseClass)
    for base in inputClass.__bases__:
        if base.__name__ == baseClass.__name__:
            return True
        else:
            if compareBaseClass(base, baseClass):
                return True
            # else continue looking
    return False


def compareClassAndBases(input, arg):
    if isinstance(arg, list):
        return bool(any([compareClassAndBases(input, a) for a in arg]))
    else:
        if type(input) is arg:
            return True
        elif input is arg:
            return True
        else:
            return compareBaseClass(input, arg)


# little utility function to wrap value as a singleton list
def listid(val):
    return [val]


# the identity function
def idfn(val):
    return copy(val)


# utility function representing a "none" function
def noneFn(x):
    return None


# returns the mapping from the entries in an array or list to their indices
def makeArrayIxMap(a):
    if isinstance(a, type(array([1]))):
        return dict(zip(a.tolist(), range(len(a))))
    else:
        return dict(zip(a, range(len(a))))


# invert an index mapping or other form of mapping
def invertMap(themap):
    """invert an index mapping or other form of mapping.

    invertMap() returns a dictionary, regardless of whether the
    argument is a dictionary or a list."""
    if isinstance(themap, dict):
        return dict(map(lambda (k,v): (v,k), themap.items()))
    elif isinstance(themap, list):
        # input domain is the position index
        return dict(zip(themap, range(len(themap))))
    elif isinstance(themap, Array):
        # input domain is the position index
        return dict(zip(themap.tolist(), range(len(themap))))
    else:
        raise TypeError, "Unsupported type for map"


# check that a sequence type is in increasing order
def isincreasing(theseq):
    try:
        v_old = theseq[0]
    except IndexError:
        raise ValueError, "Problem with sequence passed to isincreasing(). Is it empty?"
    for v in theseq[1:]:
        if v <= v_old:
            print "sequence not increasing: ", v_old, v
            return False
        v_old = v
    return True


# min and max values of a dataset
def extent(data):
    minval = min(data)
    maxval = max(data)
    if minval == maxval:
        return minval
    else:
        return [minval, maxval]


def sortedDictValues(d, onlykeys=[], reverse=False):
    """Return list of values from a dictionary in order of sorted key list.

    Adapted from original function by Alex Martelli:
     added filtering of keys."""
    if onlykeys != []:
        keys = intersect(d.keys(), onlykeys)
    else:
        keys = d.keys()
    keys.sort()
    if reverse:
        keys.reverse()
    return map(d.get, keys)

def sortedDictKeys(d, onlykeys=[], reverse=False):
    """Return sorted list of keys from a dictionary.

    Adapted from original function by Alex Martelli:
     added filtering of keys."""
    if onlykeys != []:
        keys = intersect(d.keys(), onlykeys)
    else:
        keys = d.keys()
    keys.sort()
    if reverse:
        keys.reverse()
    return keys

def sortedDictLists(d, byvalue=True, reverse=False):
    """Return (key list, value list) pair from a dictionary,
    sorted by value (default) or key.
    Adapted from an original function by Duncan Booth.
    """
    if byvalue:
        i = [(val, key) for (key, val) in d.items()]
        i.sort()
        if reverse:
            i.reverse()
        rvals = [val for (val, key) in i]
        rkeys = [key for (val, key) in i]        
    else:
        # by key
        i = d.items()
        i.sort()
        if reverse:
            i.reverse()
        rvals = [val for (key, val) in i]
        rkeys = [key for (key, val) in i]
    return (rkeys, rvals)

def sortedDictItems(d, byvalue=True, onlykeys=[]):
    """Return list of (key, value) pairs of a dictionary,
    sorted by value (default) or key.
    Adapted from an original function by Duncan Booth.
    """
    ks, vs = sortedDictLists(d, byvalue, onlykeys)
    return zip(ks,vs)
    
# ----------------------------------------------------------------------

## private versions of these utils (cannot import them from utils!)

# find intersection of two lists, sequences, etc.
def intersect(a, b):
    return filter(lambda e : e in b, a)


# find remainder of two lists, sequences, etc., after intersection
def remain(a, b):
    return filter(lambda e : e not in b, a)


# ----------------------------------------------------------------------

class Utility(object):
    """
    Utility abstract class for manipulating and analyzing dynamical systems.

    Robert Clewley, March 2005.

Subclasses of Utility could include such things as continuation tools,
dimension reduction tools, parameter estimation tools.
"""
    pass


# ----------------------------------------------------------------------
# This section adapted from scipy.interpolate
# ----------------------------------------------------------------------

_take = take
def take(a,indices,axis=0):
    x = asarray(a); y = asarray(indices)
    if shape(x) == (): x = x.flat
    if shape(y) == (): y = y.flat
    return _take(x,y,axis)

_sometrue = sometrue
def sometrue(a,axis=0):
    x = asarray(a)
    if shape(x) == (): x = x.flat
    return _sometrue(x)

def reduce_sometrue(a):
    all = a
    while len(shape(all)) > 1:
        all = sometrue(all)
    return all


class interpclass(object):
    """Abstract class for interpolators."""
    interp_axis = -1    # used to set which is default interpolation
                        # axis.  DO NOT CHANGE OR CODE WILL BREAK.


class interp0d(interpclass):
    """Design of this class based on SciPy's interp1d"""
    
    def __init__(self, x, y, axis=-1, makecopy=0, bounds_error=1,
                 fill_value=None):
        """Initialize a piecewise-constant interpolation class

        Description:
          x and y are arrays of values used to approximate some function f:
            y = f(x)
          This class returns a function whose call method uses piecewise-
          constant interpolation to find the value of new points.

        Inputs:
            x -- a 1d array of monotonically increasing real values.
                 x cannot include duplicate values. (otherwise f is
                 overspecified)
            y -- an nd array of real values.  y's length along the
                 interpolation axis must be equal to the length
                 of x.
            axis -- specifies the axis of y along which to 
                    interpolate. Interpolation defaults to the last
                    axis of y.  (default: -1)
            makecopy -- If 1, the class makes internal copies of x and y.
                    If 0, references to x and y are used. The default 
                    is to copy. (default: 0)
            bounds_error -- If 1, an error is thrown any time interpolation
                            is attempted on a value outside of the range
                            of x (where extrapolation is necessary).
                            If 0, out of bounds values are assigned the
                            NaN (#INF) value.  By default, an error is
                            raised, although this is prone to change.
                            (default: 1)
        """
        self.datapoints = (array(x, Float), array(y, Float))   # RHC -- for access from PyDSTool
        self.type = Float   # RHC -- for access from PyDSTool
        self.axis = axis
        self.makecopy = makecopy   # RHC -- renamed from copy to avoid nameclash
        self.bounds_error = bounds_error
        if fill_value is None:
            self.fill_value = NaN   # RHC -- was:   array(0.0) / array(0.0)
        else:
            self.fill_value = fill_value

        # Check that both x and y are at least 1 dimensional.
        if len(shape(x)) == 0 or len(shape(y)) == 0:
            raise ValueError, "x and y arrays must have at least one dimension."  
        # make a "view" of the y array that is rotated to the
        # interpolation axis.
        oriented_x = x
        oriented_y = swapaxes(y,self.interp_axis,axis)
        interp_axis = self.interp_axis
        len_x,len_y = shape(oriented_x)[interp_axis], \
                            shape(oriented_y)[interp_axis]
        if len_x != len_y:
            raise ValueError, "x and y arrays must be equal in length along "\
                              "interpolation axis."
        if len_x < 2 or len_y < 2:
            raise ValueError, "x and y arrays must have more than 1 entry"            
        self.x = array(oriented_x,copy=self.makecopy)
        self.y = array(oriented_y,copy=self.makecopy)


    def __call__(self,x_new):
        """Find piecewise-constant interpolated y_new = <name>(x_new).

        Inputs:
          x_new -- New independent variables.

        Outputs:
          y_new -- Piecewise-constant interpolated values corresponding to x_new.
        """
        # 1. Handle values in x_new that are outside of x.  Throw error,
        #    or return a list of mask array indicating the outofbounds values.
        #    The behavior is set by the bounds_error variable.
        ## RHC -- was   x_new = atleast_1d(x_new)
        x_new_1d = atleast_1d(x_new)
        out_of_bounds = self._check_bounds(x_new_1d)
        # 2. Find where in the orignal data, the values to interpolate
        #    would be inserted.
        #    Note: If x_new[n] = x[m], then m is returned by searchsorted.
        x_new_indices = searchsorted(self.x,x_new_1d)
        # 3. Clip x_new_indices so that they are within the range of 
        #    self.x indices and at least 1.  Removes mis-interpolation
        #    of x_new[n] = x[0]
        # RHC -- changed Int to Numeric_Int to avoid name clash with numarray
        x_new_indices = clip(x_new_indices,1,len(self.x)-1).astype(Numeric_Int)
        # 4. Calculate the region that each x_new value falls in.
        lo = x_new_indices - 1; hi = x_new_indices
        
        # !! take() should default to the last axis (IMHO) and remove
        # !! the extra argument.
        # 5. Calculate the actual value for each entry in x_new.
        y_lo = take(self.y,lo,axis=self.interp_axis)
        y_hi = take(self.y,hi,axis=self.interp_axis)
        y_new = (y_lo+y_hi)/2.
        # 6. Fill any values that were out of bounds with NaN
        # !! Need to think about how to do this efficiently for 
        # !! mutli-dimensional Cases.
        yshape = y_new.shape
        y_new = y_new.flat
        new_shape = list(yshape)
        new_shape[self.interp_axis] = 1
        sec_shape = [1]*len(new_shape)
        sec_shape[self.interp_axis] = len(out_of_bounds)
        out_of_bounds.shape = sec_shape
        new_out = ones(new_shape)*out_of_bounds
        putmask(y_new, new_out.flat, self.fill_value)
        y_new.shape = yshape
        # Rotate the values of y_new back so that they correspond to the
        # correct x_new values.
        result = swapaxes(y_new,self.interp_axis,self.axis)
        try:
            len(x_new)
            return result
        except TypeError:
            return result[0]
        return result


    def _check_bounds(self,x_new):
        # If self.bounds_error = 1, we raise an error if any x_new values
        # fall outside the range of x.  Otherwise, we return an array indicating
        # which values are outside the boundary region.  
        # !! Needs some work for multi-dimensional x !!
        below_bounds = less(x_new,self.x[0])
        above_bounds = greater(x_new,self.x[-1])
        #  Note: sometrue has been redefined to handle length 0 arrays
        # !! Could provide more information about which values are out of bounds
        # RHC -- Changed these ValueErrors to PyDSTool_BoundsErrors
        if self.bounds_error and any(sometrue(below_bounds)):
##            print "Input:", x_new
##            print "Bound:", self.x[0]
##            print "Difference input - bound:", x_new-self.x[0]
            raise PyDSTool_BoundsError, " A value in x_new is below the"\
                              " interpolation range."
        if self.bounds_error and any(sometrue(above_bounds)):
##            print "Input:", x_new
##            print "Bound:", self.x[-1]
##            print "Difference input - bound:", x_new-self.x[-1]
            raise PyDSTool_BoundsError, " A value in x_new is above the"\
                              " interpolation range."
        # !! Should we emit a warning if some values are out of bounds.
        # !! matlab does not.
        out_of_bounds = logical_or(below_bounds,above_bounds)
        return out_of_bounds


    # RHC added
    def __getstate__(self):
        d = copy(self.__dict__)
        # remove reference to Cfunc self.type
        d['type'] = _pytypename[self.type]
        return d

    # RHC added
    def __setstate__(self, state):
        self.__dict__.update(state)
        # reinstate Cfunc self.type
        self.type = _pynametype[self.type]

        

class interp1d(interpclass):    # RHC -- made this a new-style Python class
    def __init__(self, x, y, kind='linear', axis=-1,
                 makecopy = 0, bounds_error=1, fill_value=None):
        """Initialize a 1d piecewise-linear interpolation class

        Description:
          x and y are arrays of values used to approximate some function f:
            y = f(x)
          This class returns a function whose call method uses linear
          interpolation to find the value of new points.

        Inputs:
            x -- a 1d array of monotonically increasing real values.
                 x cannot include duplicate values. (otherwise f is
                 overspecified)
            y -- an nd array of real values.  y's length along the
                 interpolation axis must be equal to the length
                 of x.
            kind -- specify the kind of interpolation: 'nearest', 'linear',
                    'cubic', or 'spline'
            axis -- specifies the axis of y along which to 
                    interpolate. Interpolation defaults to the last
                    axis of y.  (default: -1)
            makecopy -- If 1, the class makes internal copies of x and y.
                    If 0, references to x and y are used. The default 
                    is to copy. (default: 0)
            bounds_error -- If 1, an error is thrown any time interpolation
                            is attempted on a value outside of the range
                            of x (where extrapolation is necessary).
                            If 0, out of bounds values are assigned the
                            NaN (#INF) value.  By default, an error is
                            raised, although this is prone to change.
                            (default: 1)
        """
        self.datapoints = (array(x, Float), array(y, Float))   # RHC -- for access from PyDSTool
        self.type = Float   # RHC -- for access from PyDSTool
        self.axis = axis
        self.makecopy = makecopy   # RHC -- renamed from copy to avoid nameclash
        self.bounds_error = bounds_error
        if fill_value is None:
            self.fill_value = NaN   # RHC -- was:   array(0.0) / array(0.0)
        else:
            self.fill_value = fill_value

        if kind != 'linear':
            raise NotImplementedError, "Only linear supported for now. Use fitpack routines for other types."

        # Check that both x and y are at least 1 dimensional.
        if len(shape(x)) == 0 or len(shape(y)) == 0:
            raise ValueError, "x and y arrays must have at least one dimension."  
        # make a "view" of the y array that is rotated to the
        # interpolation axis.
        oriented_x = x
        oriented_y = swapaxes(y,self.interp_axis,axis)
        interp_axis = self.interp_axis
        len_x,len_y = shape(oriented_x)[interp_axis], \
                            shape(oriented_y)[interp_axis]
        if len_x != len_y:
            raise ValueError, "x and y arrays must be equal in length along "\
                              "interpolation axis."
        if len_x < 2 or len_y < 2:
            raise ValueError, "x and y arrays must have more than 1 entry"            
        self.x = array(oriented_x,copy=self.makecopy)
        self.y = array(oriented_y,copy=self.makecopy)


    def __call__(self,x_new):
        """Find linearly interpolated y_new = <name>(x_new).

        Inputs:
          x_new -- New independent variables.

        Outputs:
          y_new -- Linearly interpolated values corresponding to x_new.
        """
        # 1. Handle values in x_new that are outside of x.  Throw error,
        #    or return a list of mask array indicating the outofbounds values.
        #    The behavior is set by the bounds_error variable.
        ## RHC -- was   x_new = atleast_1d(x_new)
        x_new_1d = atleast_1d(x_new)
        out_of_bounds = self._check_bounds(x_new_1d)
        # 2. Find where in the orignal data, the values to interpolate
        #    would be inserted.
        #    Note: If x_new[n] = x[m], then m is returned by searchsorted.
        x_new_indices = searchsorted(self.x,x_new_1d)
        # 3. Clip x_new_indices so that they are within the range of 
        #    self.x indices and at least 1.  Removes mis-interpolation
        #    of x_new[n] = x[0]
        # RHC -- changed Int to Numeric_Int to avoid name clash with numarray
        x_new_indices = clip(x_new_indices,1,len(self.x)-1).astype(Numeric_Int)
        # 4. Calculate the slope of regions that each x_new value falls in.
        lo = x_new_indices - 1; hi = x_new_indices
        
        # !! take() should default to the last axis (IMHO) and remove
        # !! the extra argument.
        x_lo = take(self.x,lo,axis=self.interp_axis)
        x_hi = take(self.x,hi,axis=self.interp_axis)
        y_lo = take(self.y,lo,axis=self.interp_axis)
        y_hi = take(self.y,hi,axis=self.interp_axis)
        slope = (y_hi-y_lo)/(x_hi-x_lo)
        # 5. Calculate the actual value for each entry in x_new.
        y_new = slope*(x_new_1d-x_lo) + y_lo 
        # 6. Fill any values that were out of bounds with NaN
        # !! Need to think about how to do this efficiently for 
        # !! mutli-dimensional Cases.
        yshape = y_new.shape
        y_new = y_new.flat
        new_shape = list(yshape)
        new_shape[self.interp_axis] = 1
        sec_shape = [1]*len(new_shape)
        sec_shape[self.interp_axis] = len(out_of_bounds)
        out_of_bounds.shape = sec_shape
        new_out = ones(new_shape)*out_of_bounds
        putmask(y_new, new_out.flat, self.fill_value)
        y_new.shape = yshape
        # Rotate the values of y_new back so that they correspond to the
        # correct x_new values.
        result = swapaxes(y_new,self.interp_axis,self.axis)
        try:
            len(x_new)
            return result
        except TypeError:
            return result[0]
        return result


    def _check_bounds(self,x_new):
        # If self.bounds_error = 1, we raise an error if any x_new values
        # fall outside the range of x.  Otherwise, we return an array indicating
        # which values are outside the boundary region.  
        # !! Needs some work for multi-dimensional x !!
        below_bounds = less(x_new,self.x[0])
        above_bounds = greater(x_new,self.x[-1])
        #  Note: sometrue has been redefined to handle length 0 arrays
        # !! Could provide more information about which values are out of bounds
        # RHC -- Changed these ValueErrors to PyDSTool_BoundsErrors
        if self.bounds_error and any(sometrue(below_bounds)):
##            print "Input:", x_new
##            print "Bound:", self.x[0]
##            print "Difference input - bound:", x_new-self.x[0]
            raise PyDSTool_BoundsError, " A value in x_new is below the"\
                              " interpolation range."
        if self.bounds_error and any(sometrue(above_bounds)):
##            print "Input:", x_new
##            print "Bound:", self.x[-1]
##            print "Difference input - bound:", x_new-self.x[-1]
            raise PyDSTool_BoundsError, " A value in x_new is above the"\
                              " interpolation range."
        # !! Should we emit a warning if some values are out of bounds.
        # !! matlab does not.
        out_of_bounds = logical_or(below_bounds,above_bounds)
        return out_of_bounds


    # RHC added
    def __getstate__(self):
        d = copy(self.__dict__)
        # remove reference to Cfunc self.type
        d['type'] = _pytypename[self.type]
        return d

    # RHC added
    def __setstate__(self, state):
        self.__dict__.update(state)
        # reinstate Cfunc self.type
        self.type = _pynametype[self.type]



class DomainType(object):
    def __repr__(self):
        raise NotImplementedError

    def __str__(self):
        raise NotImplementedError

    def __eq__(self, d):
        return self.__class__ == d.__class__


class ContinuousDomain(DomainType):
    def __repr__(self):
        return "Continuous Domain"

    __str__ = __repr__


class DiscreteDomain(DomainType):
    def __repr__(self):
        return "Discrete Domain"

    __str__ = __repr__


# treat these as "constants" as they are empty
global Continuous, Discrete
Continuous = ContinuousDomain()
Discrete = DiscreteDomain()



#-----------------------------------------------------------------------------
# The following code, in particular the Verbose class, was written by
# John D. Hunter as part of the front end of MatplotLib.
# (See matplotlib.sourceforge.net/license.html for details.)
#
# Copyright (c) 2002-2004 John D. Hunter; All Rights Reserved
#-----------------------------------------------------------------------------


class Verbose:
    """
    A class to handle reporting.  Set the fileo attribute to any file
    instance to handle the output.  Default is sys.stdout
    """
    levels = ('silent', 'error', 'helpful', 'debug', 'debug-annoying')
    vald = dict( [(level, i) for i,level in enumerate(levels)])

    # parse the verbosity from the command line; flags look like
    # --verbose-error or --verbose-helpful
    _commandLineVerbose = None


    for arg in sys.argv[1:]:
        if not arg.startswith('--verbose-'): continue
        _commandLineVerbose = arg[10:]

    def __init__(self, level):
        self.setLevel(level)
        self.fileo = sys.stdout
        self.erro = sys.stderr

    def setLevel(self, level):
        'set the verbosity to one of the Verbose.levels strings'

        if self._commandLineVerbose is not None:
            level = self._commandLineVerbose
        if level not in self.levels:
            raise ValueError('Illegal verbose string "%s".  Legal values are %s'%(level, self.levels))
        self.level = level

    def report(self, s, level='helpful'):
        """
        print message s to self.fileo if self.level>=level.  Return
        value indicates whether a message was issue.
        """
        if self.ge(level):
            print >>self.fileo, s
            return True
        return False

    def report_error(self, s):
        """
        print message s to self.fileo if self.level>=level.  Return
        value indicates whether a message was issued
        """
        if self.ge('error'):
            print >>self.erro, s
            return True
        return False


    def wrap(self, fmt, func, level='helpful', always=True):
        """
        return a callable function that wraps func and reports it
        output through the verbose handler if current verbosity level
        is higher than level

        if always is True, the report will occur on every function
        call; otherwise only on the first time the function is called
        """
        assert(callable, func)
        def wrapper(*args, **kwargs):
            ret = func(*args, **kwargs)

            if (always or not wrapper._spoke):
                spoke = self.report(fmt%ret, level)
                if not wrapper._spoke: wrapper._spoke = spoke
            return ret
        wrapper._spoke = False
        wrapper.__doc__ = func.__doc__
        return wrapper

    def ge(self, level):
        'return true if self.level is >= level'
        return self.vald[self.level]>=self.vald[level]


