#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      19 September 2006
#$Revision: 2.0.5 $
"""Interval class

  Robert Clewley, June 2005

Interval objects have attributes:
   name: str  (variable label, which DOES NOT have to be unique)
   typestr: 'int' or 'float' (immutable) - doesn't handle 'complex' yet
   type: type
   _loval (low endpoint val): numeric
   _hival (hi endpoint val): numeric
   _abseps: numeric
   issingleton: boolean
   _intervalstr: str
"""


# ----------------------------------------------------------------------------
# VERSION HISTORY
#
# Version 2.0.5, September 2006
#   Modified set method to allow singleton intervals to be specified using
#    [x, x] format as well as plain x, where type(x) is int or float.
#
# Version 2.0.4, September 2005
#   Added info and intersection methods to Interval.
#
# Version 2.0.3, July 2005
#   Added __copy__ method that uses pickle, to avoid invalid literal INF
#     values for floats, when copying semi- or bi-infinite intervals.
#
# Version 2.0.2, June 2005
#   Added pickling methods __getstate__ and __setstate__ to Interval.
#
# Version 2.0.1, June 2005
#   Changed all lowercase method names to camel case.
#
# ----------------------------------------------------------------------------

# Note: The integer intervals will later be used as the basis for
# supporting finitely-sampled real ranges.

## PyDSTool imports
from utils import *
from common import *
from errors import *
from utils import info as utils_info

## Other imports
from scipy import Inf, NaN, isfinite
import re, math
# override scipy (numeric) old form of numeric types
from numarray import array, Float, Int, Complex, Float64, Int32, Complex64
import copy

MIN_EXP = -15

# type identifiers
re_number = '([-+]?(\d+(\.\d*)?|\d*\.\d+)([eE][-+]?\d+)?)'

c=1  # contained const
n=0  # notcontained const
u=-1  # uncertain const


__all__ = ['notcontained', 'contained', 'uncertain', 'Interval',
           'isinterval', 'issingleton']


class IntervalMembership(int):
    """Numeric Interval membership type."""

    def __init__(self, *arg):
        val = arg[0]
        if val == -1:
            self.valstr = 'uncertain'
        if val == 0:
            self.valstr = 'notcontained'
        if val == 1:
            self.valstr = 'contained'
        if self.__int__() not in [-1, 0, 1]:
            raise ValueError, \
                  'invalid value for Numeric Interval membership type'

    def __repr__(self):
        return self.valstr

    __str__ = __repr__

    def __and__(self, v):
        """Interval `logical AND` for checking interval containment.
        e.g. if both endpoints are contained then the whole interval is."""
        sv = self.__int__()
        ov = v.__int__()
        if sv == n or ov == n:
            return notcontained
        elif sv == u or ov == u:
            return uncertain
        else: # sv == ov == c is only remaining case
            return contained

    def __rand__(self, v):
        return self.__and__(v)

    # For possible future use
    def __or__(self, v):
        """Interval `logical OR` for checking interval intersection.
        e.g. if at least one endpoint is contained then there is intersection."""
        sv = self.__int__()
        ov = v.__int__()
        if sv == c or ov == c:
            return contained
        elif sv == u and ov == u:
            return uncertain
        else:
            return notcontained

    def __ror__(self, v):
        return self.__or__(v)


global notcontained, contained, uncertain
notcontained = IntervalMembership(False)
contained = IntervalMembership(True)
uncertain = IntervalMembership(-1)



class Interval(object):
    """Numeric Interval class.

    Numeric Interval implementation for integer and float types.
    """

    # If the interval is not specified fully on initialisation then operations
    # on the object are limited to set()
    def __init__(self, name, typestr, intervalspec=None, abseps=1e-14):
        if not isinstance(name, str):
                raise PyDSTool_TypeError, 'Name must be a string'
        self.name = name
        try:
            self.type = _pynametype[typestr]
            self.typestr = typestr
        except:
            raise PyDSTool_TypeError, 'Incorrect type specified for interval'
        assert abseps > 0, "abseps argument must be positive"
        self._abseps = abseps
        self.defined = False  # default value
        if intervalspec is not None:
            self.set(intervalspec)
        else:
            # just declare the name and type
            self._loval = None
            self._hival = None


    def info(self, verboselevel=0):
        if verboselevel > 0:
            # info is defined in utils.py
            utils_info(self.__dict__, "Interval " + self.name,
                 recurseDepthLimit=1+verboselevel)
        else:
            print self.__repr__()


    # Return all the interval's relevant data in a tuple
    def __call__(self):
        if self.defined:
            return (self.name, self.type, self.typestr,\
                    (self._loval, self._hival))
        else:
            raise PyDSTool_ExistError, 'Interval undefined'


    def __eq__(self, other):
        assert isinstance(other, Interval)
        return self.defined == other.defined and \
               self.name == other.name and \
               self.type == other.type and \
               self._loval == other._loval and \
               self._hival == other._hival and \
               self._abseps == other._abseps


    def __ne__(self, other):
        return not self == other


    def __contains__(self, val):
        testresult = self.contains(val)
        if testresult == notcontained:
            return False
        elif testresult == contained:
            return True
        else:
            raise PyDSTool_UncertainValueError('cannot determine membership '
                                         '(uncertain)', val)


    def contains(self, val):
        """Report membership of val in the interval,
        returning type IntervalMembership."""

        try:
            if not self.defined:
                raise PyDSTool_ExistError, 'Interval undefined'
            if type(val) == _pytypefromtype[self.type]:
                compval = val
                if self.type == Int:
                    eps = 0
                else:
                    eps = self._abseps
            else:
                if type(val) == int:
                    compval = float(val)
                    eps = self._abseps
                elif type(val) == float:
                    # cannot ever compare a float in an integer interval,
                    # unless int(val) == val
                    if int(val) == val:
                        compval = int(val)
                        eps = 0
                    else:
                        raise PyDSTool_TypeError, 'Incorrect type of query value'
                elif not val.issingleton:
                    # catches non-intervals or non-singleton intervals
                    if not val.defined:
                        raise PyDSTool_ExistError, 'Input interval undefined'
                    if val.type != self.type and val.type == Float:
                        # meaningless to ask if float interval is contained in an
                        # integer interval!
                        raise PyDSTool_TypeError, 'Interval type mismatch'
                    if val.type == self.type == Int:
                        eps = 0
                    else:
                        eps = max(self._abseps, val._abseps)
                        try:
                            loexp = math.log(abs(self._loval), 10)
                        except OverflowError:
                            loexp = 0
                        try:
                            hiexp = math.log(abs(self._hival), 10)
                        except OverflowError:
                            hiexp = 0
                        minexpallowed = math.ceil(-MIN_EXP - max(loexp, hiexp))
                        if -math.log(eps,10) > minexpallowed:
                            eps = math.pow(10,-minexpallowed)
                    if isfinite(val._loval) or isfinite(self._loval):
                        tempIlo = val._loval >= (self._loval + eps)
                    else:
                        tempIlo = False
                    if isfinite(val._hival) or isfinite(self._hival):
                        tempIhi = val._hival <= (self._hival - eps)
                    else:
                        tempIhi = False
                    if tempIlo and tempIhi:
                        return contained
                    elif eps == 0:
                        return notcontained
                    else:
                        # having already tested for being contained, this is
                        # sufficient for uncertainty
                        if isfinite(val._loval) or isfinite(self._loval):
                            tempUlo = val._loval > (self._loval - eps)
                            tempElo = val._loval <= (self._loval - eps)
                        else:
                            tempUlo = val._loval == self._loval
                            tempElo = False
                        if isfinite(val._hival) or isfinite(self._hival):
                            tempUhi = val._hival < (self._hival + eps)
                            tempEhi = val._hival >= (self._hival + eps)
                        else:
                            tempUhi = val._hival == self._hival
                            tempEhi = False
                        if (tempUlo and not tempEhi) or (tempUhi and \
                                                         not tempElo):
                            return uncertain
                        else:
                            # else must be notcontained
                            return notcontained
                else:
                    # is a singleton interval type
                    # issingleton == True implies interval is defined
                    # Now go through same sequence of comparisons
                    if val.type == self.type:
                        compval = val.get()
                        if self.type == Int:
                            eps = 0
                        else:
                            eps = max(self._abseps, val._abseps)
                            try:
                                loexp = math.log(abs(self._loval), 10)
                            except OverflowError:
                                loexp = 0
                            try:
                                hiexp = math.log(abs(self._hival), 10)
                            except OverflowError:
                                hiexp = 0
                            minexpallowed = math.ceil(-MIN_EXP - max(loexp,
                                                                     hiexp))
                            if -math.log(eps,10) > minexpallowed:
                                eps = math.pow(10,-minexpallowed)
                    else:
                        if val.type == Int:
                            compval = val.get()
                            eps = self._abseps
                        elif val.type ==  Float:
                            # cannot ever compare a float in an integer interval
                            # unless bd values are equal to their int() versions
                            if int(val._loval) == val._loval and \
                               int(val._hival) == val._hival:
                                compval = (int(val.get(0)), int(val.get(1)))
                                eps = 0
                            else:
                                raise PyDSTool_TypeError, \
                                      'Invalid numeric type of query value'
                        else:  # unexpected (internal) error
                            raise PyDSTool_TypeError, \
                                  'Invalid numeric type of query value'
        except AttributeError:
            raise PyDSTool_TypeError, ('Expected a numeric type or a singleton '
                                 'interval. Got type '+str(type(val)))
        else:
            tempIlo = compval >= (self._loval + eps)
            tempIhi = compval <= (self._hival - eps)
            if tempIlo and tempIhi:
                return contained
            elif eps == 0: # only other possibility (no uncertainty)
                return notcontained
            else:
                # having already tested for being contained, this is
                # sufficient for uncertainty
                tempUlo = compval > (self._loval - eps)
                tempUhi = compval < (self._hival + eps)
                tempElo = compval <= (self._loval - eps)
                tempEhi = compval >= (self._hival + eps)
                if (tempUlo and not tempEhi) or (tempUhi and not tempElo):
                    return uncertain
                # else must be notcontained
                else:
                    return notcontained


    def intersect(self, other):
        if not isinstance(other, Interval):
            raise PyDSTool_TypeError("Can only intersect with other Interval "
                                     "types")
        result = None    # default, initial value if no intersection
        do_floats = self.type == Float or other.type == Float
        if self.type == Complex or other.type == Complex:
            raise TypeError, "Complex Intervals not supported"
        if do_floats:
            # mixed types or both floats results in float interval
            typestr = 'float'
        else:
            typestr = 'int'
        if self.contains(other):
            result = other
        elif other.contains(self):
            result = self
        else:
            # no total containment possible
            if other.contains(self._hival) is contained:
                # then also self.contains(other._loval)
                result = Interval('__result__', typestr, [other._loval,
                                                    self._hival])
            elif other.contains(self._loval) is contained:
                # then also self.contains(other._hival)
                result = Interval('__result__', typestr, [self._loval,
                                                    other._hival])
        return result


    def atEndPoint(self, val, bdcode):
        """val, bdcode -> Bool

        Determines whether val is at the endpoint specified by bdcode,
        to the precision of the interval's _abseps tolerance.
        bdcode can be one of 'lo', 'low', 0, 'hi', 'high', 1"""

        assert self.defined, 'Interval undefined'
        assert isinstance(val, float) or isinstance(val, int), \
               'Invalid value type'
        assert isfinite(val), "Can only test finite argument values"
        if bdcode in ['lo', 'low', 0]:
            if self.type == Float:
                return abs(val - self._loval) < self._abseps
            elif self.type == Int:
                return val == self._loval
            else:
                raise TypeError, "Unsupported value type"
        elif bdcode in ['hi', 'high', 1]:
            if self.type == Float:
                return abs(val - self._hival) < self._abseps
            elif self.type == Int:
                return val == self._hival
            else:
                raise TypeError, "Unsupported value type"
        else:
            raise ValueError, 'Invalid boundary spec code'


    # Uniformly sample the interval, returning a list of points
    # 'Strict' option forces dt to be used throughout the interval,
    # with a final step of < dt
    def uniformSample(self, dt, strict=False, avoidendpoints=False):
        """dt, strict=False, avoidendpoints=False -> point list.
        Default strict = False used for auto-selection of sample rate
          to fit interval (choice based on dt)."""

        assert self.defined
        intervalsize = self._hival - self._loval
        assert isfinite(intervalsize), "Interval must be finite"
        if dt > intervalsize:
            print "Interval size = %f, dt = %f"%(intervalsize, dt)
            raise ValueError, 'dt must be smaller than size of interval'
        if dt <= 0:
            raise ValueError, 'Must pass dt >= 0'
        if self.type == Float:
            if strict:
                # moved int() to samplist's xrange
                n = max(math.floor(intervalsize/dt),2)
                # special case: add a point at t=dt
                # otherwise just get endpoints returned!
            else: # choose automatically
                n = max(round(intervalsize/dt),2)
                dt = intervalsize/n
            samplelist = [self._loval + i*dt for i in xrange(int(n))]
            if avoidendpoints:
                samplelist.append(self._hival - 1.5*self._abseps)
                samplelist[0] = self._loval + 1.5*self._abseps
            else:
                samplelist.append(self._hival)
        elif self.type == Int:
            if not isinstance(dt, int):
                raise ValueError, ("dt must be an integer for integer "
                                        "intervals")
            if strict:
                print "Warning: 'strict' option is invalid for integer " + \
                      "interval types"
            if avoidendpoints:
                loval = self._loval+1
                hival = self._hival-1
                assert loval <= hival, ('There are no points to return with '
                                        'these options!')
            else:
                loval = self._loval
                hival = self._hival
            samplelist = range(loval, hival+1, dt)
            # extra +1 on hival because of python range() policy!
        else:
            raise TypeError, "Unsupported value type"
        return samplelist            


    # Define interval in a Interval object
    def set(self, arg):
        if type(arg) in [list, tuple, Array, NArray] and len(arg)==2:
            if arg[0] == arg[1]:
                self.set(arg[0])
            else:
                self.issingleton = False
                loval = arg[0]
                hival = arg[1]
                assert loval is not NaN and hival is not NaN, \
                       ("Cannot specify NaN as interval endpoint")
                if not loval < hival:
                    print "set() was passed loval = ", loval, \
                          " and hival = ", hival
                    raise PyDSTool_ValueError, ('Interval endpoints must be '
                                    'given in order of increasing size')
                self._intervalstr = '['+str(loval)+',' \
                                    +str(hival)+']'
                self._loval = loval
                self._hival = hival            
                self.defined = True
        elif type(arg) in [int, float]:
            assert isfinite(arg), \
                   "Singleton interval domain value must be finite"
            self.issingleton = True
            self._intervalstr = str(arg)
            self._loval = arg
            self._hival = arg
            self.defined = True
        else:
            print "Error in argument: ", arg, "of type", type(arg)
            raise PyDSTool_TypeError, \
                  'Interval spec must be a numeric or a length-2 sequence type'


    # Get the interval as a tuple or a number (for singletons),
    # or an endpoint if ix is not None
    def get(self, ix=None):
        if self.defined:
            if ix == 0:
                return self._loval
            elif ix == 1:
                return self._hival
            elif ix is None:
                if self.issingleton:
                    return self._loval
                else:
                    return [self._loval, self._hival]
            else:
                raise PyDSTool_TypeError, 'Invalid return form specified'
        else:
            raise PyDSTool_ExistError, 'Interval undefined'


    def _infostr(self, verbose=1):
        """Get info on a known interval definition."""

        if verbose > 0:
            infostr = "Interval "+self.name+"\n"
            if self.defined:
                infostr += '  ' + self.typestr+': '+\
                          self.name+' = '+self._intervalstr+ \
                          ' @ eps = '+str(self._abseps)
            else:
                infostr += '  ' + self.typestr+': '+self.name+' @ eps = '+ \
                          str(self._abseps) + " (not fully defined)"
        else:
            infostr = "Interval "+self.name
        return infostr


    def __repr__(self):
        return self._infostr(verbose=0)


    __str__ = __repr__


    def info(self, verboselevel=1):
        print self._infostr(verboselevel)


    def __copy__(self):
        pickledself = pickle.dumps(self)
        return pickle.loads(pickledself)


    def __getstate__(self):
        d = copy.copy(self.__dict__)
        # remove reference to Cfunc self.type
        d['type'] = None
        return d


    def __setstate__(self, state):
        self.__dict__.update(state)
        # reinstate Cfunc self.type
        self.type = _pynametype[self.typestr]
##        setattr(self, 'type', _pynametype[self.typestr])
        

#--------------------------------------------------------------------

# Exported utility functions

def isinterval(obj):
    """Determines whether the given obj is a Interval object."""
    return isinstance(obj, Interval)


# Check whether an interval is a singleton
def issingleton(ni):
    if ni.defined:
        return ni.issingleton
    else:
        raise PyDSTool_ExistError, 'Interval undefined'


#--------------------------------------------------------------------

## TESTING

if __name__ == '__main__':
    print float
    print Float
    print '---------- Test for Interval.py ----------'
    print "m=Interval('test1', 'float', (0,1))"
    m=Interval('test1', 'float', (0,1))
    print 'm =>   ', m
    print 'm() =>   ', m()
    print 'm.info() =>  ',
    m.info()
    print '0.1 in m =>   ', 0.1 in m
    print
    print "n=Interval('test2', 'float', (0,0.4))"
    n=Interval('test2', 'float', (0,0.4))
    try:
        print 'n in m =>  ', n in m, '\n'
    except PyDSTool_UncertainValueError, e:
        print 'UncertainValueError: ', e
    print "s=Interval('a_singleton', 'float', 0.4)"
    s=Interval('a_singleton', 'float', 0.4)
    print 's.get() =>  ', s.get()
    print "s in m =>  ", s in m, '  (works for singleton interval s)\n'
    r=Interval('another_singleton', 'float', 0.0)
    print "r=Interval('another_singleton', 'float', 0.0)"
    b=Interval('b', 'int', 0)
    print "b=Interval('b', 'int', 0)"
    try:
        print "b in m =>  ", b in m
    except PyDSTool_UncertainValueError, e:
        print 'UncertainValueError: ', e
    i=Interval('i', 'int', (0,1))
    print "i=Interval('i', 'int', (0,1))"
    print "b in i =>  ", b in i, " (true because we're comparing integer intervals)"
    
    print "\nUse the explicit `contains` method to avoid exceptions, and instead"
    print "   get a NumIntMembership type returned..."
    print "m.contains(i) =>  ", m.contains(i)
    print "m.contains(0.4) => ", m.contains(0.4)
    
    j = Interval('test3', 'float', (0,0.999999999))
    print "j = Interval('test3', 'float', (0,0.999999999))"
    print "p = m.contains(j)"
    p = m.contains(j)
    print "p is uncertain => ", p is uncertain
    
    print "\nBut don't try to compare NumIntMembership objects to booleans..."
    print "q = m.contains(0.9)"
    q = m.contains(0.9)
    print "q is True => ", q is True, " (false because q is not a boolean type)"
    print "... but can use in a statement such as 'if m.contains(0.9): ...etc.'"

    print "\nElementary `interval logic` can be performed when checking endpoints"
    print "   for interval containment."
    print "contained and notcontained => ", contained and notcontained
    print "contained and uncertain => ", contained and uncertain
    print "notcontained and notcontained => ", notcontained and notcontained

    print "\nm.uniformSample(0.09, strict=False, avoidendpoints=True) => ", \
          m.uniformSample(0.09, strict=False, avoidendpoints=True) 

    print "\nm.uniformSample(0.09, strict=False) => ", \
          m.uniformSample(0.09, strict=False) 

    print "i2=Interval('i2', 'int', (0,10))"
    i2=Interval('i2', 'int', (0,10))
    print "\ni2.uniformSample(2, strict=False, avoidendpoints=True) => ", \
          i2.uniformSample(2, strict=False, avoidendpoints=True)

    print "i3=Interval('i3', 'float', (0.,0.4))"
    i3=Interval('i3', 'float', (0.,0.4))
    print "\ni3.uniformSample(0.36, strict=False) => ", \
          i3.uniformSample(0.36, strict=False)

    print "\ni3.uniformSample(0.36, strict=False, avoidendpoints=True) => ", \
          i3.uniformSample(0.36, strict=False, avoidendpoints=True)

    print "\ni3.uniformSample(0.36, strict=True) => ", \
          i3.uniformSample(0.36, strict=True)
    


    print "\nInfinite intervals"
    print "inf1 = Interval('inf1', 'float', [0,Inf])"
    inf1 = Interval('inf1', 'float', [0,Inf])
    print "inf1.contains(inf1) => ", inf1.contains(inf1)
    print "inf2 = Interval('inf2', 'float', [-Inf,Inf])"
    inf2 = Interval('inf2', 'float', [-Inf,Inf])
    print "inf2.contains(inf2) => ", inf2.contains(inf2)
    print "inf2.contains(inf1) => ", inf2.contains(inf1)
    print "inf3 = Interval('inf3', 'float', [-Inf,0])"
    inf3 = Interval('inf3', 'float', [-Inf,0])
    print "inf3.contains(inf2) => ", inf3.contains(inf2)
    inf_int = Interval('inf3', 'int', [-Inf,0])
    print "inf_int = Interval('inf3', 'int', [-Inf,0])"
    print "inf_int.contains(inf3) => "
    try:
        inf_int.contains(inf3)
    except PyDSTool_TypeError, e:
        print " ",e
    i4 = Interval('i4', 'int', [-5,5])
    print "i4 = Interval('i4', 'int', [-5,5])"
    print "inf_int.intersect(i4) => ", inf_int.intersect(i4).get()
    print "Intersection should fail on mixed-type intervals ..."
    try:
        result1 = inf3.intersect(i4).get()
    except:
        result1 = " >> FAILURE"
    try:
        result2 = j.intersect(i2).get()
    except:
        result2 = " >> FAILURE"
    print "inf3.intersect(i4) => ", result1
    print "j.intersect(i2) => ", result2
    i5 = Interval('i5', 'int', [4,4])
    assert i5.issingleton
    print "Tests passed"
