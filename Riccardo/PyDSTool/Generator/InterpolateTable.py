# Interpolated lookup table
from __future__ import division

from allimports import *
from baseclasses import ctsGen, theGenSpecHelper
from PyDSTool.utils import *
from PyDSTool.common import *

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, transpose, shape
import math, random
from copy import copy, deepcopy

# -----------------------------------------------------------------------------

class InterpolateTable(ctsGen):
    """Data lookup table with piecewise linear or piecewise constant interpolation."""

    def __init__(self, kw):
        try:
            self._tdata = kw['tdata']
            self._xdatadict = {}
            for k, v in dict(kw['ics']).iteritems():
                self._xdatadict[str(k)] = v
            self.foundKeys = 2
            # check for other, invalid keys (but currently just ignored)
        except KeyError:
            raise PyDSTool_KeyError, 'Keywords missing in argument'
        self._tdomain = extent(self._tdata)
        self._xdomain = {}
        for x in self._xdatadict:
            self._xdomain[x] = extent(self._xdatadict[x])
        ctsGen.__init__(self, kw)
        self.funcspec = {}
        if 'vars' in kw:
            raise PyDSTool_KeyError, 'vars option invalid for interpolated table class'
        if 'auxvars' in kw:
            raise PyDSTool_KeyError, 'auxvars option invalid for interpolated table class'
        if 'tdomain' in kw:
            raise PyDSTool_KeyError, 'tdomain option invalid for interpolated table class'
        if 'xdomain' in kw:
            raise PyDSTool_KeyError, 'xdomain option invalid for interpolated table class'
        if 'pdomain' in kw:
            raise PyDSTool_KeyError, 'pdomain option invalid for interpolated table class'
        if 'ttype' in kw:
            raise PyDSTool_KeyError, 'ttype option invalid for interpolated table class'
        if 'xtype' in kw:
            raise PyDSTool_KeyError, 'xtype option invalid for interpolated table class'
        if 'method' in kw:
            if kw['method']=='linear':
                interp=interp1d
            elif kw['method']=='constant':
                interp=interp0d
            else:
                raise ValueError, "Invalid interpolation method"
            self.foundKeys += 1
        else:
            # default to piecewise linear interpolation
            interp=interp1d
        for x in self._xdatadict:
            self.funcspec[x] = Pointset({'coordarray': self._xdatadict[x],
                                         'coordtype': Float,
                                         'indepvararray': self._tdata,
                                         'indepvartype': Float,
                                         'indepvarname': 't',
                                         'coordnames': x})
        self.needKeys.extend(['tdata', 'ics'])
        self.optionalKeys.append('method')
        self.checkArgs(kw)
        self.indepvariable = Variable(listid, Interval('t_domain', 'float',
                                              self._tdomain, self._abseps),
                             Interval('t', 'float', extent(self._tdata),
                                      self._abseps), 't')
        self._register(self.indepvariable)
        for x in self._xdatadict:
            self.variables[x] = Variable(interp(copy(self._tdata),
                                          self.funcspec[x].toarray()), 't',
                                  Interval(x, 'float', self._xdomain[x]), x)
        self._register(self.variables)
        self.dimension = len(self._xdatadict)
        self.validateSpec()
        self.defined = True


    def compute(self, trajname):
        return Trajectory(trajname, [copy(v) for v in self.variables.values()],
                          self.globalt0, self.checklevel)

    computeTraj = compute


    def validateSpec(self):
        ctsGen.validateSpec(self)
        try:
            assert isoutputcts(self.indepvariable)
            for v in self.variables.values():
                assert isinstance(v, Variable)
            assert not self.inputs
        except AssertionError:
            print 'Invalid system specification'
            raise


    def __del__(self):
        ctsGen.__del__(self)




# Register this Generator with the database

symbolMapDict = {}
# in future, provide appropriate mappings for libraries math,
# random, etc. (for now it's left to FuncSpec)
theGenSpecHelper.add(InterpolateTable, symbolMapDict, 'python')
