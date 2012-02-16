# Lookup table
from __future__ import division

from allimports import *
from baseclasses import discGen, theGenSpecHelper
from PyDSTool.utils import *
from PyDSTool.common import *

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, transpose, shape
import math, random
from copy import copy, deepcopy

# -----------------------------------------------------------------------------

class LookupTable(discGen):
    """Lookup table trajectory with no interpolation.

    Independent and dependent variables may be integers or floats."""

    def __init__(self, kw):
        try:
            self._tdata = kw['tdata']
            self._xdatadict = {}
            for k, v in dict(kw['ics']).iteritems():
                self._xdatadict[str(k)] = v
            self.foundKeys = 2
            # check for other, invalid keys (but currently just ignored)
        except KeyError:
            raise PyDSTool_KeyError, 'Invalid keyword passed'
        self._tdomain = self._tdata
        discGen.__init__(self, kw)
        self.funcspec = {}
        if 'vars' in kw:
            raise PyDSTool_KeyError, 'vars option invalid for lookup table class'
        if 'auxvars' in kw:
            raise PyDSTool_KeyError, 'auxvars option invalid for lookup table class'
        if 'ttype' in kw:
            indepvartype = kw['ttype']
            assert indepvartype in _pynametype.keys(), 'Invalid ttype'
            self.foundKeys += 1
        else:
            indepvartype = 'float'
        if 'xtype' in kw:
            depvartype = kw['xtype']
            assert depvartype in _pynametype.keys(), 'Invalid xtype'
            self.foundKeys += 1
        else:
            depvartype = 'float'
        for x in self._xdatadict:
            self.funcspec[x] = Pointset({'coordarray': self._xdatadict[x],
                                         'coordtype': _pynametype[depvartype],
                                         'indepvararray': self._tdata,
                                         'indepvartype': _pynametype[indepvartype],
                                         'indepvarname': 't',
                                         'coordnames': x})
        self.needKeys.extend(['tdata', 'ics'])
        self.checkArgs(kw)
        self.indepvariable = Variable(listid, {'t_domain': self._tdomain},
                             {'t': self._tdata}, 't')
        self._register(self.indepvariable)
        for x in self._xdatadict:
            self.variables[x] = Variable(self.funcspec[x],
                                         {'t': copy(self.indepvariable.depdomain)},
                                         {x: self.funcspec[x].toarray()}, x)
        self._register(self.variables)
        self.dimension = len(self._xdatadict)
        self.validateSpec()
        self.defined = True


    def compute(self, trajname):
        if self.defined:
            self.validateSpec()
            self.clearWarnings()
            self.clearErrors()
        return Trajectory(trajname, [copy(v) for v in self.variables.values()],
                          self.globalt0, self.checklevel)


    computeTraj = compute


    def validateSpec(self):
        discGen.validateSpec(self)
        try:
            assert isoutputdiscrete(self.indepvariable)
            for v in self.variables.values():
                assert isinstance(v, Variable)
            assert not self.inputs
        except AssertionError:
            print 'Invalid system specification'
            raise


    def __del__(self):
        discGen.__del__(self)




# Register this Generator with the database

symbolMapDict = {}
# in future, provide appropriate mappings for libraries math,
# random, etc. (for now it's left to FuncSpec)
theGenSpecHelper.add(LookupTable, symbolMapDict, 'python')
