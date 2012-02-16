# Differential-delay system (incomplete)

from allimports import *
from baseclasses import ctsGen
from PyDSTool.utils import *
from PyDSTool.common import *

# Other imports
from scipy import Inf, NaN, isfinite, sometrue, alltrue
from numarray import array, arange, zeros, Float, Int, Complex, \
     Float64, Int32, Complex64, transpose, shape
import math, random
from copy import copy, deepcopy
try:
    # use pscyo JIT byte-compiler optimization, if available
    import psyco
    HAVE_PSYCO = True
except ImportError:
    HAVE_PSYCO = False

# -----------------------------------------------------------------------------

class DDEsystem(ctsGen):
    """Delay-differential equations.

    (incomplete)"""

    # spec of initial condition _interval_ will be interesting!
    # Use a Variable trajectory over that interval?
    def __init__(self, kw):
        ctsGen.__init__(self, kw)
        raise NotImplementedError


    def validateSpec(self):
        ctsGen.validateSpec(self)


    def __del__(self):
        ctsGen.__del__(self)
