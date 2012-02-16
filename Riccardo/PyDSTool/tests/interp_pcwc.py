#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      11 March 2006
#$Revision: 1.0 $
"""
    Tests for the Generator class InterpolateTable with
    piecewise-constant vs. piecewise-linear interpolation.

    Robert Clewley, March 2006.
"""

from PyDSTool import *


print '-------- Test: InterpolateTable comparing piecewise-constant and piecewise-'
print 'linear interpolation of sine function'

timeData = linspace(0, 10, 30)
sindata = sin(timeData)
xData = makeDataDict(['sinx'], [sindata])
pcwc_interp = InterpolateTable({'tdata': timeData,
                              'ics': xData,
                              'name': 'interp0d',
                              'method': 'constant',
                              'checklevel': 1,
                              'abseps': 1e-5
                              }).compute('interp')

pcwl_interp = InterpolateTable({'tdata': timeData,
                              'ics': xData,
                              'name': 'interp1d',
                              'method': 'linear',
                              'checklevel': 1,
                              'abseps': 1e-5
                              }).compute('interp')


x = linspace(0,10,300)
plot(x, pcwc_interp(x, ['sinx']))
plot(x, pcwl_interp(x, ['sinx']))
plot(x, sin(x))
