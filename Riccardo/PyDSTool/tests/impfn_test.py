#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      22 June 2005
#$Revision: 1.0 $
"""
    Tests for the Generator class.  #4.
    Test Implicit Function generator
    
    Robert Clewley, June 2005.
"""

from PyDSTool import *
try:
    from scipy import gplt
except ImportError:
    do_gnu = False
    print "Gnuplot not available - will not plot graph"
else:
    do_gnu = True
from copy import copy

print "------- Examples of implicit function as variable"
print "1D example: a half-circle using newton's method (secant method for estimating derivative)"
print "  Change sign of 'y' initial condition to solve for other half-circle"
fvarspecs = {"y": "t*t+y*y-r*r",
                "x": "t"
                }
##fnspecs = {'myauxfn': (['t'], '.5*cos(3*t)')}
dsargs=args()
##dsargs.fnspecs = fnspecs
dsargs.varspecs = fvarspecs
dsargs.algparams = {'solvemethod': 'newton',
                     'atol': 1e-4}
dsargs.xdomain = {'y': [-2,2]}
dsargs.ics = {'y': 0.75}
dsargs.tdomain = [0,2]
dsargs.pars = {'r':2}
dsargs.vars = ['y']
dsargs.checklevel = 2
dsargs.name = 'imptest'

testimp = ImplicitFnGen(dsargs)
print "Gen defined? => ", testimp.defined
print "traj1 = testimp.compute('traj1')"
traj1 = testimp.compute('traj1')
print "Gen defined now? => ", testimp.defined
print "traj1(1.5) => \n", traj1(1.5)
print "\ntraj1.dimension => ", traj1.dimension
print "traj1(0.8, ['y']) => ", traj1(0.8, ['y'])

print "\n\n2D example: a quarter-sphere with linear constraint y = 3z,"
print "using MINPACK's fsolve"
radius = 2.
print "  ... radius parameter = ", radius
fvarspecs2d = {"y": "t*t+y*y+z*z-r*r",
                  "z": "y-3.*z",
                  "x": "t"
                  }
args2d = {'varspecs': fvarspecs2d,
          'algparams': {'solvemethod': 'fsolve',
                        'atol': 1e-3},
          'xdomain': {'y': [-2,2], 'z': [-2,2]},
          'ics': {'y': 0.75, 'z': 0.9},
          'tdomain': [-2,2],
          'pars': {'r':radius},
          'vars': ['y', 'z'],
          'checklevel': 2,
          'name': 'imptest2d'
          }
testimp2d = ImplicitFnGen(args2d)

traj2 = testimp2d.compute('traj2')
p = traj2(1.5)
print "traj2(1.5) => \n", p
print "length of vector = ", sqrt(p('x')*p('x')+p('y')*p('y')+p('z')*p('z')), \
      " (should = radius parameter ", radius, ")"
print "for x = 1.5, y - 3z = ", p('y') - 3*p('z')
print "\ntraj2.dimension => ", traj2.dimension
print "isparameterized(traj2) => ", isparameterized(traj2)

print "\nTest bounds checking:"
try:
    traj2(3.)
    print "BOUNDS TEST FAILED"
except ValueError, e:
    print "traj2(3.) => Error: ", e
    print "BOUNDS TEST PASSED"

ps1 = traj2.sample(dt=0.02)

testimp2d.set(ics={'y': -0.75})
traj3 = testimp2d.compute('traj3')
ps2 = traj3.sample(dt=0.02)

if do_gnu:
    plotarray1 = transpose(ps1.toarray())
    plotarray2 = transpose(ps2.toarray())

    constraintline = array([[2.,-1,0],[2.,1,0]])
    gplt.plot3d((plotarray1,plotarray2,constraintline))
    gplt.zaxis([-2,2])
    gplt.yaxis([-2,2])
    gplt.view([80,85,1])
    #gplt.close()
