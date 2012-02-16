#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      1 Mar 2005
#$Revision: 1.0 $
"""
    Tests for the Generator class.  #5.
    This script tests:
    (1) the Explicit Function generator,
        its use of "global time", and its deletion.
    (2) Also test Implicit Function generator's deletion

    Robert Clewley, June 2005.
"""


from PyDSTool import *

ev_args = {'name': 'threshold',
           'eventtol': 1e-4,
           'eventdelay': 1e-5,
           'starttime': 0,
           'active': True,
           'term': True,
           'precise': True}
thresh_ev = Events.makePythonStateZeroCrossEvent('t', 20, 1, ev_args)

DSargs = {'tdomain': [-50, 50],
        'pars': {'speed': 1},
        'xdomain': {'s': [-1., 1.]},
        'name': 'sine',
        'globalt0': 0.4,
        'pdomain': {'speed': [0, 200]},
        'varspecs': {'s': "sin(globalindepvar(t)*speed)"},
        'events': thresh_ev}
sin_gen = ExplicitFnGen(DSargs)
sintraj1 = sin_gen.compute('sine1')
print "sintraj1.globalt0 => ", sintraj1.globalt0
print "sintraj1(0.) = ", sintraj1(0.)
print "'regular' sin(0.4) = ", sin(0.4)

sin_gen.set(pars={'speed': 2})
sintraj2 = sin_gen.compute('sine2')

print "sintraj2 independent variable domain truncated at terminal event"
print "(in the DS's local time)..."
print "sin_gen.getEventTimes() =>", sin_gen.getEventTimes()
print "sintraj2.indepdomain.get() =>", sintraj2.indepdomain.get()
print "\n"

# Simple implicit function
fvarspecs = {"y": "t*t+y*y-r*r",
             "x": "t"}
argsi = args()
argsi.varspecs = fvarspecs
argsi.algparams = {'solvemethod': 'newton',
                     'atol': 1e-4}
argsi.xdomain = {'y': [-2,2]}
argsi.ics = {'y': 0.75}
argsi.tdomain = [-2,0]
argsi.pars = {'r':2}
argsi.vars = ['y']
argsi.checklevel = 2
argsi.name = 'imptest'

testimp = ImplicitFnGen(argsi)
imptraj1 = testimp.compute('traj1')

print "Having created three trajectories containing Variable objects,"
print "here's what the Variable class has in its namespace:"
print "\ndir(Variable) =\n", dir(Variable)
dv=dir(Variable)
del sintraj1
print "\nOne trajectory has been deleted."
print "dir(Variable) =\n", dir(Variable)
assert dv != dir(Variable)
