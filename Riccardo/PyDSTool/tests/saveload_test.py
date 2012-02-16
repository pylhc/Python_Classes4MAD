#! /usr/local/bin/python
# Test pickling for saving and loading various PyDSTool objects

from PyDSTool import *
import os

# --------------------------------------------------------------------

# array
a=array([1,Inf])
b=[Inf,0]

saveObjects([a,b], 'temp_objects.pkl', True)
loadedObjs = loadObjects('temp_objects.pkl')
assert a[0]==loadedObjs[0][0]
assert a[1]==loadedObjs[0][1]
assert b[0]==loadedObjs[1][0]
print "Loading of stored arrays successful\n"

# Interval
m=Interval('test1', 'float', (-Inf,1))
s=Interval('a_singleton', 'float', 0.4)
print "m.get() => ", m.get()
saveObjects([m,s], 'temp_objects.pkl', True)
print "Saved Intervals", [obj.name for obj in [m,s]]
objs_ivals = loadObjects('temp_objects.pkl')
print "Loaded Intervals", [obj.name for obj in objs_ivals]
print " ... m.type => ", m.type
print "objs_ivals[0].get() => ", objs_ivals[0].get()

# Try loading partial list from a larger file
objs_part = loadObjects('temp_objects.pkl', ['a_singleton'])
print "Partial loading of stored file also successful\n"

# Point
x = Point(coorddict = {'x0': [1.123456789], 'x1': [-0.4], 'x2': [4000]},
       coordtype = Float64)

v = Pointset(coorddict = {'x0': 0.2, 'x1': -1.2},
            indepvardict = {'t': 0.01},
            coordtype = Float,
            indepvartype = Float)

saveObjects([x,v], 'temp_objects.pkl', True)
print "Saved Point and Pointset denoted x, v"
objs_pts = loadObjects('temp_objects.pkl')
print "Loaded Point and Pointset"
print "x => ", x
print "objs_pts[0] == x ?  => ", objs_pts[0] == x

print "\n"

# Simple Variable
var1 = Variable(Pointset(coordarray = array(range(10), Float)*0.1,
                         indepvararray = array(range(10), Float)*0.5
                        ), name='v1')
saveObjects(var1, 'temp_objects.pkl', True)
print "Saved variable object denoted var_int"
obj_var = loadObjects('temp_objects.pkl')[0]
print "Loaded variable object"
print "obj_var(1.5) => ", obj_var(1.5)

print "\n"

# Trajectory
var2 = Variable(Pointset(coordarray = array(range(10), Float)*0.25+1.0,
                         indepvararray = array(range(10), Float)*0.5
                        ), name='v2')
traj = Trajectory('traj1', [var1,var2])
saveObjects(traj, 'temp_objects.pkl', True)
print "Saved trajectory traj1"
traj_loaded = loadObjects('temp_objects.pkl')[0]
print "Loaded trajectory"
print "traj_loaded(2.0) =>\n", traj_loaded(2.0)

print "\n"

# Interpolated table Generator
xnames = ['x1', 'x2']
timeData = array([0.1, 1.1, 2.1])
x1data = array([10.2, -1.4, 4.1])
x2data = array([0.1, 0.01, 0.4])
xData = makeDataDict(xnames, [x1data, x2data])
itableArgs = {}
itableArgs['tdata'] = timeData
itableArgs['ics'] = xData
itableArgs['name'] = 'interp'
interptable = InterpolateTable(itableArgs)
itabletraj = interptable.compute('itable')
print "itabletraj(0.6) =\n", itabletraj(0.6)
saveObjects(itabletraj, 'temp_objects.pkl', True)
print "Saved interpolated table object"
obj_itab = loadObjects('temp_objects.pkl')
print "Loaded interpolated table object"
print "obj_itab[0](0.6) =\n", itabletraj(0.6)

print "\n"

# Vode object with event and external input trajectory (defined earlier)
fvarspecs = {"w": "k*w + a*itable + sin(t) + myauxfn1(t)*myauxfn2(w)",
               'aux_wdouble': 'w*2 + globalindepvar(t)',
               'aux_other': 'myauxfn1(2*t) + initcond(w)'}
fnspecs = {'myauxfn1': (['t'], '2.5*cos(3*t)'),
             'myauxfn2': (['w'], 'w/2')}
ev_args = {'name': 'threshold',
           'eventtol': 1e-4,
           'eventdelay': 1e-5,
           'starttime': 0,
           'term': True,
           }
thresh_ev = Events.makePythonStateZeroCrossEvent('w', 20, 1, ev_args)
DSargs = {'tdomain': [0.1,2.1],
          'tdata': [0.11, 2.1],
          'ics': {'w': 3.0},
          'pars': {'k':2, 'a':-0.5},
          'inputs': {'itable' : interptable.variables['x1']},
          'auxvars': ['aux_wdouble', 'aux_other'],
          'algparams': {'init_step':0.01, 'strict':False},
          'events': thresh_ev,
          'checklevel': 2,
          'name': 'ODEtest',
          'fnspecs': fnspecs,
          'varspecs': fvarspecs
          }
testODE = Vode_ODEsystem(DSargs)
odetraj = testODE.compute('testode')
print "Original events: ", testODE.showWarnings()
print "odetraj(0.6) =\n", odetraj(0.6)
saveObjects([odetraj, testODE], 'temp_objects.pkl', True)
print "Saved ODE with external input trajectory, event object, and its generator"
objs_ode = loadObjects('temp_objects.pkl')
print "Loaded ODE objects"
print "objs_ode[0](0.6) =\n", objs_ode[0](0.6)
objs_ode[1].clearWarnings()
print "objs_ode[1] warnings => ", objs_ode[1].showWarnings()
odetraj2 = objs_ode[1].compute('testode2')
print "odetraj2(0.6) =\n", odetraj2(0.6)
print "Events intact: ", objs_ode[1].showWarnings()

print "\n"

# ExplicitFnGen
args = {'tdomain': [-50, 50],
        'pars': {'speed': 1},
        'xdomain': {'s': [-1., 1.]},
        'name': 'sine',
        'globalt0': 0.4,
        'pdomain': {'speed': [0, 200]},
        'varspecs': {'s': "sin(globalindepvar(t)*speed)"}}
sin_gen = ExplicitFnGen(args)
sintraj1 = sin_gen.compute('sine1')
sin_gen.set(pars={'speed': 2})
sintraj2 = sin_gen.compute('sine2')
print "sintraj1(0.55) => ", sintraj1(0.55)
print "sintraj2(0.55) => ", sintraj2(0.55)
saveObjects([sin_gen,sintraj1,sintraj2], 'temp_objects.pkl', True)
print "Saved explicit function generator and its trajectory"
objs_sin = loadObjects('temp_objects.pkl')
print "sintraj1(0.55) == objs_sin[1](0.55) => ", objs_sin[1](0.55)
print "sintraj2(0.55) == objs_sin[2](0.55) => ", objs_sin[2](0.55)
print "Loaded objects successfully"

print "\n"

# ImplicitFnGen
fvarspecs = {"y": "t*t+y*y-r*r",
                "x": "t"
                }
argsi = {}
argsi['varspecs'] = fvarspecs
argsi['algparams'] = {'solvemethod': 'newton',
                     'atol': 1e-4}
argsi['xdomain'] = {'y': [-2,2]}
argsi['ics'] = {'y': 0.75}
argsi['tdomain'] = [-2,0]
argsi['pars'] = {'r':2}
argsi['vars'] = ['y']
argsi['checklevel'] = 2
argsi['name'] = 'imptest'

testimp = ImplicitFnGen(argsi)
traj1 = testimp.compute('traj1')
print "traj1(-0.4) =>\n", traj1(-0.4)
saveObjects([testimp, traj1], 'temp_imp.pkl', True)
print "Saved implicit function generator and its trajectory"
objs_imp = loadObjects('temp_imp.pkl')
print "Loaded objects successfully"
print "objs_imp[0]._xdomain => ", objs_imp[0]._xdomain
print "traj1(-0.4) == objs_imp[1](-0.4) =>\n", objs_imp[1](-0.4)

print "\n"

# Model
print "Tests for model in HH_model.py and HH_loaded.py"

print "\nAll tests completed successfully"
