#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      26 Aug 2006
#$Revision: 1.0 $
"""
    Embedded map generator.

    Robert Clewley, August 2006.
"""

from PyDSTool import *
from copy import copy, deepcopy
from scipy.optimize import minpack, optimize
from scipy.linalg import eigvals, det, norm

vspec = {'x':'m*(y+x-x*x*x/3)',
         'y':'n*y-x'}
pars = {'m':10, 'n': 0.1}

ev = Events.makeZeroCrossEvent("y", -1,
                            {'eventtol': 1e-5,
                             'eventdelay': 0.05,
                             'starttime': 0,
                             'term': True,
                             'name': 'y_ev'}, ['y'], [])

mapsys = Generator.Vode_ODEsystem(args(varspecs=vspec,
                                    pars=pars,
                                    algparams={'init_step':0.05, 'stiff': True},
                                    tdata=[0,60],
                                    ics={'x':0, 'y':0},
                                    name='embeddedsys',
                                    checklevel=2,
                                    events=ev))

sys=deepcopy(mapsys)

def getTraj(ics=None,pars=None,t1=50,termFlag=True,trajname='test',pplane=True):
    if ics is not None:
        sys.set(ics=ics)
    if pars is not None:
        sys.set(pars=pars)
    sys.set(tdata=[0,t1])
    sys.eventstruct.setTermFlag('y_ev', termFlag)
    traj = sys.compute(trajname)
    pts = traj.sample()
    if pplane:
        plot(pts['x'],pts['y'])
    else:
        plot(pts['t'],pts['x'])
        plot(pts['t'],pts['y'])
    return traj,pts


##traj,pts=getTraj({'x':-2,'y':0},termFlag=False)
##1/0


map_ics = {'x': 0.2}

map_args = args(name='pdc_map')
map_args.varspecs = {'x': 'res["x"]'}
map_args.vfcodeinsert_start = """
    embeddedsys.set(ics={'x': x, 'y': 0}, tdata=[0,20])
    traj=embeddedsys.compute('pdc')
##    print 't =', traj.indepdomain.get()
    ev = embeddedsys.getEvents()['y_ev']
##    print ev['t'][0]
    res = {'x':ev[0]['x']}"""
##    if ev is None:
##        raise RuntimeError('No return to cross-section!')
map_args.ignorespecial = ['res']
map_args.system = mapsys
map_args.ttype = 'int'
map_args.ics = map_ics
map_args.pars = pars
map_args.tdomain = [0,1]
map = MapSystem(map_args)

map.showSpec()

test_traj = map.compute('test', map_ics)
print test_traj.sample()

def return_map(ic):
    return map.compute('pdc_shoot', ic)(1)


# Find fixed point of map
def fp_residual(ic):
    x0 = ic[0]
##    print "\nfp_residual: x0 = ", x0
    v = return_map({'x': x0})
    print v
    return v-ic[0]


print "Finding fixed point of return map as function of initial condition"

# Equally efficient alternatives to shoot for f.p. solution
x0_guess = map_ics['x']
sol_pdc = minpack.fsolve(fp_residual, array([x0_guess]), xtol=1e-6)
print "sol_pdc = ", sol_pdc


traj,pts=getTraj({'y':0,'x':sol_pdc}, t1=10,
        termFlag=False)

map.set(ics={'x':sol_pdc})

print "Initializing PyCont to follow y=0 crossing as parameter m varied:"

pc=ContClass(map)
pcargs=args(name='FPu', type='FP-C',
  freepars = ['m'],
  StepSize = 1e-1,
  MaxNumPoints = 30, # 40 causes IndexError in PyCont/BifPoint.py, line 149
  MaxStepSize = 5e-1,
  LocBifPoints = 'all',
  verbosity = 2
 )
pc.newCurve(pcargs)

print "Running PyCont forward"

pc['FPu'].forward()

print repr(pc['FPu'].sol)
