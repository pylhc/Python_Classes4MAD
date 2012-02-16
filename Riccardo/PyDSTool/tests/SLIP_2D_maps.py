#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      10 Jul 2005
#$Revision: 1.0 $
"""
    2D SLIP in map form.

    Robert Clewley, July 2005.
"""

from PyDSTool import *
from SLIP_plot import SLIP_plot
from SLIP_2D import makeSLIP2D_Vode as makeSLIP2D
from copy import copy
from scipy.optimize import minpack, optimize
from scipy.linalg import eigvals, det, norm


# ---- Define pars and i.c.'s

k = 10.
g = 0.345
beta = 2*math.pi/5  # initial value only
pars = {'k': k, 'g': g, 'beta': beta}

z_ic = sin(beta)
y_ic = cos(beta)
info(pars, "Parameter values")
icdict = {'y': y_ic, 'z': z_ic, 'ydot': 0.8, 'zdot': .3, 'incontact': 0}
info(icdict, "Initial conditions")


# ---- Gait maps

SLIP_map = makeSLIP2D(pars, stop_at_TD=False, stop_at_LO=True)

def pdc_map(ic):
    SLIP_map.compute(trajname='pdc', force=1,
                         ics=ic, tdata=[0,12])
    # extract single point from pointset using [0] reference
    res = SLIP_map.getTrajEvents('pdc', 'liftoff')[0]
    res['incontact'] = 0 # make eligible for new i.c.
    return res


# Find fixed point of gait map parameterized by beta,
# for periodic gait
def fp_beta_fsolve(b, x0):
    SLIP_map.set(pars={'beta': b[0]})
    x0['z'] = math.sin(b[0])
    x0['y'] = math.cos(b[0])
    x1 = pdc_map(x0)
    return x1['zdot']-x0['zdot']  # Point


def fp_beta_newton(b, x0):
    SLIP_map.set(pars={'beta': b})
    x0['z'] = math.sin(b)
    x0['y'] = math.cos(b)
    x1 = pdc_map(x0)
    return x1['zdot']-x0['zdot']  # Point


icdict_pert = copy(icdict)
icpt = Point({'coorddict': icdict_pert})
print "Finding fixed point of periodic map as a function of beta"

# Equally efficient alternatives to shoot for f.p. solution
beta_pdc = minpack.fsolve(fp_beta_fsolve, beta, args=(icpt,), xtol=1e-4)
##beta_pdc = minpack.newton(fp_beta_newton, beta, args=(icpt,), tol=1e-4)
##beta_pdc = 1.21482619378
print "beta_pdc = ", beta_pdc


# update i.c. for new beta
icdict_pert['z'] = math.sin(beta_pdc)
icdict_pert['y'] = math.cos(beta_pdc)
icdict_maps = copy(icdict_pert)
icpt = Point({'coorddict': icdict_maps})
print "\nCalculating approximation to periodic maps"
SLIP_map.set(pars={'beta': beta_pdc})
states = [icpt]
for i in range(5):
    states.append(pdc_map(states[-1]))

pdcgaittraj = pointsToPointset(states, 'n', range(len(states)))


print "\nCalculating approximation to periodic trajectory"
pars['beta'] = beta_pdc
SLIP = makeSLIP2D(pars, dt=0.005)
SLIP.compute(trajname='pdc',
             tdata=[0, 12],  # t=12 is not reached if event occurs
             ics=icdict_pert,
             verboselevel=0)

print 'Plotting periodic trajectory (event dots from application of map)'
SLIP_plot(SLIP, 'pdc', 'plane')
# Plot successive (y,z) events points overlayed
nmax = pdcgaittraj.indepvararray[-1]
n=0
while n < nmax:
    LO = pdcgaittraj(n)
    plot(LO['y'],LO['z'],'go')
    if n + 1 <= nmax:
        TD = pdcgaittraj(n+1)
        plot(TD['y'],TD['z'],'ro')
        n += 2
        
show()

evsTD = SLIP.getTrajEventTimes('pdc')['touchdown']
T = evsTD[1]-evsTD[0]

def pdTmap(x):
    SLIP.compute(trajname='one_period',
                  tdata=[0, T],
                  ics=x,
                  force=True)
    return SLIP('one_period', T, ['y','z','ydot','zdot'])



print "Jacobian DP as a function of z, ydot, zdot (y contributes an e'val = 1)"
DP = diff(pdc_map, icpt, axes=['z','ydot','zdot'],
           vars=['z','ydot','zdot'], dx=2e-4)
print "DP(fp)   =\n",DP
print "|DP(fp)| =", det(DP)
evals = eigvals(DP)
print "\nDP's eigenvalues are ",evals
print "Max eigenvalue magnitude = ", max(abs(evals))

# verify beta = delta using the geometry
tTD = SLIP.getTrajEventTimes('pdc')['touchdown'][0]  # time at first touchdown
sTD = SLIP('pdc',tTD)  # system state at t=TD
delta=math.atan(-sTD('ydot')/sTD('zdot'))  # angle of velocity vector from horizontal
theta_qv=beta_pdc-delta
print "Geometry shows that beta = delta, as theta_qv = beta-delta ..."
print "theta_qv =", theta_qv

# Check against Eq. (14)
vmag_TD = math.sqrt(sum(SLIP('pdc',tTD,['ydot','zdot'])**2)) # vel magnitude
tLO = SLIP.getTrajEventTimes('pdc')['liftoff'][1]
vmag_LO = math.sqrt(sum(SLIP('pdc',tLO,['ydot','zdot'])**2))
vmag_LOtheory = math.sqrt(vmag_TD**2+2*g*(sin(beta_pdc)-SLIP('pdc',tLO,'z')))

print "Using Eq. (14), difference in velocity magnitudes v_LO given by"
print "simulation vs. conservation of energy =", vmag_LOtheory - vmag_LO
show()
