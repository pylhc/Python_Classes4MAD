#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      10 Jul 2005
#$Revision: 1.0 $
"""
    Find SLIP periodic orbit.

    Robert Clewley, July 2005.
"""

from PyDSTool import *
from SLIP_plot import SLIP_plot
from SLIP_2D import *

from scipy.optimize import minpack, optimize
from time import clock
from copy import copy

# ---- Find periodic gait

# ---- Define pars and i.c.'s

k = 10.
g = 0.345
beta = 2*math.pi/5
##beta = 1.25305395342

pars = {'k': k, 'g': g, 'beta': beta}
info(pars, "Parameter values")

z_ic = sin(beta)
y_ic = math.sqrt(1-z_ic*z_ic)
ics = {'y': y_ic, 'z': z_ic, 'ydot': .6, 'zdot': .3}
info(ics, "Initial conditions")

# Choose makeSLIP2D_Dopri or makeSLIP2D_Vode as alternatives...
SLIP = makeSLIP2D_Radau(pars)

# ---- Find periodic gait -----------------------------------------

# residual fn for searching in beta
def residual_fn_beta(x):
    ics['incontact'] = 0
    print "Trying beta =", x[0]
    try:
        SLIP.compute(pars={'beta': x[0]}, force=True,
                         trajname='par_est', tdata=[0, 12],
                         ics=ics)
    except PyDSTool_ExistError:
        # beta chosen such that no eligible generators found
        print " ... arbitrarily setting cost to be 1000"
        return 1000
    # time of touchdown and liftoff
    evs = SLIP.getTrajEventTimes('par_est')
    numTDs = len(evs['touchdown'])
    numLOs = len(evs['liftoff'])
    if numTDs >= 2 and numLOs >= 3:
        TDev = evs['touchdown'][1]
        LOev = evs['liftoff'][1]
        if LOev < TDev:
            LOev = evs['liftoff'][2]
        z3 = SLIP('par_est',LOev,'z')
        print "z3=",z3
        delta = math.asin(z3)
        print "delta =",delta
        # Dpsi = pi-beta-delta
        # target Dpsi = pi-2*beta, i.e. beta=delta
        cost = (beta-delta)**2
        print "cost =", cost
    else:
        raise RuntimeError, "Not enough events found"
    return cost

# use optimizer with boundary constraints
##beta_opt = minpack.fsolve(residual_fn_beta, beta, xtol=Dpsi_tol)

# -----------------------------------------------------------------

def residual_fn_ydot(x):
    ics['ydot'] = x
    ics['incontact'] = 0
##    print "\n** Trying ydot =", x
    try:
        SLIP.compute(force=True,
                         trajname='par_est', tdata=[0, 20],
                         ics=ics)
    except PyDSTool_ExistError:
        # ydot chosen such that no eligible generators found
        print " ... arbitrarily setting cost to be 1000"
        return 1000
    # time of touchdown and liftoff
    evs = SLIP.getTrajEventTimes('par_est')
    numTDs = len(evs['touchdown'])
    numLOs = len(evs['liftoff'])
    if numTDs >= 2 and numLOs >= 3:
        TDev = evs['touchdown'][1]
        LOev = evs['liftoff'][1]
        if LOev < TDev:
            LOev = evs['liftoff'][2]
        z3 = SLIP('par_est',LOev,'z')
##        print "z3=",z3
        delta = math.asin(z3)
##        print "delta =",delta
        # Dpsi = pi-beta-delta
        # target Dpsi = pi-2*beta, i.e. beta=delta
        cost = (beta-delta)**2
        print "cost =", cost
    else:
        print evs
        print "Found %i TD events, %i LO events"%(numTDs, numLOs)
        raise RuntimeError, "Not enough events found"
    return cost

Dpsi_tol = 5e-3
print "Finding periodic gait by varying initial y velocity,"
print "to tolerance in Dpsi of", Dpsi_tol, "\n"

# use optimizer with boundary constraints
ydot_opt = optimize.fminbound(residual_fn_ydot, 0.5, 0.75, xtol=Dpsi_tol)
#ydot_opt = 0.666463229031
print "ydot for periodic gait = ", ydot_opt


# ---- Compute periodic trajectory

icdict_pdc = copy(ics)
icdict_pdc['ydot'] = ydot_opt
icdict_pdc['incontact'] = 0

print "Computing trajectory...\n"
start = clock()
SLIP.compute(trajname='pdc',
                 tdata=[0, 12],
                 ics=icdict_pdc,
                 verboselevel=0)   # optional
print '... finished in %.3f seconds.\n' % (clock()-start)

print 'Plotting periodic trajectory'

SLIP_plot(SLIP, 'pdc', 'plane')
show()
