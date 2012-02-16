#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      11 Oct 2005
#$Revision: 1.0 $
"""
    Version of 2D spring-loaded inverted pendulum (SLIP)
    using the VODE/Dopri/Radau integrators.

    Robert Clewley, July 2005.
"""

# Adapted from Justin Seipel's Matlab code

from PyDSTool import *
from copy import copy

# --------------------------------------------------------------------------

def getStanceArgs(grav_in_stance, algparams, liftoff_ev, pars, abseps):
    if grav_in_stance:
        grav_str = "-g"
    else:
        grav_str = ""
    return {'pars': pars,
            'fnspecs': {'zeta': (['y', 'z'], 'sqrt(y*y+z*z)'),
##                          initial conditions cannot be accessed in event
##                          'yrel': (['v'], 'v-initcond(y)+cos(beta)')},
                          'yrel': (['y','yic'], 'y-yic-cos(beta)')},
            'varspecs': {'y':    "ydot",
                            'ydot': "k*(1/zeta(yrel(y,initcond(y)),z)-1)*yrel(y,initcond(y))",
                            'z':    "zdot",
                            'zdot': "k*(1/zeta(yrel(y,initcond(y)),z)-1)*z"+grav_str,
                            'zetaval': "zeta(yrel(y,initcond(y)),z)",
                            'incontact': "0"
                            },
            'auxvars': ['zetaval'],
            'xdomain': {'y': [0, Inf], 'z': [0,Inf], 'zetaval': [0,1],
                        'incontact': 1},
            'ics': {'zetaval': 1., 'incontact': 1},
            'algparams': algparams,
            'events': liftoff_ev,
            'abseps': abseps,
            'name': 'stance'}


def getFlightArgs(algparams, touchdown_ev, pars, abseps):
    return {'pars': pars,
           'fnspecs': {'zeta': (['y', 'z'], 'sqrt(y*y+z*z)')},
           'varspecs': {'y':    "initcond(ydot)",
                           'ydot': "0",
                           'z':    "zdot",
                           'zdot': "-g",
                           'zetaval': "0.",
                           'incontact': "0"},
           'xdomain': {'y': [0, Inf], 'z': [0,Inf], 'zetaval': [0,1],
                       'incontact': 0},
           'ics': {'zetaval': 1., 'incontact': 0},
           'algparams': algparams,
           'events': touchdown_ev,
           'abseps': abseps,
           'name': 'flight'}


def makeSLIP2D_Vode(pars, dt=0.1, abseps=1e-4, grav_in_stance=True,
                    stop_at_TD=False, stop_at_LO=False):
    allgen_names = ['stance', 'flight']
    stance_args, flight_args = makeDS_parts(pars, dt, abseps, grav_in_stance)
    stanceDS = Generator.Vode_ODEsystem(stance_args)
    flightDS = Generator.Vode_ODEsystem(flight_args)

    return makeSLIPModel(stanceDS, flightDS, stop_at_TD, stop_at_LO)


def makeSLIP2D_Dopri(pars, dt=0.1, abseps=1e-4, grav_in_stance=True,
                     stop_at_TD=False, stop_at_LO=False):
    allgen_names = ['stance', 'flight']
    stance_args, flight_args = makeDS_parts(pars, dt, abseps, grav_in_stance,
                                            'c')
    stanceDS = Generator.Dopri_ODEsystem(stance_args)
    flightDS = Generator.Dopri_ODEsystem(flight_args)

    return makeSLIPModel(stanceDS, flightDS, stop_at_TD, stop_at_LO)


def makeSLIP2D_Radau(pars, dt=0.1, abseps=1e-4, grav_in_stance=True,
                     stop_at_TD=False, stop_at_LO=False):
    stance_args, flight_args = makeDS_parts(pars, dt, abseps, grav_in_stance,
                                            'c')
    stanceDS = Generator.Radau_ODEsystem(stance_args)
    flightDS = Generator.Radau_ODEsystem(flight_args)

    return makeSLIPModel(stanceDS, flightDS, stop_at_TD, stop_at_LO)


def makeSLIPModel(stanceDS, flightDS, stop_at_TD, stop_at_LO):
    allgen_names = ['stance', 'flight']
    toflight_map=EvMapping({"xdict['incontact']": "0"})

    tostance_map=EvMapping({"xdict['incontact']": "1"})

    if stop_at_TD:
        flightDS_info = makeGenInfoEntry(flightDS, allgen_names,
                                    [('touchdown', 'terminate')])
    else:
        flightDS_info = makeGenInfoEntry(flightDS, allgen_names,
                                     [('touchdown', ('stance', tostance_map))])
    if stop_at_LO:
        stanceDS_info = makeGenInfoEntry(stanceDS, allgen_names,
                                    [('liftoff', 'terminate')])
    else:
        stanceDS_info = makeGenInfoEntry(stanceDS, allgen_names,
                                     [('liftoff', ('flight', toflight_map))])
    genInfoDict = makeGenInfo([stanceDS_info, flightDS_info])

    SLIP = Model.Model({'name': 'SLIP', 'genInfo': genInfoDict})
    SLIP.forceIntVars(['zetaval', 'incontact'])
    return SLIP


def makeDS_parts(pars, dt, abseps, grav_in_stance=True,
                 targetlang='python'):
    # In contrast to VODE version, variable zetaval is a full ODE variable
    # with a RHS = 0, and not an auxiliary variable that is only an
    # algebraic function of the ODE variables. This difference is due to
    # zetaval not being excluded from the (default) 'vars' keyword in
    # the Dopri/Radau Generator initialization arguments.
    algparams = {'init_step': dt, 'max_step': dt}
    liftoff_args = {'eventtol': abseps/10,
               'eventdelay': abseps*10,
               'eventinterval': abseps*10,
               'active': True,
               'term': True,
               'precise': True,
               'name': 'liftoff'}
    liftoff_ev = Events.makeZeroCrossEvent('zeta(yrel(y,initcond(y)),z)-1', 1,
                          liftoff_args, ['y','z'], ['beta'],
                          fnspecs={'zeta': (['y', 'z'], 'sqrt(y*y+z*z)'),
                                     'yrel': (['y','yic'], 'y-yic-cos(beta)')},
                          targetlang=targetlang)

    stance_args = getStanceArgs(grav_in_stance, algparams, liftoff_ev, pars,
                                abseps)

    touchdown_args = {'eventtol': abseps/10,
               'eventdelay': abseps*10,
               'eventinterval': abseps*10,
               'active': True,
               'term': True,
               'precise': True,
               'name': 'touchdown'}
    touchdown_ev = Events.makeZeroCrossEvent('z-sin(beta)', -1, touchdown_args,
                                             ['z'], ['beta'],
                                             targetlang=targetlang)

    flight_args = getFlightArgs(algparams, touchdown_ev, pars, abseps)
    return (stance_args, flight_args)

