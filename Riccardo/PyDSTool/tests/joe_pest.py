#! /usr/local/bin/python
"""A parameter estimation problem inspired by Joe Tien's research.

   Requires files 'ttraj_prebotfast', 'xtraj_prebotfast',
   't_prefastout', and 'x_prefastout'. ttraj_prebotfast (time points)
   and xtraj_prebotfast (state variables) contain a reference (target)
   trajectory, with corresponding parameter values h=0.3, gl=40,
   vl=-38., gnap=2.8. t_prefastout and x_prefastout contain the
   initial guess.  Here I'm optimizing over two variables, gl and vl.
   The parameter values for the initial guess are h=0.3, gl=38,
   vl=-35, gnap=2.8.

   Robert Clewley, March 2005.
"""

from PyDSTool import *
from PyDSTool.ParamEst import LMpest

import random
from copy import copy
from time import clock


# ------------------------------------------------------------


def makeTestDS(name, par_args, ic_args, tdomain, evs=None):

    # Variables: v, n
    # Parameters: h, gl, vl, gnap

    fnspecs = {'m': (['v'], '1.0/(1.0+exp((v-th_m)/s_m))'),
      'an': (['v'], '1.0/(2.0*tau_n)*exp(-(v-th_n)/(2.0*s_n))'),
      'bn': (['v'], '1.0/(2.0*tau_n)*exp((v-th_n)/(2.0*s_n))'),
      'mnap': (['v'], '1.0/(1.0+exp((v-th_mnap)/s_mnap))')
      }

    vstr = '-1.0/C*(gna*m(v)*m(v)*m(v)*(1-n)*(v-vna)+gk*n*n*n*n*(v-vk)+' +\
                       'gnap*mnap(v)*h*(v-vna)+gl*(v-vl))'
    nstr = 'an(v)*(1.0-n)-bn(v)*n'

    DSargs = args()
    DSargs.tdomain = tdomain
    DSargs.pars = par_args
    DSargs.varspecs = {'v': vstr, 'n': nstr}
    DSargs.fnspecs = fnspecs
    DSargs.xdomain = {'v': [-130, 70], 'n': [0,1]}
    DSargs.algparams = {'init_step':0.02}
    DSargs.ics = ic_args
    DSargs.name = name
    return Generator.Vode_ODEsystem(DSargs)


# ------------------------------------------------------------

pars = {'C': 21.0, 'gna': 28.0, 'gk': 11.2, 'vna': 50.0, 'th_m': -34.0,
      's_m': -5.0, 'vk': -85.0, 'th_n': -29.0, 's_n': -4.0, 'tau_n': 10.0,
      'th_mnap': -40.0, 's_mnap': -6.0}

tdomain = [0,11]
icdict = {'v': -24.228, 'n': 0.43084}

# make reference system
# in data system version, these "known" pars are not used except to
# show the goal values to the user. they are not part of the param. est. calc.
est_pars_ref = {'h':0.3, 'gl': 40, 'vl': -38., 'gnap': 2.8}
pars_ref = copy(pars)
pars_ref.update(est_pars_ref)
##
##refDS = makeTestDS('ref', pars_ref, icdict, tdomain)
##reftraj = refDS.compute('ref')

allValDict = utils.importPointset('xtraj_prebotfast',
                                  t='ttraj_prebotfast')
xDataDict = utils.makeDataDict(['v', 'n'], allValDict['vararray'])
tmesh = allValDict['t'].tolist()
refDS = Generator.LookupTable({'tdata': tmesh,
                                'ics': xDataDict,
                                'name': 'ref'
                                })
reftraj = refDS.compute('ref')

# make test system
est_pars_test = {'h':0.3, 'gl': 38, 'vl': -35, 'gnap': 2.8}
pars_test = copy(pars)
pars_test.update(est_pars_test)
testDS = makeTestDS('ref', pars_test, icdict, tdomain)
testModel = embed(testDS, icdict)

est_parnames = ['gl', 'vl']

# parameter estimation
print 'Starting parameter estimation'
print 'Goal pars are gl = ', est_pars_ref['gl'], ' vl = ', est_pars_ref['vl']
pest_pars = LMpest(freeParams=est_parnames,
                 depVarNames=['v'],
                 testModel=testModel,
                 goalTraj=reftraj,
                 trange=tdomain,
                )

start = clock()
pestData_par = pest_pars.run(parDict={'ftol':3e-3,
                                      'xtol':1e-3},
                             tmesh=tmesh,
                             verbose=True)
print '  ... finished in %.3f seconds.\n' % (clock()-start)

## Finish preparing plots
print '\nPreparing plots'
disp_dt = 0.05
##plotData_goal = reftraj.sample(['v'], disp_dt)
goalleg = "v original"
plotData_par = testModel.sample('par_est', ['v'], disp_dt, precise=True)

pylab.ylabel('v')
pylab.xlabel('t')
##goalline=plot(plotData_par['t'], plotData_goal['v'])
goal_v = reftraj(tmesh, 'v').toarray()
goalline = plot(tmesh, goal_v, 'ok')
estline = plot(plotData_par['t'], plotData_par['v'])
estleg = 'v estimated'

pylab.legend([goalline, estline],
             [goalleg, estleg],
             'lower left')
show()
