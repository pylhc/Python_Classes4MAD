#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      1 Mar 2005
#$Revision: 1.0 $
"""
    Parameter Estimation tests #4.

    Variant of #3 but with sodium conductance pars being optimized

    Robert Clewley, March 2005.
"""

# PyDSTool imports
from PyDSTool import *
from PyDSTool.ParamEst import LMpest, BoundMin
import HH_model

from time import clock

# ----------------------------------------------------------------

par_args_HH_goal = {'gna': 100, 'gk': 80, 'gl': 0.1,
            'vna': 50, 'vk': -100, 'vl': -67,
            'I': 1.3, 'C': 1.0}
ic_args_HH = {'v':-70, 'm': 0, 'h': 1, 'n': 0}
HH_goal = HH_model.makeHHneuron('goalHH', par_args_HH_goal, ic_args_HH)
##                        extra_terms='-0.04*(sin(9.1*t)*cos(2.6*t)+sin(5.1119*t+2))*(v-60)')
HH_goal.set(tdata=[0,20])
goaltraj = HH_goal.compute('goalHHtraj')

HH_event_args = args(name='HH_zerothresh',
               eventtol=1e-3,
               eventdelay=1e-3,
               starttime=0,
               active=True,
               precise=True)
HH_thresh_ev = Events.makePythonStateZeroCrossEvent('v', -30, 1, 
                                   HH_event_args, goaltraj.variables['v'])


result = HH_thresh_ev.searchForEvents((0, 20))
HH_spike_t = result[0][0]
print "HH spike time found at ", HH_spike_t

## Set up test HH model
par_args_HH_test = {'gna': 95, 'gk': 81, 'gl': 0.12,
            'vna': 50, 'vk': -99, 'vl': -67.2,
            'I': 1.32, 'C': 1.0}
# Note that these params are not the same as that for goal, even though we're not
# optimizing them

DS_event_args = args(name='threshold',
           eventtol=5e-3,
           eventdelay=1e-3,
           starttime=0,
           active=True,
           term=False,
           precise=True)
thresh_ev = Events.makePythonStateZeroCrossEvent('v', -30, 1, DS_event_args)
HH_test = HH_model.makeHHneuron('testHH', par_args_HH_test, ic_args_HH,
                                thresh_ev)
HH_test.set(tdata=[0,20], algparams={'atol':1e-9,'rtol':1e-8,
                                               'min_step':1e-5})

# Make model out of HH DS
HH_test_model = embed(HH_test, ic_args_HH)

# time sample for parameter estimation only (not for plotting)
dt = 0.4


## Parameter estimation
print 'Estimating pars gna and vl for fit to non-identical HH cell'
print 'Goal values are gna, gl =', par_args_HH_goal['gna'], ', ', \
            par_args_HH_goal['gl'], ' ...'


def residual_fn(s, p, tmesh, HH_spike_t):
    for i in xrange(s.numFreePars):
        s.modelArgs[s.parTypeStr[i]][s.freeParNames[i]] = p[i]
    s.testModel.compute(**s.modelArgs)
    try:
        evtime = s.testModel.getTrajEventTimes('par_est', 'threshold')[0]
    except ValueError:
        evtime = 20
##    gl = s.testModel.pars['gl']
##    gna = s.testModel.pars['gna']
##    gl_adjust = 10*(0.0001/pow(gl,2) + 0.0001/pow(gl-0.15,2))
##    gna_adjust = 0.000001*pow(gna-90,4)
##    #  300*(1/pow(gna-40,2) + 1/pow(gna-140,2))
    sp_adjust = (evtime-HH_spike_t)**2
    # the following "adjusted" has not been vectorized
##    adjusted = sum([pow(gna_adjust + gl_adjust + sp_adjust + \
##                    8/(1+sp_adjust) *(s.goalTraj(t,s.depVarNames) - \
##                                    s.testModel(s.trajname,
##                                                          t,
##                                                          s.depVarNames)),2) \
##                                                            for t in tmesh])
##    print gna, gl, gna_adjust, gl_adjust, sp_adjust
##    print "adjusted error = ", adjusted, "\n"
##    return gna_adjust + gl_adjust + sp_adjust + \
##           8/(1+sp_adjust) * array([s.goalTraj(t,s.depVarNames)[0] - \
##                                    s.testModel(s.trajname,
##                                                          t,
##                                                          s.depVarNames)[0] \
##                                                            for t in tmesh])
    return s.goalTraj(tmesh,s.depVarNames).toarray() - \
                        s.testModel(s.trajname, tmesh, s.depVarNames).toarray()

pest_pars = LMpest(freeParams=['gna', 'gl'],
                 depVarNames=['v'],
                 testModel=HH_test_model,
                 goalTraj=goaltraj,
                 trange=[0,20],
                 residual_fn=residual_fn,
                 usePsyco=True
                )

t0=clock()
pestData_par = pest_pars.run(parDict={'ftol':1e-4,
                                      'xtol':1e-4,
                                      'args':(HH_spike_t,)},
                             verbose=True)
print '... finished in %.4f seconds\n'%(clock()-t0)


## Finish preparing plots
print '\nPreparing plots'
disp_dt = 0.05
plotData_goal = goaltraj.sample(['v'], disp_dt)
goalleg = "v original"
plotData_par = HH_test_model.sample('par_est', ['v'], disp_dt)

pylab.ylabel('v')
pylab.xlabel('t')
goalline=plot(plotData_goal['t'], plotData_goal['v'])
estline = plot(plotData_par['t'], plotData_par['v'])
estleg = 'v estimated'

pylab.legend([goalline, estline],
             [goalleg, estleg],
             'lower left')
show()

print 'Done'
