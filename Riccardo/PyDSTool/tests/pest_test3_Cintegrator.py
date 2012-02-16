#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      1 Mar 2005
#$Revision: 1.0 $
"""
    Parameter Estimation tests #3.

    Robert Clewley, March 2005.
"""

# PyDSTool imports
from PyDSTool import *
from PyDSTool.ParamEst import LMpest, BoundMin
import HH_model_Cintegrator as HH_model

# Other imports
from time import clock

# ----------------------------------------------------------------

par_args_HH_goal = {'gna': 100, 'gk': 80, 'gl': 0.1,
            'vna': 50, 'vk': -100, 'vl': -67,
            'I': 1.3, 'C': 1.0}
ic_args_HH = {'v':-68, 'm': 0, 'h': 1, 'n': 0}

HH_event_args = args(name='HH_zerothresh',
               eventtol=1e-2,
               eventdelay=1e-3,
               starttime=0,
               active=True)
HH_thresh_ev = Events.makeZeroCrossEvent('v', 1, HH_event_args, ['v'],
                                         targetlang='c')

HH_goal = HH_model.makeHHneuron('goalHH', par_args_HH_goal, ic_args_HH,
            evs=HH_thresh_ev,
            extra_terms='-0.04*(sin(9.1*t)*cos(2.6*t)+sin(5.1119*t+2))*(v-60)')

HH_goal.set(tdata=[0,15])

goaltraj = HH_goal.compute('goalHHtraj')


HH_spike_t = HH_goal.getEventTimes()['HH_zerothresh'][0]
print "HH spike time found at ", HH_spike_t


## Set up test HH model
par_args_HH_test = {'gna': 100, 'gk': 80, 'gl': 0.12,
            'vna': 50, 'vk': -100, 'vl': -70,
            'I': 1.34, 'C': 1.0}
# Note that I is not the same as that for goal, even though we're not
# optimizing this parameter. Increasing I from original 1.3 to 1.34
# causes slow convergence. Without sp_adjust having a factor of 15 I=1.34
# won't converge.

DS_event_args = args(name='threshold',
           eventtol=1e-2,
           eventdelay=1e-3,
           starttime=0,
           active=True,
           term=False,
           precise=True)
thresh_ev = Events.makeZeroCrossEvent('v', 1, DS_event_args, ['v'],
                                      targetlang='c')
HH_test = HH_model.makeHHneuron('testHH', par_args_HH_test, ic_args_HH,
                                thresh_ev)

HH_test.set(tdata=[0,15], algparams={'atol':1e-5,'rtol':1e-5})

# Make model out of HH DS
HH_test_model = embed(HH_test, ic_args_HH)

# time sample for parameter estimation only (not for plotting)
dt = 0.4


## Parameter estimation
print 'Estimating pars gl and vl for fit'
print 'Goal values are vl, gl =', par_args_HH_goal['vl'], ', ', \
            par_args_HH_goal['gl'], ' ...'


def residual_fn(s, p, tmesh, HH_spike_t):
    for i in xrange(s.numFreePars):
        s.modelArgs[s.parTypeStr[i]][s.freeParNames[i]] = p[i]
    s.testModel.compute(**s.modelArgs)
    try:
        evtime = s.testModel.getTrajEventTimes('par_est', 'threshold')[0]
    except ValueError:
        evtime = 20
    gl = s.testModel.pars['gl']
    vl = s.testModel.pars['vl']
    gl_adjust = 30*(0.0001/gl**2 + 0.0001/(gl-0.15)**2)
    vl_adjust = 0.001*(vl+66)**4
    sp_adjust = 20*(evtime-HH_spike_t)**2
    vec = vl_adjust + gl_adjust + sp_adjust + \
           1/(1+10*sp_adjust) * (s.goalTraj(tmesh,s.depVarNames).toarray() - \
                       s.testModel(s.trajname, tmesh, s.depVarNames).toarray())
    adjusted = sum(vec**2.0)
##    print vl, gl, vl_adjust, gl_adjust, sp_adjust
    print "adjusted error = ", adjusted, "\n"
    return vec


pest_pars = LMpest(freeParams=['vl', 'gl'],
                 depVarNames=['v'],
                 testModel=HH_test_model,
                 goalTraj=goaltraj,
                 trange=[0,15],
                 residual_fn=residual_fn
                )

start = clock()
pestData_par = pest_pars.run(parDict={'ftol':1e-4,
                                      'xtol':1e-4,
                                      'args':(HH_spike_t,)},
                             verbose=True)
print '... finished in %.3f seconds\n' % (clock()-start)


## Finish preparing plots
print '\nPreparing plots'
disp_dt = 0.05
plotData_goal = goaltraj.sample(['v'], disp_dt, precise=True)
goalleg = "v original"
plotData_par = HH_test_model.sample('par_est', ['v'], disp_dt, precise=True)

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
