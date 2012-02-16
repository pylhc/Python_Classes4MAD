#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      1 Mar 2005
#$Revision: 1.0 $
"""
    Parameter Estimation tests #1.

    Robert Clewley, February 2005.
"""

# PyDSTool imports
from PyDSTool import *
from PyDSTool.ParamEst import LMpest

# Other imports
import time


# ----------------------------------------------------------------


## Initialize some values
trange = [0,2.5]
xic = {'w':3.0}
goalpars = {'k':2, 'a':25}

## Prepare sine-wave input to represent 'noise' in the goal trajectory
fnArgs = args(varspecs={'s': "sin(t*speed)"},
          xdomain={'s': [-1, 1]},
          pars={'speed': 18},
          name='sine')
sin_input=Generator.ExplicitFnGen(fnArgs)
sin_input_traj = sin_input.compute('noise')


## Prepare goal ODE trajectory
print 'Preparing goal trajectory from perturbed system (with k = ', \
    goalpars['k'], ') ...'

xfn_goalstr = "50-k*w+a*(2-t)*sin_input"

goalDSargs = args(algparams={'init_step':0.02, 'strictdt':True},
              tdata=trange,
              pars=goalpars,
              varspecs={'w':xfn_goalstr},
              xdomain={'w':[0, 1000]},
              inputs={'sin_input' : sin_input_traj.variables['s']},
              checklevel=1,
              ics=xic,
              name='ODEtest'
              )
goalODE = Generator.Vode_ODEsystem(goalDSargs)
goaltraj = goalODE.compute('goal')
print '... finished.\n'


## Get plot data for goal orbit
tplotData = arange(min(trange), max(trange), 0.02)
wplotData = goaltraj(tplotData,['w'])
goalleg = 'w'+' for k = '+str(goalpars['k'])


## Parameter estimation
print 'Estimating parameter k for fit (assuming initial condition for w)'
print 'Goal value is k = ', goalpars['k'], " ..."

xfn_str = "50-k*w"
testDSargs = args(algparams={'init_step':0.02, 'strictopt':True},
              varspecs={'w':xfn_str},
              xdomain={'w':[0, 1000]},
              tdata=trange,
              pars={'k':0.1},
              checklevel=2,
              ics=xic,
              name='testODEsol'
              )
testODE = Generator.Vode_ODEsystem(testDSargs)

genInfoEntry = makeGenInfoEntry(testODE, ['testODEsol'])
genInfoDict = makeGenInfo([genInfoEntry])
modelArgs_par = {'genInfo': genInfoDict,
             'name': 'test_model_par',
             'ics': xic}

testModel_par = Model.Model(modelArgs_par)
pest_par = LMpest(freeParams=['k'],
                 depVarNames=['w'],
                 testModel=testModel_par,
                 goalTraj=goaltraj,
                 trange=trange
                )

start_time = time.clock()
pestData_par = pest_par.run(verbose=True)
bestFitModel_par = pestData_par['sys_fit']
print '... finished in %.4f seconds\n' % (time.clock()-start_time)


## Initial condition estimation
print "Estimating initial condition for w (assuming k is correct)"
print "Goal value is w(0) = ", xic['w'], " ..."
modelArgs_ic = {'genInfo': genInfoDict,
             'name': 'test_model_ic',
             'ics': {'w': 0.0}}
testModel_ic = Model.Model(modelArgs_ic)
pest_ic = LMpest(freeParams=['w'],
                 depVarNames=['w'],
                 testModel=testModel_ic,
                 goalTraj=goaltraj,
                 trange=trange
                )
pestData_ic = pest_ic.run(verbose=True)
bestFitModel_ic = pestData_ic['sys_fit']
print '... finished'


## Finish preparing plots
print '\nPreparing plots'
west_plotData_par = bestFitModel_par.sample('par_est', dt=0.02,
                                                           tlo=min(trange),
                                                           thi=max(trange),
                                                precise=True)
west_plotData_ic = bestFitModel_ic.sample('par_est', dt=0.02,
                                                         tlo=min(trange),
                                                         thi=max(trange),
                                                precise=True)

pylab.ylabel('w')
pylab.xlabel('t')
goalline=plot(tplotData, wplotData['w'])
estline_par = plot(west_plotData_par['t'], west_plotData_par['w'])
estleg_par = 'w'+' (estimated) for k = '+str(pestData_par['pars_fit']['k'])
estline_ic = plot(west_plotData_ic['t'], west_plotData_ic['w'])
estleg_ic = 'w for estimated initial condition at w(0) = ' + \
                    str(pestData_ic['pars_fit']['w'])
pylab.legend([goalline, estline_par, estline_ic],
             [goalleg, estleg_par, estleg_ic],
             'lower right')
show()

print 'Done'
