#! /usr/local/bin/python
#Author:    Robert Clewley
#Date:      1 Mar 2005
#$Revision: 1.0 $
"""
    Parameter Estimation tests #2.

    Robert Clewley, March 2005.
"""

# PyDSTool imports
from PyDSTool import *
from PyDSTool.ParamEst import LMpest
import HH_model, IF_squarespike_model


# ----------------------------------------------------------------

par_args_HH = {'gna': 100, 'gk': 80, 'gl': 0.1,
            'vna': 50, 'vk': -100, 'vl': -67,
            'I': 1.35, 'C': 1.0}
# deliberately set I not quite 1.3, as used for IF neuron
ic_args_HH = {'v':-70.0, 'm': 0, 'h': 1, 'n': 0}
HH = HH_model.makeHHneuron('goalHH', par_args_HH, ic_args_HH)

HH.set(tdata=[0,15])

HH_traj = HH.compute('test')
HH_sampleData = {}
HH_sampleData['t'] = []
HH_sampleData['v'] = []
sample_dt = 0.06
count = 0
countlim = 5
print "Generating non-uniform samples from HH orbit..."
tsamples = arange(0, 14, sample_dt)
vsamples = HH_traj(tsamples, ['v']).toarray()
for i in range(len(tsamples)):
    t = tsamples[i]
    v = vsamples[i]
    if v > -57:
        HH_sampleData['t'].append(t)
        HH_sampleData['v'].append(v)
    else:
        # reduce sample rate for non-spiking region
        count += 1
        if count == countlim:
            HH_sampleData['t'].append(t)
            HH_sampleData['v'].append(v)
            count = 0
print "... done"
        

tableArgs = {'tdata': HH_sampleData['t'],
             'ics': {'v': HH_sampleData['v']},
             'name': 'HH_data'}
HH_DataTable = Generator.LookupTable(tableArgs)
tmesh_par = HH_sampleData['t']

par_args_linear = {'I': 1.3, 'gl': 0.1, 'vl': -67, 'threshval': -60, 'C': 1.0}
par_args_spike = {'splen': 1.0}

## Parameter estimation for firing threshold
icdict = {'v': -70.0, 'excited': 0}
IFmodel_thr = IF_squarespike_model.makeIFneuron('IF_thr_fit', par_args_linear,
                                par_args_spike, icdict=icdict)

# un-fitted IF trajectory
IFmodel_thr.compute(trajname='orig',
                        tdata=[0,tmesh_par[-1]],
                        ics={'v':-70, 'excited':0})
orig_pdata = IFmodel_thr.sample('orig', ['v'], 0.1)


HH_event_args = {'name': 'HH_zerothresh',
               'eventtol': 1e-2,
               'eventdelay': 1e-3,
               'starttime': 0,
               'active': True}
HH_thresh_ev = Events.makePythonStateZeroCrossEvent('v', 0, 1, HH_event_args,
                                               HH_traj.variables['v'])


result = HH_thresh_ev.searchForEvents((0, 15))
HH_spike_t = result[0][0]
print "HH spike time found at ", HH_spike_t

def residual_fn_thr(s, p, tmesh, HH_spike_t):
    s.modelArgs[s.parTypeStr][s.freeParNames[0]] = p
    s.testModel.compute(**s.modelArgs)
    # time of spike
    tparts = s.testModel.getTrajTimePartitions('par_est')
    if len(tparts) == 1:
        raise RuntimeError, "IF event not found"
    cost = abs(HH_spike_t-tparts[1][1])
    print "cost =", cost
    return cost


pest_thr = BoundMin(freeParams=['threshval'],
                 depVarNames=['v'],
                 testModel=IFmodel_thr,
                 goalTraj=HH_DataTable,
                 trange=[0,15],
                 residual_fn=residual_fn_thr
                )

pestData_thr = pest_thr.run(parConstraints=[-63,-58],
                            xtol=5e-3,
                            tmesh=tmesh_par,
                            extra_args = (HH_spike_t,),
                            verbose=True)


## Parameter estimation for spike length
print "\nParam est. for spike length ..."
if not pestData_thr['success']:
    raise RuntimeError, "Failure: will not continue"

thresh_fit = pestData_thr['pars_fit']['threshval']

par_args_linear = {'I': 1.3, 'gl': 0.1, 'vl': -67, 'threshval': thresh_fit,
                   'C': 1.0}
par_args_spike = {'splen': 1.0}

HH_datatable_traj = HH_DataTable.compute('goaltraj')

# find closest (t, v) point for i.c. near spike
ic_not_found = True
tmesh_ic = []
for i in xrange(len(HH_sampleData['t'])):
    t = HH_sampleData['t'][i]
    if t >= 7.0 and t < 11:
        tmesh_ic.append(HH_sampleData['t'][i])
        if ic_not_found:
            t_ic = t
            v_ic = HH_datatable_traj(t, ['v'])
            ic_not_found = False
    if t >= 11:
        break


IFmodel_splen = IF_squarespike_model.makeIFneuron('IF_splen_fit', par_args_linear,
                                par_args_spike, icdict={'v':-70, 'excited':0})

## test IF trajectory
IFmodel_splen.compute(trajname='test', tdata=[0,t_ic])
IF_ic = IFmodel_splen('test', t_ic, ['v'])

IFmodel_splen.set(ics={'v': IF_ic})
pest_splen = BoundMin(freeParams=['splen'],
                 depVarNames=['v'],
                 testModel=IFmodel_splen,
                 goalTraj=HH_datatable_traj,
                 trange=[t_ic,12]
                )

pestData_splen = pest_splen.run(xtol=0.01, parConstraints=[0.2,1.0],
                                tmesh=tmesh_ic, verbose=True)

IFmodel_splen.set(pars={'splen': pestData_splen['pars_fit']['splen'],
                        'threshval': thresh_fit})

IFmodel_splen.compute(trajname='disp',
                      tdata=[0,15],
                      ics={'v':-70, 'excited':0})

## Plot data
print "Acquiring plot data"
origline=plot(orig_pdata['t'], orig_pdata['v'])
origleg = "Un-fitted IF orbit"
IF_sampleData = []
for t in HH_sampleData['t']:
    IF_sampleData.append(IFmodel_splen('disp', t, ['v']))
pylab.ylabel('w')
pylab.xlabel('t')
goalline=plot(HH_sampleData['t'], HH_sampleData['v'], 'bo')
goalleg = 'HH reference'
estline_splen = plot(HH_sampleData['t'], IF_sampleData, 'k-',
                         linewidth=2)
estleg_splen = 'IF spike thresh & width fitted'
pylab.legend([goalline, estline_splen, origline],
             [goalleg, estleg_splen, origleg],
             'lower left')
show()
print 'Done'
