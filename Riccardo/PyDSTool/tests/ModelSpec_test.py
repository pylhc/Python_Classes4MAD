#! /usr/local/bin/python
# ModelSpec_test1.py
# Robert Clewley, September 2005

from PyDSTool import *
from PyDSTool.Toolbox.neuralcomp import *
from copy import copy

targetGen = 'Dopri_ODEsystem'

# -------------------------------------------------------------------------

if targetGen == 'Vode_ODEsystem':
    targetlang='python'
else:
    targetlang='c'

# -------------------------------------------------------------------------

v = Var(voltage)
ma = 0.32*(v+54)/(1-Exp(-(v+54.)/4))
mb = 0.28*(v+27)/(Exp((v+27.)/5)-1)
ha = .128*Exp(-(50.+v)/18)
hb = 4/(1+Exp(-(v+27.)/5))
channel_Na1 = makeChannel_rates('Na', voltage, 'm', False, ma, mb, 3,
                                'h', False, ha, hb, 1, vrev=50, g=100)

na = .032*(v+52)/(1-Exp(-(v+52.)/5))
nb = .5*Exp(-(57.+v)/40)
channel_K1 = makeChannel_rates('K', voltage, 'n', False, na, nb, 4, vrev=-100, g=80)

channel_Ib1 = makeBiasChannel('Ib', 2.1)
channel_Lk1 = makeChannel_rates('Lk', vrev=-67, g=0.1)

HHcell1 = makeSoma('cell1', channelList=[channel_Lk1, channel_Ib1,
                   channel_Na1, channel_K1], C=1.0)
channel_Ib2 = makeBiasChannel('Ib', 2.2)
channel_Lk2 = copy(channel_Lk1)
channel_Na2 = copy(channel_Na1)
channel_K2 = copy(channel_K1)
##HHcell2 = makeSoma('cell2', channelList=[channel_Lk2, channel_Ib2,
##                   channel_Na2, channel_K2], C=1.0)
HHcell2 = copy(HHcell1)
HHcell2.rename('cell2')

syn12 = connectWithSynapse('s12', 'inh', HHcell1, HHcell2, g=1)
syn21 = connectWithSynapse('s21', 'inh', HHcell2, HHcell1, g=0.5)

HHnet = makePointNeuronNetwork('HHnet_test1', [HHcell1, HHcell2, syn12, syn21])
print "Compiled network specification"

# ensure HHnet has a copy of its flattened spec for use in making event
# (otherwise it is not necessary)
HHnet.flattenSpec()

# try building an event that picks out when RHS of cell1's Voltage eqn is 0
# i.e. when dV/dt=0

stat_ev_args = {'name': 'cell1_stat',
               'eventtol': 1e-3,
               'eventdelay': 1e-3,
               'starttime': 0,
               'term': False
                }
# stationary event => v = 0
stat_ev = Events.makeZeroCrossEvent(HHnet.flatSpec['vars']['cell1_V'],
                        0, stat_ev_args, targetlang=targetlang,
                        flatspec=HHnet.flatSpec)

alg_args = {'init_step':0.15}
ic_args_net = {'cell1.V':-70.0, 'cell1.Na.m': 0,
               'cell1.Na.h': 1, 'cell1.K.n': 0,
               'cell2.V':-75.0, 'cell2.Na.m': 0,
               'cell2.Na.h': 1, 'cell2.K.n': 0,
               's12.s_cell1_cell2': 0, 's21.s_cell2_cell1': 0}
modelC_net = ModelConstructor('HHnet_model',
                          generatorspecs={'HHnet_test1': {'modelspec': HHnet,
                                                'target': targetGen,
                                                'algparams': alg_args}},
                          indepvar=('t',[0,100]),
                          parvalues={'cell1.chan_s21.g': 0.3,
                                     'cell2.chan_s12.g': 0.35},
                          eventtol=1e-5)
modelC_net.addEvents('HHnet_test1', stat_ev)
HHmodel_net = modelC_net.getModel()
print "Instantiated model for target ODE solver "+targetGen

verboselevel = 2
print "Computing trajectory using verbosity level %d..."%verboselevel
HHmodel_net.compute(trajname='test',
                     tdata=[0, 100],
                     ics=ic_args_net,
                     verboselevel=verboselevel)

print "Compiling data for plot"
v_dat = HHmodel_net.sample('test', dt=0.15)
pylab.figure()
v1line = plot(v_dat['t'], v_dat['cell1.V'])
v2line = plot(v_dat['t'], v_dat['cell2.V'])

show()
