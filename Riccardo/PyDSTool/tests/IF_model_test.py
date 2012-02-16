#! /usr/local/bin/python
"""Integrate-and-Fire neuron model test.

The model exhibits a square-pulse spike on reaching a firing
threshold. The model utilizes an 'absolute refractory period' and a
discontinuous state reset from refractoriness using features of the
hybrid system simulator.

Uses events to detect threshold crossing accurately, and separate
"vector fields" are simulated for sub-threshold and spiking behaviour,
using individual Generators. The Generators are brought together using
the Model class as a hybrid dynamical system.

    Robert Clewley, March-August 2005.
"""

from PyDSTool import *
from time import clock

# ---------------------------------------------------------------------------

print '-------- Model test'
allgen_names = ['leak', 'spike']

# 'excited' is an internal variable of the model, and is used to
# ensure that the.compute() method can determine which Generator
# to start the calculation with
leak_event_args = {'name': 'threshold',
                   'eventtol': 1e-3,
                   'eventdelay': 1e-5,
                   'starttime': 0,
                   'active': True,
                   'term': True,
                   'precise': True}
leak_thresh_ev = Events.makePythonStateZeroCrossEvent('V', -60, 1,
                                                 leak_event_args)
leak_args = {'pars': {'I': 1.3, 'gl': 0.1, 'vl': -67},
              'xdomain': {'V': [-100,50], 'excited': 0.},
              'varspecs': {'V':       "I - gl*(V-vl)",
                              'excited': "0"
                              },
              'algparams': {'init_step': 0.02},
              'events': leak_thresh_ev,
              'abseps': 1.e-7,
              'name': 'leak'}
DS_leak = Generator.Vode_ODEsystem(leak_args)


# spike length parameter 'splen' must be contained within 'tdomain' in
# order to get a fully-formed square-pulse `spike`
spike_args = {'tdomain': [0.0, 1.5],
                'varspecs': {'V': "if(t<splen,50,-95)", 'excited': "1."},
                'pars': {'splen': 0.75},
                'xdomain': {'V': [-97, 51], 'excited': 1.},
                'name': 'spike'}
DS_spike = Generator.ExplicitFnGen(spike_args)

# test discrete state mapping that is used at changes of vector field
# after terminal events
# must set excited to 0 in order to meet bounds requirements of IF gen
epmapping = EvMapping({"xdict['V']": "xdict['V']+15",
                       "xdict['excited']": "0"})

# build model object from individual DS's
DS_leak_info = makeGenInfoEntry(DS_leak, allgen_names, [('threshold', 'spike')])
DS_spike_info = makeGenInfoEntry(DS_spike, allgen_names,
                                 [('time', ('leak', epmapping))])
genInfoDict = makeGenInfo([DS_leak_info, DS_spike_info])

IFmodel = Model.Model({'name': 'IF_fit', 'genInfo': genInfoDict})
icdict = {'V': -80., 'excited': 0.}

print "Computing trajectory...\n"
start = clock()
IFmodel.compute(trajname='onespike',
                  tdata=[0, 30],
                  ics=icdict)
print '... finished in %.3f seconds.\n' % (clock()-start)

print 'Preparing plot to show non-identity mapping of epoch state transitions'
plotData = IFmodel.sample('onespike', ['V'], 0.02)
pylab.ylabel('V')
pylab.xlabel('t')
vline=plot(plotData['t'], plotData['V'])

print "\n\nInformation about Model's generators:\n"
info(IFmodel.query('generators'))

show()
