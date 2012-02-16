#! /usr/local/bin/python
"""Library function to build Integrate and fire model with square-pulse
spike, as a hybrid system.

   Robert Clewley, March 2005.
"""

from PyDSTool import *
from time import clock

# ----------------------------------------------------------------

def makeLinearLeak(name, rhs, par_args, inputs):
    rhs_full = {'v': "(I - gl*(v-vl))/C",
                'excited': "0"}
    if rhs is not None:
        assert isinstance(rhs, dict)
        rhs_full.update(rhs)
##    rhs_full['testaux'] = "globalindepvar(t)-50"
    for parname in ['threshval', 'vl', 'gl', 'I', 'C']:
        assert parname in par_args, "Essential pars missing"
    DS_event_args = {'name': 'threshold',
               'eventtol': 1e-2,
               'eventdelay': 1e-3,
               'starttime': 0,
               'term': True
                }
    thresh_ev = Events.makePythonStateZeroCrossEvent('v', "p['threshval']",
                                                  1, DS_event_args)

    DS_args = {'pars': par_args,
               'varspecs': rhs_full,
               'xdomain': {'v': [-120,50], 'excited': 0.},
               'algparams': {'init_step': 0.1},
               'events': thresh_ev,
               'name': name}

    if inputs != {}:
        DS_args['inputs'] = inputs

    return Generator.Vode_ODEsystem(DS_args)


def makeSpike(name, par_args):
    # spike length parameter 'splen' must be contained within 'tdomain' in
    # order to get a fully-formed square-pulse `spike`
    DS_spike_args = {'tdomain': [0.0, 1.5],
            'varspecs': {'v': "if(t<splen,48,-97)", 'excited': "1"},
##                               'testaux': "globalindepvar(t)-50"}
            'pars': {'splen': par_args['splen']},
            'xdomain': {'v': [-98, 51], 'excited': 1.},
            'name': name}
    return Generator.ExplicitFnGen(DS_spike_args)


def makeIFneuron(name, par_args_linear, par_args_spike, rhs=None, inputs={},
                 icdict=None):
    allDSnames = ['linear', 'spike']

    epmapping = EvMapping({"xdict['excited']": "0"})
    # must set excited to 0 in order to meet bounds requirements of IF gen


    DS_linear = makeLinearLeak('linear', rhs, par_args_linear, inputs)
    DS_spike = makeSpike('spike', par_args_spike)
    DS_linear_info = makeGenInfoEntry(DS_linear, allDSnames,
                                           [('threshold', 'spike')])
    DS_spike_info = makeGenInfoEntry(DS_spike, allDSnames,
                                          [('time', ('linear', epmapping))])
    genInfoDict = makeGenInfo([DS_linear_info, DS_spike_info])

    # 'excited' is an internal variable of the model, and is used to
    # ensure that the.compute() method can determine which DS
    # to start the calculation with
    mod_args = {'name': name,
               'genInfo': genInfoDict}
    if icdict is not None:
        mod_args['ics'] = icdict

    IFmodel = Model.Model(mod_args)
    IFmodel.forceIntVars(['excited'])
    return IFmodel


# ----------------------------------------------------------------

if __name__=='__main__':
    print '-------- IF model test 1'

    par_args_linear = {'I': 1.3, 'gl': 0.1, 'vl': -67, 'threshval': -65, 'C': 1}
    par_args_spike = {'splen': 0.75}

    IFmodel = makeIFneuron('IF_fit', par_args_linear, par_args_spike)
    icdict = {'v': -80, 'excited': 0}

    start = clock()
    print 'Computing trajectory...'
    IFmodel.compute(trajname='onespike',
                        tdata=[0, 60],
                        ics=icdict,
                        verboselevel=0)
    print '\n... finished in %.3f seconds.\n' % (clock()-start)

    IFmodel.set(pars={'I': 1.0, 'threshval': -60})
    print 'Recomputing trajectory with new params...'
    IFmodel.compute(trajname='twospike',
                        tdata=[0, 60],
                        ics=icdict)


    print 'Preparing plot'
    plotData = IFmodel.sample('onespike', ['v'], 0.05)
    plotData2 = IFmodel.sample('twospike', ['v'], 0.05)
    pylab.ylabel('v')
    pylab.xlabel('t')
    vline=plot(plotData['t'], plotData['v'])
    vline2=plot(plotData2['t'], plotData2['v'])

    # For testing, should try this version with auxiliary variable
##    plotData = IFmodel.sample('onespike', ['v', 'testaux'], 0.05)
##    plotData2 = IFmodel.sample('twospike', ['v', 'testaux'], 0.05)
##    pylab.ylabel('v, testaux')
##    pylab.xlabel('t')
##    vline=plot(plotData['t'], plotData['v'])
##    vline2=plot(plotData2['t'], plotData2['v'])
##    aline=plot(plotData['t'], plotData['testaux'])

    print "Last point of hybrid trajectory: "
    print "IFmodel.getEndPoint('onespike') -->\n", IFmodel.getEndPoint('onespike')

    print "First point of hybrid trajectory: "
    print "IFmodel.getEndPoint('onespike', 0) -->\n", IFmodel.getEndPoint('onespike', 0)

    show()

    print "Testing IF hybrid model as mapping ..."
    num_parts = len(IFmodel.getTrajTimePartitions('twospike'))
    eventvals = IFmodel('twospike', range(0, num_parts+1), asmap=True)
    for i in range(0,num_parts+1):
        print "(t, v) at event(%i) = (%.4f, %.4f)" % (i, eventvals(i)('v'),
                                              eventvals(i)('t'))
    print "\nAlternative access to explicit event info using getTrajEvents(trajname) method:\n"
    print IFmodel.getTrajEvents('twospike')
