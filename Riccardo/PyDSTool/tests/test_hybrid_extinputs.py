#! /usr/local/bin/python
"""
    Test hybrid system with external inputs and mixed classes of
    integrator.

    Robert Clewley, September 2006
"""
from PyDSTool import *
from copy import copy

# ------------------------------------------------------------

timeData = linspace(0, 100, 200)
sindata = sin(timeData)+2*cos(timeData/2+pi/5)
xData = makeDataDict(['in'], [sindata])
my_input = InterpolateTable({'tdata': timeData,
                              'ics': xData,
                              'name': 'interp1d',
                              'method': 'linear',
                              'checklevel': 1,
                              'abseps': 1e-5
                              }).compute('interp')


DSargs = {}
DSargs['varspecs'] = {'x': 'in+1', 'inval': 'in'}
DSargs['vars'] = ['x']  # inval is an aux variable
DSargs['algparams'] = {'init_step': 0.01, 'refine': 0, 'max_step': 0.01,
                       'rtol': 1e-4, 'atol': 1e-4}
DSargs['checklevel'] = 2
DSargs['ics'] = {'x': 0}
DSargs['inputs'] = {'in': my_input.variables['in']}
DSargs['name'] = 'gen1'
DSargs['events'] = Events.makeZeroCrossEvent('x-5', 1,
                                    {'name': 'ev_togen2',
                                     'eventtol': 1e-4,
                                     'precise': True,
                                     'term': True}, ['x'], targetlang='c')

gen1 = Generator.Radau_ODEsystem(DSargs)

DSargs2=copy(DSargs)
DSargs2['name'] = 'gen2'
DSargs2['varspecs'] = {'x': '-abs(in)/2', 'inval': 'in'}
DSargs2['algparams'] = {'init_step': 0.01,
                       'rtol': 1e-4, 'atol': 1e-4}
DSargs2['events'] = Events.makeZeroCrossEvent('x+5', -1,
                                    {'name': 'ev_togen1',
                                     'eventtol': 1e-4,
                                     'precise': True,
                                     'term': True}, ['x'])
gen2 = Generator.Vode_ODEsystem(DSargs2)

allnames = ['gen1','gen2']
g1Info = makeGenInfoEntry(gen1, allnames, [('ev_togen2', 'gen2')])
g2Info = makeGenInfoEntry(gen2, allnames, [('ev_togen1', 'gen1')])
genInfoDict = makeGenInfo([g1Info, g2Info])

m = Model.Model(name='test_inputs', genInfo=genInfoDict)


print 'Integrating...'
m.compute(trajname='test', tdata=[0,100], ics={'x':5.1})

plotData = m.sample('test')

print "Testing output"
diffs = []
for t in range(0, 100, 10):
    diffs.append(abs(my_input(t)-m('test',t)['inval']))
assert all(array(diffs)<1e-2), \
       "Error in hybrid system's version of external input too large"

print 'Preparing plot'

pylab.ylabel('x')
pylab.xlabel('t')
xline=plot(plotData['t'], plotData['x'])
iline=plot(plotData['t'], plotData['inval'])
iline_true=plot(timeData, sindata)
show()
