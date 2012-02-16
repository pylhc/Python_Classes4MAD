#! /usr/local/bin/python

from PyDSTool import *

objs = loadObjects('temp_HH.pkl')
HH = objs[0]
HHtraj1 = objs[1]
plotData1 = HHtraj1.sample(dt=0.05, eventdata=HH.getEventTimes())
evt1=HH.getEventTimes()['thresh_ev'] # HH was last used to compute HHtraj1 in HH_model.py

HH.set(tdata=[0, 10], pars={'C': 1.2})
HHtraj2 = HH.compute('test2')
plotData2 = HHtraj2.sample(dt=0.05, eventdata=HH.getEventTimes())
evt2=HH.getEventTimes()['thresh_ev']

yaxislabelstr = 'v'
pylab.ylabel(yaxislabelstr)
pylab.xlabel('t')

vline1=plot(plotData1['t'], plotData1['v'])
vline2=plot(plotData2['t'], plotData2['v'])

plot(evt1, HHtraj1(evt1, 'v'), 'ro')
plot(evt2, HHtraj2(evt2, 'v'), 'ro')
show()
