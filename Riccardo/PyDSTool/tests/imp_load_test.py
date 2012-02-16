#! /usr/local/bin/python

from PyDSTool import *

objs = loadObjects('temp_imp.pkl')
impgen = objs[0]
imptraj = objs[1]
print "Loaded objects successfully"
print "impgen._xdomain => ", impgen._xdomain
print "imptraj(-0.4) =>\n", imptraj(-0.4)

impgen.set(pars={'r':10.}, xdomain={'y': [-10,10]})
imptraj2 = impgen.compute('test2')
print "imptraj2(-0.4) =>\n", imptraj2(-0.4)
