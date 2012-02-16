#! /usr/local/bin/python

# Test object deletion, and corresponding deletion of dynamically created
# class methods that are no longer referenced

from PyDSTool import *

var1 = Variable(Pointset(coordarray = array(range(10), Float)*0.1,
                        indepvararray = array(range(10), Float)*0.5
                        ), name='v1')

print var1
del var1
try:
    print var1
    raise RuntimeError, "Variable deletion unsuccessful"
except NameError:
    print "Simple variable deletion successful"
