#! /usr/local/bin/python
# PyDSTool/Generator/__init__.py
#Author:    Robert Clewley
#Date:      18 Jul 2006
#$Revision: 1.1.4 $
"""Trajectory generator classes.

   Robert Clewley, September 2005
"""

# -----------------------------------------------------------------------------
# Global version history
#
# 1.1.4, July 2006
#   Added support for supplying list of Parameter symbolic quantity definitions
#     in lieu of dictionary of name -> value pairs.
#
# 1.1.3, May 2006
#   Renamed 'computeTraj' and 'setPars' methods to 'compute' and 'set',
#     retaining 'computeTraj' as an alias for backwards compatibility.
#
# 1.1.2, Mar 2006
#   Added genDB registry object, findGenSubClasses.
#   Changed initialization arguments to accept ModelSpec Quantities.
#   Changed 'xdatadict' key to 'ics', and allow Point to be used.
#
# 1.1.1, Sep 2005
#   Added export of all concrete Gen class names for use by ModelConstructor.
#   Added Radau integrator.
#
# Local version histories for sub-classes ends at v1.1, 25 Jun 05 when this
#    sub-package split off so that each class has its own module and version
#    history inside it.
#
# 1.1, June 2005
#  Added functionality for ImplicitFnGen (without event handling)
#  Changed 'time' attribute to 'indepvariable' to be consistent with
#   ImplicitFnGen
#  Renamed _integparams attribute to _algparams
#  Moved copyVarDict utility to common.py
#  ODE -- self.indepvardomain no longer adjusted in Generator after terminal
#   event during integration: truncation of domain only passed to trajectory
#  Backwards integration in Dopri now supported
#
# 1.0.4, June 2005
#  Added support for multiple event tolerances, rtols etc.
#  Added 'pdomain' for optional parameter domain specifications -> becomes
#   the dictionary attribute parameterDomains
#
# 1.0.3, June 2005
#  Renamed methods / attributes with underscore separators to camel case names
#  Added some leading underscores to names intended to be private
#
# -----------------------------------------------------------------------------

from baseclasses import *
from ODEsystem import *
from Vode_ODEsystem import *
from Dopri_ODEsystem import *
from Radau_ODEsystem import *
from ADMC_ODEsystem import *
from ExplicitFnGen import *
from ImplicitFnGen import *
from LookupTable import *
from InterpolateTable import *
from MapSystem import *


def findGenSubClasses(superclass):
    """Find all Generator sub-classes of a certain class, e.g. ODEsystem."""
    assert isinstance(superclass, str), \
           "findGenSubClasses requires a string as the name of the class to search for subclasses."
    subclasslist = []
    sc = eval(superclass)
    for x in theGenSpecHelper.gshDB.keys():
        if compareClassAndBases(theGenSpecHelper.gshDB[x].genClass,sc):
            subclasslist.append(x)
    return subclasslist


# ----------------------------------------------------------------------

## TESTING

if __name__ == '__main__':
    import Events
    print '-------- Test: InterpolateTable'
    xnames = ['x1', 'x2']
    timeData = array([0.1, 1.1, 2.1])
    x1data = array([10.2, -1.4, 4.1])
    x2data = array([0.1, 0.01, 0.4])
    xData = makeDataDict(xnames, [x1data, x2data])
    itableArgs = {}
    itableArgs['tdata'] = timeData
    itableArgs['ics'] = xData
    itableArgs['name'] = 'interp'
    interptable = InterpolateTable(itableArgs)
    itabletraj = interptable.compute('itable')
    print "itabletraj(0.4, 'x1') = ", itabletraj(0.4, 'x1')
    print "                         (x1 should be 6.72)"
    print "itabletraj(1.1) = ", itabletraj(1.1)

    print '-------- Test: LookupTable'
    ltableArgs = copy(itableArgs)
    ltableArgs['name'] = 'lookup'
    lookuptable = LookupTable(ltableArgs)
    ltabletraj = lookuptable.compute('ltable')
    print "ltabletraj(1.1) = ", ltabletraj(1.1)
    print "             (x2 should be 0.01)\n"
    print "ltabletraj(0.4) ="
    try:
        ltabletraj(0.4)
    except ValueError, e:
        print "... error occurred: ", e
        print "(because 0.4 is not in the tdata array)"


    print '\n-------- Test: ODE system'

    fvarspecs = {"w": "k*w + a*itable + sin(t) + myauxfn1(t)*myauxfn2(w)",
               'aux_wdouble': 'w*2 + globalindepvar(t)',
               'aux_other': 'myauxfn1(2*t) + initcond(w)'}
    fnspecs = {'myauxfn1': (['t'], '2.5*cos(3*t)'),
                 'myauxfn2': (['w'], 'w/2')}
    DSargs = {'tdomain': [0.1,2.1],
              'pars': {'k':2, 'a':-0.5},
              'inputs': {'itable' : interptable.variables['x1']},
              'auxvars': ['aux_wdouble', 'aux_other'],
              'algparams': {'init_step':0.01, 'strict':False},
              'checklevel': 2,
              'name': 'ODEtest',
              'fnspecs': fnspecs,
              'varspecs': fvarspecs
              }
    testODE = Vode_ODEsystem(DSargs)
    print "params set => ", testODE.pars
    print "DS defined? => ", testODE.defined
    print "testODE.set(...)"
    testODE.set(ics={'w':3.0},
                    tdata=[0.11,2.1])
    print "testtraj = testODE.compute('test1')"
    testtraj = testODE.compute('test1')
    print "DS defined now? => ", testODE.defined
    print "testtraj(0.5) => ", testtraj(0.5)
    print "testODE.showWarnings() => "
    print testODE.showWarnings()
    print "\ntestODE.indepvariable.depdomain => ", testODE.indepvariable.depdomain
    print "testtraj(0.2, 'aux_other') => ", \
        testtraj(0.2, 'aux_other')

    print "\nNow adding a terminating co-ordinate threshold event..."
    ev_args = {'name': 'threshold',
               'eventtol': 1e-4,
               'eventdelay': 1e-5,
               'starttime': 0,
               'active': True,  # = default
               'term': True,
               'precise': True  # = default
               }
    thresh_ev = Events.makePythonStateZeroCrossEvent('w', 20, 1, ev_args)
    testODE.eventstruct.add(thresh_ev)
    print "Recomputing trajectory:"
    print "traj2 = testODE.compute('test2')"
    traj2 = testODE.compute('test2')
    print "\ntestODE.showWarnings() => "
    print testODE.showWarnings()
    print "\ntestODE.indepvariable.depdomain => ", testODE.indepvariable.depdomain


    print '\n-------- Test: Explicit functional trajectory'
    print "Explicit functional trajectory 'sin_gen' computes sin(t*speed)"

    print "Make 'xdomain' argument smaller than known limits for sine wave: [-1.001, 0.7]"
    args = {'tdomain': [-50, 50],
            'pars': {'speed': 1},
            'xdomain': {'s': [-1., 0.7]},
            'name': 'sine',
            'globalt0': 0.4,
            'pdomain': {'speed': [0, 200]},
            'varspecs': {'s': "sin(globalindepvar(t)*speed)"}}
    sin_gen = ExplicitFnGen(args)
    sintraj = sin_gen.compute('sinewave')
    print "Using alternative form of call (specific names in specific order)..."
    print "sintraj(0.0, ['s'], checklevel=2) = ", sintraj(0.0, ['s'], checklevel=2)
    try:
        print "sintraj(0.8, checklevel=2) \n ", sintraj(0.8, checklevel=2)
    except ValueError, e:
        print "... raised error: ", e
    print "Set limits properly now, to [-1., 1.] ..."
    print "sin_gen.set({'xdomain': {'s': [-1., 1.]}})"
    sin_gen.set(xdomain={'s': [-1., 1.]})
    sintraj2 = sin_gen.compute('sinewave2')
    print "sintraj2(0.8, checklevel=2) = ", sintraj2(0.8, checklevel=2)
