#! /sw/bin/python2.3
#Author:    Drew LaMar
#Date:      3 June 2006
#$Revision: 0.3.0 $
""" EXAMPLE: Lotka-Volterra (3D)

    Drew LaMar, March 2006
"""

from PyDSTool import *

pars = {'b1': 5., 'b2': 12., 'b3': 12.,
        'a11': 5., 'a12': 5., 'a13': 0.4,
        'a21': 6., 'a22': 6., 'a23': 6.,
        'a31': 13., 'a32': 3., 'a33': 3.}

icdict = {'x1': 0.6, 'x2': 0.31304348, 'x3': 1.08695652}

# Set up model
x1str = 'x1*(b1 - a11*x1 - a12*x2 - a13*x3)'
x2str = 'x2*(b2 - a21*x1 - a22*x2 - a23*x3)'
x3str = 'x3*(b3 - a31*x1 - a32*x2 - a33*x3)'

DSargs = args(name='LotkaVolterra')
DSargs.pars = pars
DSargs.varspecs = {'x1': x1str, 'x2': x2str, 'x3': x3str}
DSargs.ics = icdict

testDS = Generator.Vode_ODEsystem(DSargs)

# Set up continuation class
PyCont = ContClass(testDS)

PCargs = args(name='EQ1', type='EP-C')
PCargs.freepars = ['a13']
PCargs.StepSize = 1e-2
PCargs.MaxNumPoints = 50
PCargs.MaxStepSize = 1e-1
PCargs.LocBifPoints = 'all'
PCargs.verbosity = 2
PyCont.newCurve(PCargs)

print 'Computing curve...'
start = clock()
PyCont['EQ1'].forward()
print 'done in %.3f seconds!' % (clock()-start)

PCargs.name = 'HO1'
PCargs.type = 'H-C1'
PCargs.initpoint = 'EQ1:H1'
PCargs.freepars = ['a13','a23']
PCargs.MaxNumPoints = 10
PCargs.MaxStepSize = 0.1
PCargs.LocBifPoints = 'all'
PyCont.newCurve(PCargs)

print 'Computing Hopf curve...'
start = clock()
PyCont['HO1'].backward()
print 'done in %.3f seconds!' % (clock()-start)

# Plot
PyCont.display(('a13','x2'))
