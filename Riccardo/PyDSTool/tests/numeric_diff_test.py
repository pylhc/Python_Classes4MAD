#! /usr/local/bin/python
"""Testing numerical differentiation using diff.

   Robert Clewley, September 2005.
"""

from PyDSTool import *


print "******** Test scalar f\n"

def f1(x):
    if type(x) in [int, float]:
        return 3*x*x
    else:
        print type(x)
        raise TypeError

def f2(x):
    assert type(x) in [Array, NArray] and len(x) == 4
    return sqrt(sum((x-array([1.,0.,1.,1.]))**2))

print "f1: R -> R"
df1 = diff(f1, 0.5)
print "df1 = ", df1, " ... should be 3."
assert df1 - 3 < 1e7

print "\n"

print "f2: R^4 -> R"
x0 = array([0., 0., 0., 0.])
print "f2(x0) = ", f2(x0)
df2 = diff(f2, x0)
print "diff(f2, x0) = ", simplifyMatrixRepr(df2)
dx1 = 0.25
dx2 = 0.25
print "Compare the 1st order Taylor expansion of f2 at x0+dx = x0+[dx1,dx2,0,0] in the first argument,",
print simplifyMatrixRepr(f2(x0) + diff(f2, x0, vars=[0])*dx1),
print ", to the value of f2(x0+dx) =",
print f2(array([dx1,dx2,0.,0.])+x0), "\n\n"

print "3-dimensional f3:"
def f3(y):
    assert len(y) == 3
    assert type(y) in [Array, NArray]
    return y**2-array([1.,0.,1.])
y0 = array([3., 2., 1.])
df3 = diff(f3, y0)
print simplifyMatrixRepr(df3)

print "\n"
print "\n******** Test vector f\n"

print "f4 : R^3 -> R^3"
def f4(y):
    assert len(y) == 3
    assert type(y) in [Array, NArray]
    return array([y[0]*y[2], y[0]*5, y[2]*0.5])
df4 = diff(f4, y0, axes=[1])
print "diff(f4, y0, axes=[1]) =", simplifyMatrixRepr(df4)
dy = array([0.2, 0., -.2])
print "\nComparing 1st order Taylor series for nearby point to actual value:"
print "Compare ", simplifyMatrixRepr(f4(y0+dy)), " to ",
print simplifyMatrixRepr(mat(f4(y0)) + diff(f4, y0)*dy)

print "\n\n******** Test Point f\n"

x0=Var('x0')
x1=Var('x1')
x2=Var('x2')
f5_x0 = Fun(x0*x2, [x0,x1,x2], 'f5_x0')
f5_x1 = Fun(x0*5, [x0,x1,x2], 'f5_x1')
f5_x2 = Fun(x2**0.5, [x0,x1,x2], 'f5_x2')
y0pt = Point({'coordarray': y0, 'coordnames': ['x0', 'x1', 'x2']})
y1pt = Point({'coordarray': array([3.1, 2., .94]), 'coordnames': ['x0', 'x1', 'x2']})

# could also have defined F directly from f5_x[i] definitions
F=Fun([f5_x0(x0,x1,x2),f5_x1(x0,x1,x2),f5_x2(x0,x1,x2)],[x0,x1,x2], 'F')

print "f5: R^3 -> R^3 defined as:"
print "[x0, x1, x2] |-> [",f5_x0(x0,x1,x2),',',f5_x1(x0,x1,x2),',',f5_x2(x0,x1,x2),']',"\n"
print "Symbolic differentiation of f5 evaluated at y0, by calling"
print "Diff([f5_x0(x0,x1,x2),f5_x1(x0,x1,x2),f5_x2(x0,x1,x2)],[x0,x1,x2]).eval(y0pt):"
print Diff([f5_x0(x0,x1,x2),f5_x1(x0,x1,x2),f5_x2(x0,x1,x2)],[x0,x1,x2]).eval(y0pt)

def f5(z):
    assert isinstance(z, Point)
    return Point({'coorddict': {'x0': z('x0')*z('x2'),
                                'x1': z('x0')*5,
                                'x2': z('x2')**0.5}})
df5 = diff(f5, y0pt, axes=['x1','x2'])
print "\ndiff(f5, y0pt, axes=['x1','x2']):"
print simplifyMatrixRepr(df5)
print "\nComparing 1st order Taylor series for nearby point to actual value:"
print array(f5(y1pt)), "vs.", 
print array(f5(y0pt)) + simplifyMatrixRepr(diff(f5, y0pt))*(y1pt-y0pt)
