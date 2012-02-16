#! /usr/local/bin/python
"""Tests on symbolic differentiation using ModelSpec Quantity objects.

    Robert Clewley, September 2005.
"""

from PyDSTool import *

print "Showing the variety of ways that symbolic Diff() can be used:"

f1 = '[-3*x**2+2*(x+y),-y/2]'
f2 = ['-3*x**2+2*(x+y)','-y/2']
x=Var('x')
y=Var('y')
f3 = [-3*Pow(x,2)+2*(x+y),-y/2] 
f4 = ['-3*x**2.+2*(x+y)','-y/2.']
xx = QuantSpec('dummy','x')

print "The following strings should all be identical"
print Diff(f1,'x')
print Diff(f2,'x')
print Diff(f3,'x')
print Diff(f1,x)
print Diff(f2,x)
print Diff(f3,x)
print Diff(f3,xx)
print "\n"
print Diff(f1, ['x','y'])
print Diff(f1, [x,y])
print Diff(f1, [xx,y])
print Diff(f2, ['x','y'])
print Diff(f2, [x,y])
print Diff(f3, ['x','y'])
print Diff(f3, [x,y])

print "--------------------------\n"

print "Now some more complex tests..."
t=Var('t')
s=Var('s')

assert str(Diff(Pow((t*5),2),t)) != '0'

p=Par('3.','p')
f = Fun(QuantSpec('f', str(2.0+s-10*(t**2)+Exp(p))), ['s','t'])
f_0=Fun(QuantSpec('f_0', str(Diff(f(s,t),t))), ['t'])
print 2*f.eval(s=3,t=t)
print Diff('-10*Pow(t,2)','t')
print Diff(2*f.eval(s=3,t=t), t)
print Diff(3+t*f.eval(s=3,t=t),t)
print Diff(3+t*f(s,t),t).eval(s=3,t=1,p=p)
print Diff(3+t*f(s,t),t).eval(s=3,t=1,p=p())
assert Diff(str(f(s,t)),'t') == Diff(f(s,t),t)
q1=Diff(f(s,t),t)
q2=Diff(str(f(s,t)),t)
assert q1 == q2
q1.difference(q2)
print "\n"
print Diff(f(t,s),t)
print Diff(2*f(3,t*5), t)
assert str(Diff(2*f(3,t*5), t)) != str(0)
assert f(s,t) != f(t,s)

print f(s,t).eval()
q=f(s,t)
print q.eval()

print Diff('g(s)',s)
print Diff('g(s)',s).eval()
dg_dt=Fun(QuantSpec('g_0', '2-Sin(t/2)'),['t'])
assert str(Diff('g(t)',t).eval()) != 'g_0(t)'
print "\n\n"
print Diff('g(s)',s)
print Diff('g(s)',s).eval()


g=Fun('',[t],'g') # declare empty function
assert str(g(t)) == 'g(t)'
print Diff(g(s),s).eval()

assert eval(str(Diff('pow(1,2)*t','t'))) == 1
assert eval(str(Diff(Pow(1,2)*t,t))) == 1
assert str(Diff(Sin(Pow(t,1)),t)) == 'Cos(t)' 

q=QuantSpec('q','-0+3+pow(g(x)*h(y,x),1)*1')
print Diff(q,'x')
assert str(Diff(q,'x')) == 'g_0(x)*h(y,x)+g(x)*h_1(y,x)'
print Diff(q,'x').eval()
assert str(Diff(q,'x').eval()) == '(2-Sin(x/2))*h(y,x)+g(x)*h_1(y,x)'

p0=Var('p0')
p1=Var('p1')

pv=Var([p0,p1], 'p')
print pv()
print pv.eval()

u=Var('Pi/(2*Sin(Pi*t/2))','u')
assert u.eval(t=1).tonumeric() == pi/2

#------------------
print "\nSymbolic vector tests"

q0=Var(p0+3,'q0')
q1=Var(Diff(1+Sin(Pow(p0,3)+q0),p0),'q1')

qv=Var([q0,q1], 'q')
print qv()
print qv.eval()

v=Var('v')
w=Var('w')
f=Var([-3*Pow((2*v+1),3)+2*(w+v),-w/2], 'f')
    
df = Diff(f, [v,w])
print df
dfe = df.eval(v=3,w=10).tonumeric()
print dfe
assert isinstance(dfe, Array)
assert isinstance(df.fromvector(), list)

y0=Var('y0')
y1=Var('y1')
y2=Var('y2')
t=Var('t')

ydot0=Fun(-0.04*y0 + 1e4*y1*y2, [y0, y1, y2], 'ydot0')
ydot2=Fun(3e7*y1*y1, [y0, y1, y2], 'ydot2')
ydot1=Fun(-ydot0(y0,y1,y2)-ydot2(y0,y1,y2), [y0, y1, y2], 'ydot1')

F = Fun([ydot0(y0,y1,y2),ydot1(y0,y1,y2),ydot2(y0,y1,y2)], [y0,y1,y2], 'F')
assert F.dim == 3
assert str(Diff(F,[y0,y1,y2])) == '[[-0.04,10000*y2,10000*y1],[0.040000000000000001,(-10000*y2)-30000000*2*y1,-10000*y1],[0,60000000*y1,0]]'

jac=Fun(Diff(F,[y0,y1,y2]), [t, y0, y1, y2], 'Jacobian')
assert jac(t, 0.1,y0+1,0.5).eval(y0=0) == jac(t, 0.1,1+y0,0.5).eval(y0=0)
assert jac(t, 0.1,y0,0.5) == jac(t, 0.1,0+y0,0.5)

f1 = Fun([-3*x**3+2*(x+y),-y/2], [x,y], 'f1')
f2 = ['-3*x**3+2*(x+y)','-y/2']
f3 = [-3*x**3.+2*(x+y),-y/2.]
print "\n\nVector-valued function f(x,y) =", f1
print "The function string can be passed to Diff in various ways..."
print str(f1)
print str(f2)
print str(f3)
print "\nThe following outputs are for Diff(f,'x') for each of these forms"
print "They should all be the same (except for some may contain decimal points)"
x=Var('x')
y=Var('y')
f4 = [-3*Pow((2*x+1),3)+2*(x+y),-y/2] 
xx = QuantSpec('dummy','x')
f5=Var([-3*Pow((2*x+1),3)+2*(x+y),-y/2], 'f5')

assert Diff(f1,x) == Diff(f1,'x')
print Diff(f1,x)
print Diff(f3,x)
print Diff(f3,xx)
print Diff(f4,x)
print Diff(f4,xx)
print "\nExamples of Jacobian Diff(f, [x,y])..."
assert Diff(f1, [x,y]) == Diff(f1, ['x','y']) == Diff(f1(x,y), [x,y])
print Diff(f2, ['x','y'])
print Diff(f3, ['x','y'])
print Diff(f1, [xx,y])
print Diff(f1, [xx,'y'])
print Diff(f2, [x,y])
print Diff(f3, [x,y]), "\n"
print Diff(f4, [x,y])
df5 = Diff(f5, [x,y])
print df5
print df5.eval(x=3,y=10).tonumeric()
print df5.eval(x=3,y=10).fromvector(0)
print df5.fromvector(0)
assert isinstance(df5.fromvector(), list)
a = df5.fromvector(0).eval(x=3,y=10).tonumeric()
b = df5.eval(x=3,y=10).tonumeric()[0]
assert a[0]==b[0] and a[1]==b[1]
