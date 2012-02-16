from PyDSTool import *

qa = Var(['a*3', 'b'], 'array_test')
assert str(qa.eval(a=1)) == '[3,b]'

testq = QuantSpec('d', 'a')
testq.simplify()
assert testq() == 'a'
assert str(testq.eval(a=3)) == '3'
q = QuantSpec('q', 'zeta(yrel(y,initcond(y)),z)-1')
print q.eval({})
assert 'initcond' in str(q.eval({}))
q2=QuantSpec('q','Exp(-spikeTable+b)/k')
assert 'spikeTable' in q2.freeSymbols

#    x = Var('x')
#    print x.isDefined()
#    xs = QuantSpec('x', '1-rel - 2*x + cos(z) + 2e10', 'RHSfuncSpec')
#    x.bindSpec(xs)
#    print x.isDefined(),"\n"
x = Var(QuantSpec('x', '1-rel - 2*x + cos(z) + 2e10', 'RHSfuncSpec'))
p = Par('p')
az = Var(QuantSpec('z', 'myfunc(0,z)+abs(x+1)', 'RHSfuncSpec'))
w = Var('x-1/w[i]', 'w[i,0,1]', specType='RHSfuncSpec')

class myLeaf1(LeafComponent):
    pass
class myLeaf2(LeafComponent):
    pass
class myNode(Component):
    pass

c = myLeaf1('leaf1', targetLangs,
                 compatibleGens=['AGenerator'], compatibleNodes=[myNode])

assert c.isDefined() == False
c.add(x)
print c.freeSymbols, c.isDefined()
c.add(az)
print c.freeSymbols, c.isDefined()
c.add(w)
print c.freeSymbols, c.isDefined()

print "\n"

c.compileFuncSpec()
print c.funcSpecDict

print "\n"

q = Par('qpar')
y = Var(QuantSpec('rel', 'v+p'), domain=[0,1])
g = Fun(QuantSpec('qfunc', '-1.e-05+sin(qpar)*(10.e-5-xtol)'), ['xtol'])
d = myLeaf2('leaf2', targetLangs,
                 compatibleGens=['AGenerator'], compatibleNodes=[myNode])
##    d.add(y)
q_dummy = Var(QuantSpec('q_notpar', '-2+sin(30)'))
g_dummy = Fun(QuantSpec('qfunc_dummy', 'sin(q_notpar)*(10.e-5-xtol)'), ['xtol'])
d.add([q_dummy, g_dummy])  # will delete these later
d.add([q,g])

d2 = myLeaf2('leaf3', targetLangs,
                 compatibleGens=['AGenerator'], compatibleNodes=[myNode])
d2.add([q,g])

v = Var(QuantSpec('v', 'v * myfunc(rel,v) - sin(p)*t', 'RHSfuncSpec'))
# p is a global parameter so this is ok in v
f = Fun(QuantSpec('myfunc', '2.0+s-t+exp(p)'), ['s','t'])
# t is just a local argument here, so it won't clash with its
# occurrence in v (which we'll see is declared as a global
# when we call flattenSpec()).
ipar = Par('ipar')
z = Var('z[i]+v/(i*ipar)', 'z[i,0,5]', specType='RHSfuncSpec')
a = myNode('sys1', targetLangs, [myLeaf1, myLeaf2],
                  compatibleGens=['AGenerator','BGenerator'])
a.add([f,p,y])
print a.isDefined(True)
a.add(c)
print a.freeSymbols, a.isDefined(), a.isComplete()
a.add(d)
print a.freeSymbols, a.isDefined(), a.isComplete()
a.add(d2)
print a.freeSymbols, a.isDefined(), a.isComplete()
a.add(v)
print "Added v"
print a.freeSymbols, a.isDefined(), a.isComplete()
print "Removed v"
a.remove(v)
print a.freeSymbols, a.isDefined(), a.isComplete()
a.add([z,ipar])
print a.freeSymbols, a.isDefined(), a.isComplete()
print "\na._registry -->  "
print a._registry
print "Re-added v"
a.add(v)
print a.freeSymbols, a.isDefined(), a.isComplete()
print "\nv in a -->", v in a

print "\n"
try:
    a.compileFuncSpec()
    test = True
    print "\nWas able to add a non-parameter as a free name to an"
    print "auxiliary function definition! This should not have happened!"
except TypeError, errinfo:
    print "\nSuccessfully could not add a non-parameter as a free name"
    print "in an auxiliary function"
    print "Error was: ", errinfo
    test = False
assert not test
a.remove(['leaf2.qfunc_dummy', 'leaf2.q_notpar'])

print "---------  sys1: funcSpecDict ---------------------"
a.compileFuncSpec()
info(a.funcSpecDict)

print "\n\n-------------  Flatten spec with unravelling\n"
print "\n\ninfo(a.flattenSpec()) --> \n"
info(a.flattenSpec(globalRefs=['t']), "Model specification")
print "\n\n-------------  Flatten spec with no unravelling\n"
print "\n\ninfo(a.flattenSpec(False, globalRefs=['t'])) --> \n"
info(a.flattenSpec(False, globalRefs=['t']), "Model specification")

print "\n\nDemos for functions (results are strings):\n"
h = f(p, -x)
z = QuantSpec('zero','0')
print "h = f(p, -x) --> ", h
print "z = QuantSpec('zero','0') --> ", z
print "f(g(3)*1,h) --> ", f(g(3)*1,h)
print "f(g(p),h) --> ", f(g(p),h)
print "f(g(p),0*h) --> ", f(g(p),0*h)
print "f(g(x),h+z) --> ", f(g(x),h+z)
# e is the math constant, but it doesn't evaluate to a float!
print "f(g(x()),(e+h)/2) --> ", f(g(x()),(e+h)/2)
print "f(g(x()),-h) --> ", f(g(x()),-h)
print "f(g(x()),.5-h+0) --> ", f(g(x()),.5-h+0)
print "Sin(pi+q) --> ", Sin(pi+q)
qsin=QuantSpec('qsin','zv-sin(beta)')
assert str(qsin.eval()) == 'zv-sin(beta)'

print "\n\nDemos for local scope evaluation and **:\n"
print "q=Var('xv+1','qv')"
print "x=Var('3','xv')"
q=Var('xv+1','qv')
x=Var('3','xv')
sc1 = str(q.eval()) == '4'
print "q.eval() == 4? ", sc1
assert sc1
print "a=x/q"
a=x/q
sc2 = str(a) == 'xv/qv'
print "a == xv/qv? ", sc2
assert sc2
sc3 = str(a.eval())=='0.75'
print "a.eval() == 0.75? ", sc3
assert sc3
sc4 = str(a.eval(xv=5))=='5/qv'
print "a.eval(xv=5) == 5/q? ", sc4
assert sc4
sc5 = str(a.eval(xv=5,qv=q()))=='0.83333333333333337'
print "a.eval(xv=5,qv=q()) == 0.83333333333333337? ", sc5
assert sc5
sc6 =str(a.eval({'xv': 10, 'qv': q()}))=='0.90909090909090906'
print "a.eval({'xv': 10, 'qv': q()}) == 0.90909090909090906? ", sc6
assert sc6

print "qs=QuantSpec('qsv','xsv+1')"
print "xs=QuantSpec('xsv','3')"
qs=QuantSpec('qsv','xsv+1')
xs=QuantSpec('xsv','3')
qse = qs.eval()
qt1 = str(qse) == '4'
print "qs.eval() == 4? ", qt1
assert qt1
assert qse.tonumeric() == 4
print "as = xs/qs"
as=xs/qs
qt2 = str(as) == '3/(xsv+1)'
print "as == 3/(xsv+1)? ", qt2
assert qt2
qt3 = str(as.eval()) == '0.75'
print "as.eval() == 0.75? ", qt3
assert qt3
ps = as**xs
print "ps = as**xs"
qt4 = str(ps) == 'Pow(3/(xsv+1),3)'
print "ps == Pow(3/(xsv+1),3)? ", qt4
assert qt4
qt5 = str(ps.eval()) == str(0.75**3)
print "ps.eval() == 0.421875? ", qt5
assert qt5

print "sq=QuantSpec('sv','sin(xsv)')"
print "s2q=QuantSpec('s2v','Sin(xv)')"
sq=QuantSpec('sv','sin(xsv)')
s2q=QuantSpec('s2v','Sin(xv)')
print "sq.eval() --> ", sq.eval()
print "s2q.eval() --> ", s2q.eval()
assert sq.eval().tonumeric() == s2q.eval().tonumeric()
assert sq[:] == ['sin','(','xsv',')']

print "\n\nDemos for multiple quantity definitions:\n"
mp=QuantSpec('p','a + 3*z[4*i-2]') 
m=Var(mp, 'z[i,2,5]', specType='RHSfuncSpec')
v=Var('3*z[i-1]+z4-i', 'z[i,1,5]', specType='RHSfuncSpec')
print "mp=QuantSpec('p','a + 3*z[4*i-2]')"
print "m=Var(mp, 'z[i,2,5]', specType='RHSfuncSpec')"
print "v=Var('3*z[i-1]+z4-i', 'z[i,1,5]', specType='RHSfuncSpec')"
print "v[3] -->", v[3]
assert str(v[3])=='z3'
print "v.freeSymbols -->", v.freeSymbols
assert v.freeSymbols == ['z0']
print "\nModelSpec a already contains 'z0', which was defined as part of"
print "a multiple quantity definition, so check that attempting to add"
print "v to a results in an error ..."
try:
    a.add(v)
    print "v was successfully added -- test failed"
except:
    print "v was rejected - test passed"
print "\nTest of eval method, e.g. on a function f(s,t)..."
print "f.eval(s='1', t='t_val') -->", f.eval(s='1', t='t_val')
print "f.eval(s=1, t='t_val', p=0.5) -->", f.eval(s=1, t='t_val', p=0.5)
print "\nTesting convertPowers():"
cp_tests = ["phi1dot^m3", "1+phi1dot^m3*s",
            "phi1dot**m3", "1+phi1dot**m3*s",
            "sin(x^3)**4", "(2/3)^2.5", "3^cos(x)-pi",
            "3^(cos(x)-pi)", "2^(sin(y**p))"]
for spec in cp_tests:
    print spec, " --> ", convertPowers(spec)


qc=QuantSpec('t', "a+coot+b/'coot'")
assert str(qc.eval()) == 'a+coot+b/"coot"'
coot=QuantSpec('coot', "1.05")
assert str(qc.eval()) == 'a+1.05+b/"coot"'

print "\nTest of function calling with argument names that clash with"
print "bound names inside the function."
x0=Var('x0')
x1=Var('x1')
x2=Var('x2')
F=Fun([x0*x2,x0*5,x2**0.5], [x0,x1,x2], 'F')
print "F=Fun([x0*x2,x0*5,x2**0.5], [x0,x1,x2], 'F')"
print "F(3,2,Sin(x0))) = [3*Sin(x0),15,Pow(Sin(x0),0.5)] ..."
print "  ... even though x0 is a bound name inside definition of F"
assert str(F(3,2,Sin(x0)))=='[3*Sin(x0),15,Pow(Sin(x0),0.5)]'
