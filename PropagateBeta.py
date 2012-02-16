from math import *
from simplex import Simplex

L=22.965000
b=1.
KL=sqrt(0.008730196766) * 6.37
K=sqrt(0.008730196766)

def alf1(b,L,w):
    return (L-w)/b


def bet1(b,L,w):
    return b+(L-w)**2/b

def bet2(b,L,w,KL,K):
    bet1=b+(L-w)**2/b
    alf1=(L-w)/b
    return bet1*cos(KL)**2+alf1*2*cos(KL)*sin(KL)/K+1./b*(sin(KL)/K)**2

def beth2(b,L,w,KL,K):
    bet1=b+(L-w)**2/b
    alf1=(L-w)/b
    return bet1*cosh(KL)**2+alf1*2*cosh(KL)*sinh(KL)/K+1./b*(sinh(KL)/K)**2


def abh(b,L,w,KL,K):
    bet1=b+(L-w)**2/b
    alf1=(L-w)/b
    KL2=2*KL
    sinhc=sinh(KL2)/KL2
    res=0.5*bet1*(1+sinhc)+alf1*sinh(KL)**2/KL/K+(sinhc-1)/(2.*b*K**2)
    return res


def ab(b,L,w,KL,K):
    bet1=b+(L-w)**2/b
    alf1=(L-w)/b
    KL2=2*KL
    sinc=sin(KL2)/KL2
    res=0.5*bet1*(1+sinc)+alf1*sin(KL)**2/KL/K+(1-sinc)/(2.*b*K**2)
    return res


#print bet1(b,L,0), bet2(b,L,0,KL,K), beth2(b,L,0,KL,K)
#print -alf1(b,L,0)
print "abh, ab"
print abh(b,L,0,KL,K) , ab(b,L,0,KL,K)
print abh(b,L,0,KL,K)*6.37*7.5e-06/4/pi , ab(b,L,0,KL,K)*6.37*7.5e-06/4/pi


dbx= 2*(1/tan(2*pi*0.31)*(1-cos(2*pi*0.002322782978))+sin(2*pi*0.002322782978))/(6.37*7.5e-06)

dby= -2*(1/tan(2*pi*0.32)*(1-cos(-2*pi*0.00295657854))+sin(-2*pi*0.00295657854))/(6.37*7.5e-06)

print dbx, dby

import sys
sys.exit()

def chi2(d):
    global L,KL,K,dqx,dqy
    b=d[0]
    w=d[1]
    c2= (ab(b,L,w,KL,K)-dbx)**2+(abh(b,L,-w,KL,K)-dby)**2
    print c2
    return c2

 

    
print "guess:",chi2([1.,0])
print ab(b,L,0,KL,K), dbx
print abh(b,L,0,KL,K),dby

import sys
#sys.exit()

deltafamilies1=[1.7,-0.6]
incr=[0.01,0.001]
s = Simplex(chi2, deltafamilies1, incr)
values, err, iter = s.minimize(epsilon=0.00000001,maxiters=10000)
print
print 'args = ', values
print 'error = ', err
print 'iterations = ', iter

print "guess:",chi2(values)
b=values[0]
w=values[1]
print ab(b,L,w,KL,K), dbx
print abh(b,L,w,KL,K),dby


