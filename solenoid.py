from numpy import *
from scipy.integrate import quad
from math import sqrt

def Kintegrand(m,theta):
    return 1./sqrt(1.-m*sin(theta)**2)

def K(m):
    Kint= lambda theta: Kintegrand(m,theta)
    return quad(Kint, 0., pi/2)[0]


def Eintegrand(m,theta):
    return sqrt(1.-m*sin(theta)**2)

def E(m):
    Eint= lambda theta: Eintegrand(m,theta)
    return quad(Eint, 0., pi/2)[0]


def Pintegrand(n,m,theta):
    return 1./sqrt(1.-m*sin(theta)**2)/(1.-n*sin(theta)**2)

def P(n,m):
    Pint= lambda theta: Pintegrand(n,m,theta)
    return quad(Pint, 0., pi/2)[0]



def h2(a,rho):
    return 4*a*rho/(a+rho)**2

def k2(a,rho,x):
    return 4*a*rho/((a+rho)**2+x**2)

def bracketBrho(a,rho,x):
    kk2=k2(a,rho,x)
    kk=sqrt(kk2)
    return (kk2-2.)/kk*K(kk2)+2./kk*E(kk2)

def bracketBz(a,rho,x):
    kk2=k2(a,rho,x)
    kk=sqrt(kk2)
    hh2=h2(a,rho)
    return x*kk*(K(kk2)+(a-rho)/(a+rho)*P(hh2,kk2))


def Brho(a,rho,z,muIover4pi,L):
    xplus=z+L/2.
    xminus=z-L/2.
    return muIover4pi/L*sqrt(a/rho)*(bracketBrho(a,rho,xplus)-bracketBrho(a,rho,xminus))
       
         
def Bz(a,rho,z,muIover4pi,L):
    xplus=z+L/2.
    xminus=z-L/2.
    return -muIover4pi/L/2./sqrt(a*rho)*(bracketBz(a,rho,xplus)-bracketBz(a,rho,xminus))

def B(rho,z):
    a=6.3/2.
    muIover4pi=-4./0.4488815
    L=12.5
    return Brho(a,rho, z, muIover4pi, L),  Bz(a,rho, z, muIover4pi, L) 

for i in range(0,100):
    print i*0.001, i*0.1, B(0.001,i*0.1)
