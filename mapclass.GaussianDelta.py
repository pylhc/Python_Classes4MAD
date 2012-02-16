#!/usr/bin/env pythonafs

from Numeric import *
from math import *
from string import split
from LinearAlgebra import *
import sys

#################
def gammln( xx):
#################
    cof=[76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5]
    y=x=xx
    tmp=x+5.5
    tmp -= (x+0.5)*log(tmp)
    ser=1.000000000190015;
    for c in cof:
        y=y+1
        ser += c/y        
    return -tmp+log(2.5066282746310005*ser/x)

#########################
class Map:
#########################
    "MAP coefficients from madx-PTC output (fort.18)"

    def __init__(self, order=100, filename='fort.18'): 
        ietall=0
        self.order=order
        xyzd=['x', 'px', 'y', 'py', 'd','s' ]
        for line in open(filename):
            if ("etall" in line) :
                ietall=ietall+1
                exec "self."+xyzd[ietall-1]+"=[]"
            sline=split(line)
            if (len(sline)==8):
                a=[float(sline[1]), int(sline[3]), int(sline[4]), int(sline[5]), int(sline[6]), int(sline[7]) ]
                if ((a[1]+a[2]+a[3]+a[4]+a[5]) <= self.order ):
                    exec "self."+xyzd[ietall-1]+".append("+str(a)+")"
        print "Initialized map with # of coefficients in x,px,y,py:",len(self.x), len(self.px), len(self.y), len(self.py)


    def xf(self, vx,vpx,vy,vpy,vd):
        suma=0
        for  coeff in self.x:
            suma += coeff[0]*vx**coeff[1]*vpx**coeff[2]*vy**coeff[3]*vpy**coeff[4]*vd**coeff[5]
        return suma

    def pxf(self, vx,vpx,vy,vpy,vd):
        suma=0
        for  coeff in self.px:
            suma += coeff[0]*vx**coeff[1]*vpx**coeff[2]*vy**coeff[3]*vpy**coeff[4]*vd**coeff[5]
        return suma

    def yf(self, vx,vpx,vy,vpy,vd):
        suma=0
        for  coeff in self.y:
            suma += coeff[0]*vx**coeff[1]*vpx**coeff[2]*vy**coeff[3]*vpy**coeff[4]*vd**coeff[5]
        return suma

    def pyf(self, vx,vpx,vy,vpy,vd):
        suma=0
        for  coeff in self.py:
            suma += coeff[0]*vx**coeff[1]*vpx**coeff[2]*vy**coeff[3]*vpy**coeff[4]*vd**coeff[5]
        return suma

    def f(self, i):
        return [self.xf(i[0],i[1],i[2],i[3],i[4]),  self.pxf(i[0],i[1],i[2],i[3],i[4]), self.yf(i[0],i[1],i[2],i[3],i[4]),  self.pyf(i[0],i[1],i[2],i[3],i[4]) ]


    def offset(self, xory, i):
        sx=0
        if 'x' in xory:
            mapxory=self.x
        else:
            mapxory=self.y
        for coeff1 in mapxory:
            jj=coeff1[1]
            kk=coeff1[2] 
            ll=coeff1[3] 
            mm=coeff1[4] 
            nn=coeff1[5]
            if ((jj/2==jj/2.) & (kk/2==kk/2.) & (ll/2==ll/2.) & (mm/2==mm/2.) & (nn/2==nn/2.)):
                sigmaprod=pow(i[0], jj)*pow(i[1], kk)*pow(i[2], ll)*pow(i[3], mm)*pow(i[4], nn)
                if (sigmaprod >0):
                    Gammasumln=gammln(0.5+jj/2.)+gammln(0.5+kk/2.)+gammln(0.5+ll/2.)+gammln(0.5+mm/2.)+gammln(0.5+nn/2.)
                    factor=pow(2, (jj+kk+ll+mm+nn)/2.)/pow(pi, 2.5)
                    sx += coeff1[0]*factor*exp(Gammasumln)*sigmaprod  
        return sx


    def sigma(self, xory, i):
        sx2=0
        exec 'mapxory=self.'+xory
        for coeff1 in mapxory:
            for coeff2 in mapxory:
                jj=coeff1[1] + coeff2[1]
                kk=coeff1[2] + coeff2[2]
                ll=coeff1[3] + coeff2[3]
                mm=coeff1[4] + coeff2[4]
                nn=coeff1[5] + coeff2[5]
                if (coeff1 >= coeff2):
                    countfactor=2.0
                    if (coeff1 == coeff2):
                        countfactor=1.0      
                    if ((jj/2==jj/2.) & (kk/2==kk/2.) & (ll/2==ll/2.) & (mm/2==mm/2.) & (nn/2==nn/2.)):
                        sigmaprod=pow(i[0], jj)*pow(i[1], kk)*pow(i[2], ll)*pow(i[3], mm)*pow(i[4], nn)
                        if (sigmaprod >0):
                            Gammasumln=gammln(0.5+jj/2.)+gammln(0.5+kk/2.)+gammln(0.5+ll/2.)+gammln(0.5+mm/2.)+gammln(0.5+nn/2.)
                            factor= countfactor*pow(2, (jj+kk+ll+mm+nn)/2.)/pow(pi, 2.5)
                            sxt = coeff1[0]*coeff2[0]*factor*exp(Gammasumln)*sigmaprod   
                            #print jj,kk,ll,mm,nn,":",coeff1[1],coeff1[2],coeff1[3],coeff1[4],coeff1[5], ":" ,sxt
                            sx2 += sxt             
        return sx2

    def generatelist(self,xory,i):
        sx2=0
        if 'x' in xory:
            mapxory=self.x
        else:
            mapxory=self.y
        exec "self.list"+xory+"=[]"
        for coeff1 in mapxory:
            for coeff2 in mapxory:
                jj=coeff1[1] + coeff2[1]
                kk=coeff1[2] + coeff2[2]
                ll=coeff1[3] + coeff2[3]
                mm=coeff1[4] + coeff2[4]
                nn=coeff1[5] + coeff2[5]
                if (coeff1 >= coeff2):
                    countfactor=2.0
                    if (coeff1 == coeff2):
                        countfactor=1.0                    
                    if ((jj/2==jj/2.) & (kk/2==kk/2.) & (ll/2==ll/2.) & (mm/2==mm/2.) & (nn/2==nn/2.)):
                        sigmaprod=pow(i[0], jj)*pow(i[1], kk)*pow(i[2], ll)*pow(i[3], mm)*pow(i[4], nn)
                        if (sigmaprod >0):
                            Gammasumln=gammln(0.5+jj/2.)+gammln(0.5+kk/2.)+gammln(0.5+ll/2.)+gammln(0.5+mm/2.)+gammln(0.5+nn/2.)
                            factor=countfactor*pow(2, (jj+kk+ll+mm+nn)/2.)/pow(pi, 2.5)
                            sxt = coeff1[0]*coeff2[0]*factor*exp(Gammasumln)*sigmaprod
                            elist=[-abs(sxt),sxt,coeff1[1], coeff1[2],coeff1[3],coeff1[4], coeff1[5], coeff2[1], coeff2[2],coeff2[3],coeff2[4], coeff2[5]]
                            exec "self.list"+xory+".append(elist)"
        exec "self.list"+xory+".sort()"
        return


sigmaFFSstart=[3.9e-6,5.75e-8, 3.76e-7, 1.773e-8,0.01]

#map=Map(4)
#map.generatelist('x',sigmaFFSstart)
#map.generatelist('y',sigmaFFSstart)
#print map.listx[0:9]
#print
#print map.listy[0:20]

#i=[3.9e-6,5.75e-8,3.76e-7,1.773e-8,0]
#print map.offset('x',i), sqrt(map.sigma('x',i)-map.offset('x',i)**2)

#map=Map(1)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
#map=Map(2)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
#map=Map(3)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
#map=Map(4)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
#map=Map(5)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
#map=Map(6)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
#map=Map(7)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
#map=Map(8)
#print sqrt(map.sigma('x',sigmaFFSstart)), sqrt(map.sigma('y',sigmaFFSstart))
