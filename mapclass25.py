'''
.. module:: mapclass25

MAPCLASS is conceived to optimize
the non-linear aberrations  of the
Final Focus System of CLIC.

Written February 2006

.. moduleauthor:: Rogelio Tomas <Rogelio.Tomas@cern.ch>

'''
#!/usr/bin/env python2.5

from numpy import *
#from math import *
from string import split
#from LinearAlgebra import *
import sys

#################
def gammln( x):
#################
    cof=[76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5]
    y=x
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
    '''
    MAP coefficients from madx-PTC output

    :param int order: Calculate map up to this order
    :param string filename: Input filename
    '''

    def __init__(self, order=6, filename='fort.18'): 
        ietall=0
        self.order=order
        xyzd=['x', 'px', 'y', 'py', 'd','s' ]
        for line in open(filename):
            if ("etall" in line) :
                ietall=ietall+1
                setattr(self,xyzd[ietall-1],[])
                setattr(self,xyzd[ietall-1]+'r',1.0*zeros([order+1]*5))
            sline=split(line)
            if (len(sline)==8):
                a=[float(sline[1]), int(sline[3]), int(sline[4]), int(sline[5]), int(sline[6]), int(sline[7]) ]
                if ((a[1]+a[2]+a[3]+a[4]+a[5]) <= self.order ):
                    getattr(self,xyzd[ietall-1]).append(a)
                    getattr(self,xyzd[ietall-1]+"r")[a[1:]] = float(sline[1])
        print "Initialized map with # of coefficients in x,px,y,py:",len(self.x), len(self.px), len(self.y), len(self.py)

    def comp(self, map2):
        if (len(self.x) < len(map2.x)):
            print "Self map has fewer elements than map2!!"
            print "This gives a wrong result"
        chi2=0
        for v in self.x:
            if v[5]==0: chi2+=(v[0]-map2.xr[v[1],v[2],v[3],v[4],v[5]])**2
        for v in self.y:
            if v[5]==0: chi2+=(v[0]-map2.yr[v[1],v[2],v[3],v[4],v[5]])**2
        for v in self.py:
            if v[5]==0: chi2+=(v[0]-map2.pyr[v[1],v[2],v[3],v[4],v[5]])**2
        for v in self.px:
            if v[5]==0: chi2+=(v[0]-map2.pxr[v[1],v[2],v[3],v[4],v[5]])**2
        return chi2



    def compc(self, map2):
        if (len(self.x) < len(map2.x)):
            print "Self map has fewer elements than map2!!"
            print "This gives a wrong result"
        chi2=0
        for v in self.x:
            chi2+=(v[0]-map2.xr[v[1],v[2],v[3],v[4],v[5]])**2
        for v in self.y:
            chi2+=(v[0]-map2.yr[v[1],v[2],v[3],v[4],v[5]])**2
        for v in self.py:
            chi2+=(v[0]-map2.pyr[v[1],v[2],v[3],v[4],v[5]])**2
        for v in self.px:
            chi2+=(v[0]-map2.pxr[v[1],v[2],v[3],v[4],v[5]])**2
        return chi2




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
        '''
        Calculate the beam offset

        :param string xory: Which coordinate to calculate for (x,y,px, or py)
        :param list i: Size of beam in sigma [x,px,y,py]
        '''
        sx=0
        mapxory=getattr(self,xory)
        for coeff1 in mapxory:
            jj=coeff1[1]
            kk=coeff1[2] 
            ll=coeff1[3] 
            mm=coeff1[4] 
            nn=coeff1[5]
            if ((jj/2==jj/2.) & (kk/2==kk/2.) & (ll/2==ll/2.) & (mm/2==mm/2.) & (nn/2==nn/2.)):
                sigmaprod=pow(i[0], jj)*pow(i[1], kk)*pow(i[2], ll)*pow(i[3], mm)*pow(i[4]/2., nn)
                if (sigmaprod >0):
                    Gammasumln=gammln(0.5+jj/2.)+gammln(0.5+kk/2.)+gammln(0.5+ll/2.)+gammln(0.5+mm/2.)
                    factor=pow(2, (jj+kk+ll+mm)/2.)/pow(pi, 2.)/(nn+1)
                    sx += coeff1[0]*factor*exp(Gammasumln)*sigmaprod  
        return sx


    def sigma(self, xory, i):
        '''
        Calculate the beam size in sigma.

        :param string xory: Which coordinate to calculate for (x,y,px, or py)
        :param list i: Size of beam in sigma [x,px,y,py]
        '''
        sx2=0
        mapxory=getattr(self,xory)
        for coeff1 in mapxory:
            for coeff2 in mapxory:
                if (coeff1[1:] >= coeff2[1:]):
                    countfactor=2.0
                    if (coeff1[1:] == coeff2[1:]):
                        countfactor=1.0
                    jj=coeff1[1] + coeff2[1]
                    kk=coeff1[2] + coeff2[2]
                    ll=coeff1[3] + coeff2[3]
                    mm=coeff1[4] + coeff2[4]
                    nn=coeff1[5] + coeff2[5]
                    if ((jj/2==jj/2.) & (kk/2==kk/2.) & (ll/2==ll/2.) & (mm/2==mm/2.) & (nn/2==nn/2.)):
                        sigmaprod=pow(i[0], jj)*pow(i[1], kk)*pow(i[2], ll)*pow(i[3], mm)*pow(i[4]/2., nn)
                        if (sigmaprod >0):
                            Gammasumln=gammln(0.5+jj/2.)+gammln(0.5+kk/2.)+gammln(0.5+ll/2.)+gammln(0.5+mm/2.)
                            factor= countfactor*pow(2, (jj+kk+ll+mm)/2.)/pow(pi, 2.)/(nn+1)
                            sxt = coeff1[0]*coeff2[0]*factor*exp(Gammasumln)*sigmaprod   
                            #print jj,kk,ll,mm,nn,":",coeff1[1],coeff1[2],coeff1[3],coeff1[4],coeff1[5], ":" ,sxt
                            sx2 += sxt             
        return sx2

    def generatelist(self,xory,i):
        sx2=0
        mapxory=getattr(self,xory)
        lxory='list'+xory
        setattr(self,lxory,[])
        for coeff1 in mapxory:
            for coeff2 in mapxory:
                if (coeff1 >= coeff2):
                    jj=coeff1[1] + coeff2[1]
                    kk=coeff1[2] + coeff2[2]
                    ll=coeff1[3] + coeff2[3]
                    mm=coeff1[4] + coeff2[4]
                    nn=coeff1[5] + coeff2[5]
                    if ((jj%2==0) & (kk%2==0) & (ll%2==0) & (mm%2==0) & (nn%2==0)):
                        sigmaprod=pow(i[0], jj)*pow(i[1], kk)*pow(i[2], ll)*pow(i[3], mm)*pow(i[4]/2., nn)
                        if (sigmaprod >0):
                            Gammasumln=sum([gammln(0.5+jj/2.),gammln(0.5+kk/2.),gammln(0.5+ll/2.),gammln(0.5+mm/2.)])
                            factor=pow(2, (jj+kk+ll+mm)/2.)/pow(pi, 2.)/(nn+1)
                            if coeff1 > coeff2:
                                factor*=2
                            sxt = coeff1[0]*coeff2[0]*factor*exp(Gammasumln)*sigmaprod
                            elist=[-abs(sxt),sxt,coeff1[1], coeff1[2],coeff1[3],coeff1[4], coeff1[5], coeff2[1], coeff2[2],coeff2[3],coeff2[4], coeff2[5]]
                            getattr(self,lxory).append(elist)
        getattr(self,lxory).sort()
        return


############################################
#### some examples of usage #####################
############################################
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
