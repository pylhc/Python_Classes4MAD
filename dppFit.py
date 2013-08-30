#--- Chromatic beta fit(nth order), reads linx/y files
#--- dict {bpm: (location, 1/beta*dbeta/dp)}
#--- Jun 1, 2009, R. Calaga

import os, sys
import Numeric as npy
import LinearAlgebra as LA
LAinv=LA.generalized_inverse
from metaclass import twiss


def _polyFit(xData,yData,m):
    #--- p(x)=c[0]+c[1]x+...+c[m]x^m
    a = npy.zeros((m+1,m+1))*1.0
    b = npy.zeros((m+1))*1.0
    s = npy.zeros((2*m+1))*1.0
    for i in range(len(xData)):
        temp = yData[i]
        for j in range(m+1):
            b[j] = b[j] + temp
            temp = temp*xData[i]
        temp = 1.0
        for j in range(2*m+1):
            s[j] = s[j] + temp
            temp = temp*xData[i]
    for i in range(m+1):
        for j in range(m+1):
            a[i,j] = s[i+j]
    return _invert(a,b)

def _invert(a,b,tol=1.0e-9):
    return npy.dot(LAinv(a,tol),b)


def _stDev(c,xData,yData):
    def evalPoly(c,x):
        m=len(c)-1; p=c[m]
        for j in range(m):
            p = p*x + c[m-j-1]
        return p
    n=len(xData)-1; m=len(c)-1; sigma=0.0
    for i in range(n+1):
        p = evalPoly(c,xData[i])
        sigma = sigma + (yData[i] - p)**2
    sigma = npy.sqrt(sigma/(n - m))
    return sigma

def _findFile():
    return filter(lambda p: '.py' not in p,sys.argv)

class dppFit:
    def __init__(self, order=1):
        #--- DBDP= (1/beta)*(dbeta/dp)
        self.order=order;self.DPDB={};nm=_findFile()
        DPP=[];QXY=[];xdat=[];ydat=[];VAR={};SC={};
        print '#- Poly Fit, Order:',self.order
        for k in nm:
            aa=twiss(k);iter=0;
            DPP.append(aa.DPP);
            try: QXY.append(aa.Q1)
            except: QXY.append(aa.Q2)
            for j in aa.NAME: #-- make dict, dpp vs amp
                try: dPP=aa.DPP
                except: dPP=0.0
                if j not in VAR:
                    VAR[j]=[];SC[j]=aa.S[aa.indx[j]]
                VAR[j].append([dPP,aa.AMPX[iter]]);iter+=1

        #--- fit dpp vs amp and norm by amp
        for k in VAR.iteritems():
            k1=npy.array(k[1])
            if len(k1)>2:
                coeff=_polyFit(k1[:,0],k1[:,1],self.order)
                err=_stDev(coeff, k1[:,0], k1[:,1])
                self.DPDB[k[0]]=SC[k[0]],coeff[1]/coeff[0],err
                print SC[k[0]],coeff[1]/coeff[0],err
            #else: print k[0],'not fitted'

        #-- calculate tune & chroms
        if len(DPP)>1:
            coeff=_polyFit(DPP,QXY,self.order)
            #print "Tune:", coeff[0]
            #print "Chrom:", coeff[1],'+/-',\
            #    _stDev(coeff, DPP, QXY)


if __name__ == "__main__":
    #--- usage, reads linx/y files
    dppFit(order=1)

