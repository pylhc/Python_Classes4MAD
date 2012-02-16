#-- class for fitting algorithms 
#-- Polynomial fitting (least squares type) - see usage below
#-- non-linear least sq - NOT IMPLEMENTED
#--- R. Calaga, July 9, 2009

import os, sys
try:
    import numpy as npy; PINV=npy.linalg.pinv; 
    tp=npy.transpose; INV=npy.linalg.inv;
    from metaclass25 import twiss;
except ImportError:
    try:
        import Numeric as npy
        from LinearAlgebra import generalized_inverse as PINV
        from metaclass import twiss; tp=npy.transpose
    except ImportError: print "No Numpy -or- Numeric";sys.exit()

    


class FIT:
    '''
    R = S.R;                % The "R" from A = QR
    d = (R'*R)\eye(n+1);    % The covariance matrix
    d = diag(d)';           % ROW vector of the diagonal elements
    MSE = (S.normr^2)/S.df; % variance of the residuals
    x = inv(A'*inv(V)*A)*A'*inv(V)*B
    mse = B'*(inv(V) - inv(V)*A*inv(A'*inv(V)*A)*A'*inv(V))*B./(m-n)
    S = inv(A'*inv(V)*A)*mse
    stdx = sqrt(diag(S))
    '''
    
    def __init__(self, X,Y,m,tol=1.0e-15):
        if len(X) != len(Y):  raise ValueError, 'unequal length'
        self.xD=X; self.yD=Y; self.m=m; self.tol=tol; 
        self.err=[]; self.dof=len(X)-m-1
    
    def rmatrix(self):
        #--- p(x)=c[0]+c[1]x+...+c[m]x^m
        a = npy.zeros((self.m+1,self.m+1))*1.0
        b = npy.zeros((self.m+1))*1.0
        s = npy.zeros((2*self.m+1))*1.0
        for i in range(len(self.xD)):
            temp=self.yD[i]
            for j in range(self.m+1):
                b[j]+=temp; temp=temp*self.xD[i]
            temp=1.0
            for j in range(2*self.m+1):
                s[j]+=temp; temp=temp*self.xD[i]
        for i in range(self.m+1):
            for j in range(self.m+1):
                a[i,j]=s[i+j]
        return a,b
        
    def help(self):
        print 'a=FIT(xd,yd,2,tol=1e-14)'
        print 'a=linreg()'
        print 'Using SVD: coeff,errors=a.SVDFIT()' 
        print 'Using QR: coeff,errors=a.QRFIT()' 

    def QRFIT(self):
        self.a,self.b=self.rmatrix()
        print "Using PolyFit with QR Decomposition"
        Q,R=npy.linalg.qr(self.a)
        d=INV(npy.dot(INV(R),R))
        d=npy.diag(npy.dot(d,npy.eye(len(R))))
        c=npy.dot(INV(R),npy.dot(INV(Q),self.b))      
        #dc = t*sqrt(1+sumsq (A/s.R, 2))*resid/sqrt(sef.dof)                
        e=self.error(c)
        return c,e

    def SVDFIT(self):
        #--- x = inv(At*A)*(At*b)
        self.a,self.b=self.rmatrix()
        print 'Using PolyFit with SVD, Tol=',self.tol
        c=npy.dot(PINV(self.a,self.tol),self.b)
        e=self.error(c)
        return c,e

    def var(self,x):
        return npy.mean(npy.abs(x-npy.mean(x))**2) 
    
#    def error(self, c):
#        resid=npy.dot(self.a,c)-self.b
#        rms=npy.sqrt(npy.inner(resid,resid)/len(resid))
#        SS=(YtY-npy.dot(tp(c),XtY))/self.dof
#        Rsq=1-SS/(YtY-len(self.b)*npy.mean(self.b)**2)
#        sig=self.residual(c)
#        return err
    
    def error(self,c):
        YtY=npy.dot(tp(self.b),self.b)
        XtX=npy.dot(tp(self.a),self.a)
        XtY=npy.dot(tp(self.a),self.b)
        def evalPoly(c,x):
            m=len(c)-1; p=c[m]
            for j in range(m): p=p*x+c[m-j-1]
            return p    
        n=len(self.xD)-1; m=len(c)-1; sigma=0.0
        for i in range(n+1):
            p = evalPoly(c,self.xD[i])
            sigma+=(self.yD[i]-p)**2
        sig=npy.sqrt(sigma/(n-m))
        Rsq=1-sig/(YtY-len(self.b)*npy.mean(self.b)**2)
        print 'R-Squared:', Rsq
        err=npy.sqrt(npy.diag(sig**2*INV(XtX)))
        return err

    def linreg(self):
        """
        Summary
               Linear regression of y = ax + b
        Usage
               [coeff], [err] = FIT.linreg()
               Returns coefficients to the regression line "y=ax+b" 
               from x[] and y[], and R^2 Value
        """        
        print "Using linear regression"
        N=len(self.xD); Sx=Sy=Sxx=Syy=Sxy=0.0
        for x, y in map(None, self.xD, self.yD):
            Sx = Sx + x; Sy = Sy + y
            Sxx = Sxx + x*x; Syy = Syy + y*y
            Sxy = Sxy + x*y
        det = Sxx*N - Sx*Sx
        a, b = (Sxy*N - Sy*Sx)/det, (Sxx*Sy - Sx*Sxy)/det
        meanerror = residual = 0.0
        for x, y in map(None, self.xD, self.yD):
           meanerror = meanerror + (y - Sy/N)**2
           residual = residual + (y - a * x - b)**2
        RR = 1 - residual/meanerror
        print "R-Squared:", RR
        ss = residual / (N-2)
        Var_a, Var_b = ss * N / det, ss * Sxx / det
        return [a,b],[npy.sqrt(Var_a),npy.sqrt(Var_b)]


if __name__ == "__main__":
    #import numpy as npy
    from random import random; 
    xd=npy.arange(10)
    yd=map(lambda x: 1.2*x+4+(0.5-random())*5e-4, xd)
    
    #--- polyFit(xData, yData, polynomial_order, tolerance)
    a=FIT(xd,yd,1,tol=1e-14)
    cf,ce=a.QRFIT() #--- or QRFIT
    for j in range(len(cf)): 
        print 'x^'+str(j)+': ',cf[j],'+/-',ce[j]

    #-- plot
    #import pylab
    #ydf=map(lambda x: cf[2]*x**2+cf[1]*x+cf[0], xd)
    #pylab.scatter(xd,yd);pylab.plot(xd,ydf);pylab.show()
    
    print '\n'
    cf,ce=a.linreg() 
    print 'x^0:',cf[1], '+/-', ce[1]
    print 'x^1:',cf[0], '+/-', ce[0]
    
    
