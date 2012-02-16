#!/usr/bin/env python

from Numeric import *
from string import split
from LinearAlgebra import *
import sys

################
def factorial(n):
#################
    f = 1
    while (n >= 2):
        f, n = f * n, n - 1
    return f

###################
def phase(mcv):
###################
     c=[]
     for i in range(len(mcv)):
         c.append(math.atan2(mcv[i].imag, mcv[i].real))
     return array(c)

#########################
class twiss:
#########################
    "Twiss parameters from madx output (with free choice of select items)"

    def __init__(self, filename): 
        
        for line in open(filename):
            if ("@ " in line and "%le " in line) :
                label=split(line)[1]
                exec "self."+label+"= "+str(float(split(line)[3]))

            if ("* " in line) :
                alllabels=split(line)
                #print "alllabels",len(alllabels)
                for j in range(1,len(alllabels)):
                    exec "self."+alllabels[j]+"= []"
                            
            if ("$ " in line) :
                alltypes=split(line)
                
            if ("@" not in line and "*" not in line and "$ " not in line) :
                values=split(line)
                for j in range(0,len(values)):
                    if ("%le" in alltypes[j+1]):                      
                      exec "self."+alllabels[j+1]+".append("+str(float(values[j]))+")"
                    if ("s" in alltypes[j+1]):
                      exec "self."+alllabels[j+1]+".append("+values[j]+")"

        for j in range(1,len(alllabels)):
            if ("%le" in alltypes[j]):  
                exec "self."+alllabels[j]+"= array(self."+alllabels[j]+")"           


    def fterms(self):
        self.f3000= []
        self.f2100= []
        self.f1020= []
        self.f1002= []
        self.f1011= []
        for i in range(0,len(self.S)):
            phix = self.MUX-self.MUX[i]
            phiy = self.MUY-self.MUY[i]
            for j in range(0,i):
                phix[j] += self.Q1
                phiy[j] += self.Q2
            dumm=-sum(self.K2L*self.BETX**1.5*e**(3*complex(0,1)*2*pi*phix))/48.
            self.f3000.append(dumm)
            dumm=-sum(self.K2L*self.BETX**1.5*e**(complex(0,1)*2*pi*phix))/16.
            self.f2100.append(dumm)
            dumm=sum(self.K2L*self.BETX**0.5*self.BETY*e**(complex(0,1)*2*pi*(phix+2*phiy)))/16.
            self.f1020.append(dumm)
            dumm=sum(self.K2L*self.BETX**0.5*self.BETY*e**(complex(0,1)*2*pi*(phix-2*phiy)))/16.
            self.f1002.append(dumm)   
            dumm=sum(self.K2L*self.BETX**0.5*self.BETY*e**(complex(0,1)*2*pi*(phix)))/8.
            self.f1011.append(dumm)

        self.f3000=array(self.f3000)
        self.f2100=array(self.f2100)
        self.f1020=array(self.f1020)
        self.f1002=array(self.f1002) 
        self.f1011=array(self.f1011)

        
    def ftermsMatrix(self):
        self.Mft= []
        self.f2100= []
        self.f1020= []
        self.f1002= []
        for i in range(0,len(self.S)):
            phix = self.MUX-self.MUX[i]
            phiy = self.MUY-self.MUY[i]
            for j in range(0,i):
                phix[j] += self.Q1
                phiy[j] += self.Q2
            dumm=-self.BETX**1.5*e**(3*complex(0,1)*2*pi*phix)/48.
            self.Mft.append(dumm.real)
            self.Mft.append(dumm.imag)
            
#            dumm=-sum(self.K2L*self.BETX**1.5*e**(complex(0,1)*2*pi*phix))/16.
#            self.f2100.append(dumm)
#            dumm=sum(self.K2L*self.BETX**0.5*self.BETY*e**(complex(0,1)*2*pi*(phix+2*phiy)))/16.
#            self.f1020.append(dumm)
#            dumm=sum(self.K2L*self.BETX**0.5*self.BETY*e**(complex(0,1)*2*pi*(phix-2*phiy)))/16.
#            self.f1002.append(dumm) 

    def ChromMatrix(self):
        self.Chrom = []
        self.Chrom.append(self.BETX*self.DX/4/pi)
        self.Chrom.append(-self.BETY*self.DX/4/pi)
        self.Chrom=array(self.Chrom)

    def Cmatrix(self):
        self.C = []
        self.gamma = []
        self.f1001 = []
        self.f1010 = []
        J = reshape(array([0,1,-1,0]),(2,2))
        for j in range(0,len(self.S)):
            R = array([[self.R11[j],self.R12[j]],[self.R21[j],self.R22[j]]])
            C = matrixmultiply(-J,matrixmultiply(transpose(R),J))
            C = (1/sqrt(1+determinant(R)))*C

            g11 = 1/sqrt(self.BETX[j])
            g12 = 0
            g21 = self.ALFX[j]/sqrt(self.BETX[j])
            g22 = sqrt(self.BETX[j])
            Ga = reshape(array([g11,g12,g21,g22]),(2,2))

            g11 = 1/sqrt(self.BETY[j])
            g12 = 0
            g21 = self.ALFY[j]/sqrt(self.BETY[j])
            g22 = sqrt(self.BETY[j])
            Gb = reshape(array([g11,g12,g21,g22]),(2,2))
            C = matrixmultiply(Ga, matrixmultiply(C, inverse(Gb)))
            gamma=1-determinant(C)
            self.gamma.append(gamma)
            C = ravel(C)
            self.C.append(C)
            self.f1001.append(((C[0]+C[3])*1j + (C[1]-C[2]))/4/gamma)
            self.f1010.append(((C[0]-C[3])*1j +(-C[1]-C[2]))/4/gamma)
        self.f1001=array(self.f1001)
        self.f1010=array(self.f1010)




# Take care which TWISS or TWISSERR
x=twiss('twissErr')

x.fterms()
f=open()

x.ftermsMatrix()
#print x.Mft

x.ChromMatrix()
#print x.Chrom

#Split matrixes into two to impose the chromaticity
nx=len(x.Mft[1])
F2=transpose(transpose(x.Mft)[nx-2:nx])
F_2=transpose(transpose(x.Mft)[0:nx-2])

Q2=transpose(transpose(x.Chrom)[nx-2:nx])
Q_2=transpose(transpose(x.Chrom)[0:nx-2])

#Target chromaticiy
Q=array([40.44,29.75])
Q2_1=inverse(Q2)
Q2_1Q_2=matrixmultiply(Q2_1,Q_2)
Q2_1Q=matrixmultiply(Q2_1,Q)
K_2=matrixmultiply(generalized_inverse(F_2-matrixmultiply(F2,Q2_1Q_2)), -matrixmultiply(F2,Q2_1Q))

K2=Q2_1Q-matrixmultiply(Q2_1Q_2,K_2)
K=concatenate([K_2,K2])
print len(K)


x.K2L=K

x.fterms()

print abs(x.f3000)

print matrixmultiply(x.Chrom,K)
