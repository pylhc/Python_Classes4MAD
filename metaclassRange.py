#!/usr/bin/env python

from Numeric import *
from string import split
from LinearAlgebra import *
import sys


#########################
class twiss:
#########################
    "Twiss parameters from madx output (with free choice of select items)"

    def __init__(self, filename, srange=None): 
        
        for line in open(filename):
            if ("@ " in line and "%le " in line) :
                label=split(line)[1]
                exec "self."+label+"= "+str(float(split(line)[3]))

            if ("* " in line) :
                alllabels=split(line)
                print alllabels,len(alllabels)
                for j in range(1,len(alllabels)):
                    exec "self."+alllabels[j]+"= []"
                    if (alllabels[j]=='S'):
                        sj=j-1

                    
            if ("$ " in line) :
                alltypes=split(line)

            values=split(line)
            isinrange=1
            if ("@" not in line and "*" not in line and "$ " not in line) :
                if (srange):
                    if (float(values[sj])>srange[0] and float(values[sj])<srange[1]):
                        isinrange=1
                    else:
                        isinrange=0
            
            if ("@" not in line and "*" not in line and "$ " not in line and isinrange==1) :
            
                for j in range(0,len(values)):
                    if ("%hd" in alltypes[j+1]):                      
                      exec "self."+alllabels[j+1]+".append("+str(int(values[j]))+")"
                    
                    if ("%le" in alltypes[j+1]):                      
                      exec "self."+alllabels[j+1]+".append("+str(float(values[j]))+")"
                    if ("s" in alltypes[j+1]):
                      exec "self."+alllabels[j+1]+".append("+values[j]+")"

        for j in range(1,len(alllabels)):
            if (("%le" in alltypes[j]) | ("%hd" in alltypes[j])  ):  
                exec "self."+alllabels[j]+"= array(self."+alllabels[j]+")"           


    def fterms(self):
        self.f3000= []
        self.f2100= []
        self.f1020= []
        self.f1002= []
        self.f20001= []
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
            dumm=sum((self.K1L-2*self.K2L*self.DX)*self.BETX*e**(2*complex(0,1)*2*pi*phix))/16.
            self.f20001.append(dumm)

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

    def beatMatrix(self):
        self.RM = []
        for j in range(0,len(self.S)):
            self.RM.append(-self.BETX*cos(2*pi*(self.Q1-2*abs(self.MUX[j]-self.MUX)))/sin(2*pi*self.Q1))
        for j in range(0,len(self.S)):
            self.RM.append(-self.BETY*cos(2*pi*(self.Q2-2*abs(self.MUY[j]-self.MUY)))/sin(2*pi*self.Q2))
        self.RM=array(self.RM)



# Read the twiss class from the twiss file
#x=twiss('twiss')
# use it as:
#print x.Q1, x.Q2, x.BETX[0]


# run beaMatrix for example:
#x.beatMatrix()
#print x.RM[0]


# BETA-BEAT CORRECTION
# first compute the response matrix by:
# x.beatMatrix()
# Define targetbeat as an array containing the desired changed in Dbeta/beta (x,y)
# targetbeat=Dbeta/beta
# dkl gives the required integrated strengths by: 
# dkl=matrixmultiply(generalized_inverse(x.RM,0.003),targetbeat)


#Want to explore the singular values?:
#svd=singular_value_decomposition(x.RM)


# Computing SEXTUPOLAR RESONANCE TERMS:
# x.fterms()
# The fterms are arrays evaluated at all the elements:
# print x.f3000 , x.f2100  , x.f1020, x.f1002


# COUPLING
# Compute the Cmatrix, gamma, f1001 and f1010 from the Twiss file at all elements
# x.Cmatrix()
# print x.C[0]   (four components of C at the first elements)
# print x.f1001
# ...


