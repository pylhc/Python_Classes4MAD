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

#x=twiss('twiss.4Cmatrix')
#print x.Q1, x.Q2, x.ENERGY, x.LENGTH, x.BETX[1], x.NAME[0], x.MUX[1]
#x.fterms()
#print x.f3000

x.Cmatrix()

def aver(arr) :
    return sum(abs(arr))/len(arr)

def aver2(arr) :
    return  sum(abs(arr*arr))/len(arr)

af1001=aver(x.f1001)
af1010=aver(x.f1010)
a2f1001=aver2(x.f1001)
a2f1010=aver2(x.f1010)

print "avef1001=",af1001,";"
print "avef1010=",af1010,";"
print "emitratio=",8*(af1001**2+af1010**2),";"


f = open('f1xxx.dat','w')
for i in range(len(x.f1001)):
    f.write(str(x.S[i])+' '+\
            str(x.f1001[i].real)+' '+str(x.f1001[i].imag)+' '+\
            str(x.f1010[i].real)+' '+str(x.f1010[i].imag)+' '+\
            str(x.DY[i])+'\n')
f.close()





#sqrt(a2f1001-af1001**2), af1010, sqrt(a2f1010-af1010**2), 8*(af1001**2+af1010**2), 16*(af1001*sqrt(a2f1001-af1001**2)+af1010*sqrt(a2f1010-af1010**2))
