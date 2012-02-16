from Numeric import *
import sys
import pickle
from string import split
from math import acos

#########################
def list_to_float(arg):
#########################
        vals = []
        for token in arg:
            vals.append(float(token))
        return vals

#########################
class BPM:
#########################
    "a class for BPMs"
    def __init__(self, cutline):
        self.name=cutline[1]
        self.location=float(cutline[2])
        self.data=array(list_to_float(cutline[3:]))


#########################
class rhicdata:
#########################
    "Read the sdds asci RHIC data into an object"

    def __init__(self, filename):
        self.H=[]
        self.V=[]
        self.byname={}
        for line in open(filename):
            if ("#" in line) :
                self.title=line
            else:
                cutline=split(line)
                if (int(cutline[0])==0): #Horizontal BPM
                    self.H.append(BPM(cutline))
                    self.byname[cutline[1]]=self.H[-1].data
                else:                            #Vertical BPM
                   self.V.append(BPM(cutline))
                   self.byname[cutline[1]]=self.H[-1].data


    def correlation(self, istart=20, iend=2000, sigma=0):
        self.phaseadv=[]
        for i in range(1, len(self.H)):
            s1=self.H[i-1].data[istart:iend]
            s2=self.H[i].data[istart:iend]
            s1=s1-average(s1)
            s2=s2-average(s2)
            sum1=sum(s1*s1)/len(s1)
            sum2=sum(s2*s2)/len(s1)
            sum12=sum(s2*s1)/len(s1)
            if sum1*sum2 != 0:
		    self.phaseadv.append(acos(sum12/sqrt((sum1-sigma**2)*(sum2-sigma**2))))
            else:
		    self.phaseadv.append('nan')
        
