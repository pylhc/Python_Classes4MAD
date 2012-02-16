#!/usr/bin/env python

import sys, gzip
try:
    from numpy import *
except ImportError:
    try:
        from Numeric import *
        from LinearAlgebra import *
    except ImportError: sys.exit()

from string import split
from string import replace

def replac(nm): #-- added aug12_08 calaga
    nm=replace(nm,"-","_")
    nm=replace(nm,".","_")
    return nm


#########################
class YASPtwiss:
#########################
    "Read orbit data from YASP output in the same way as metaclass"
    def __init__(self, filename=''): 


        self.HNAME=[]
        self.HCO=[]
        self.Hindx={}
        self.VNAME=[]
        self.VCO=[]
        self.Vindx={}
        self.HcorrNAME=[]
        self.Hcorr=[]
        self.VcorrNAME=[]
        self.Vcorr=[]

        if filename=='':
            return

        
        if '.gz' in filename:
            f=gzip.open(filename, 'rb')
        else:
            f=open(filename, 'r')



        #fspl=filename.split('/')[-1].split('.')
        #self.exciter=fspl[1]+fspl[2]
        #self.plane=fspl[3]
        #if "minus" in fspl[4]:
        #  self.str=-1.0
        #else:
        #  self.str=1.0
        


        corr=0
        count=0
        for line in f:
            cl=line.split()
            if ("@ " in line) :
                label=replac(cl[1])
		try:
                	exec "self."+label+"= '"+str((cl[3]))+"'"
		except:
			print "Problem parsing:", line
			#print "Ignored"

            if "# CORRECTOR" in line:
                corr=1

            #if ("* " in line) :
            #    alllabels=split(line)
            #    print "alllabels",len(alllabels)
            #    for j in range(1,len(alllabels)):
            #        exec "self."+alllabels[j]+"= []"



            
            if corr==0 and cl[1]=='H' and cl[8]=="OK":
                
                self.HCO.append(float(cl[3])*1e-3)
                #self.Hindx[replace(cl[0], '"', '')]=len(self.HNAME)-1
                #self.Hindx[replace(cl[0], '"', '').upper()]=len(self.HNAME)-1
                
                self.HNAME.append(cl[0])
                self.Hindx[self.HNAME[-1]]=len(self.HNAME)-1
                self.Hindx[self.HNAME[-1].upper()]=len(self.HNAME)-1

                
            if corr==0 and cl[1]=='V' and cl[8]=="OK":
                self.VCO.append(float(cl[3])*1e-3)
                
                self.VNAME.append(cl[0])
                self.Vindx[self.VNAME[-1]]=len(self.VNAME)-1
                self.Vindx[self.VNAME[-1].upper()]=len(self.VNAME)-1

            if corr==1 and cl[1]=='H':
                self.Hcorr.append(cl[4])
                self.HcorrNAME.append(cl[0])
                self.Hindx[self.HcorrNAME[-1]]=len(self.HcorrNAME)-1
                self.Hindx[self.HcorrNAME[-1].upper()]=len(self.HcorrNAME)-1
            
            if corr==1 and cl[1]=='V':
                self.Vcorr.append(cl[4])
                self.VcorrNAME.append(cl[0])
                self.Vindx[self.VcorrNAME[-1]]=len(self.VcorrNAME)-1
                self.Vindx[self.VcorrNAME[-1].upper()]=len(self.VcorrNAME)-1 
        self.HCO=array(self.HCO)
        self.VCO=array(self.VCO)




