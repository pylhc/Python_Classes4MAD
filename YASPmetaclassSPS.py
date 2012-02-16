#!/usr/bin/env python

from Numeric import *
from string import split
from string import replace
from LinearAlgebra import *
import sys
import gzip

#########################
class YASPtwiss:
#########################
    "Read orbit data from YASP output in the same way as metaclass"
    def __init__(self, filename, dictionary={}): 

        if len(dictionary) ==0:
          print "If SPS: I need the dictionary for YASPtwiss"
          print "If LHC: I continue"
          #sys.exit()
        
        if '.gz' in filename:
            f=gzip.open(filename, 'rb')
        else:
            f=open(filename, 'r')



        fspl=filename.split('/')[-1].split('.')
        self.exciter=fspl[1]+fspl[2]
        self.plane=fspl[3]
        if "minus" in fspl[4]:
          self.str=-1.0
        else:
          self.str=1.0
        

        self.HNAME=[]
        self.HCO=[]
        self.Hindx={}
        self.VNAME=[]
        self.VCO=[]
        self.Vindx={}
        corr=0
        count=0
        for line in f:
            cl=line.split()
            if ("@ " in line and "%le" in line) :
                label=cl[1]
		try:
                	exec "self."+label+"= "+str(float(cl[3]))
		except:
			print "Problem parsing:", line
			print "Ignored"

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
                if len(dictionary)==0:
                  clname=cl[0]
                else:
                  clname=dictionary[cl[0]]
                
                self.HNAME.append(clname)
                self.Hindx[self.HNAME[-1]]=len(self.HNAME)-1
                self.Hindx[self.HNAME[-1].upper()]=len(self.HNAME)-1

                
            if corr==0 and cl[1]=='V' and cl[8]=="OK":
                self.VCO.append(float(cl[3])*1e-3)
                
                if len(dictionary)==0:
                    clname=cl[0]
                else:
                    clname=dictionary[cl[0]]
                
                self.VNAME.append(clname)
                self.Vindx[self.VNAME[-1]]=len(self.VNAME)-1
                self.Vindx[self.VNAME[-1].upper()]=len(self.VNAME)-1

            #if corr==1 and len(cl)>4 and 'STR' not in cl[4] and abs(float(cl[4])) > 0:
               #count=count+1
               #excit=cl[0]
               #excit=excit.replace(".","")
               #self.exciter=excit
               #self.plane=cl[1]
               #self.str=float(cl[4])
               #if count > 1:
                   #print "Error: YASP file using 2 correctors:", filename
                   #sys.exit()
