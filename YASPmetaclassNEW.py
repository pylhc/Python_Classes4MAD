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
        self.HCOerr=[]
        self.Hindx={}
        self.VNAME=[]
        self.VCO=[]
        self.VCOerr=[]
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
                self.HCOerr.append(0.0)
                #self.Hindx[replace(cl[0], '"', '')]=len(self.HNAME)-1
                #self.Hindx[replace(cl[0], '"', '').upper()]=len(self.HNAME)-1
                
                self.HNAME.append(cl[0])
                self.Hindx[self.HNAME[-1]]=len(self.HNAME)-1
                self.Hindx[self.HNAME[-1].upper()]=len(self.HNAME)-1

                
            if corr==0 and cl[1]=='V' and cl[8]=="OK":
                self.VCO.append(float(cl[3])*1e-3)
                self.VCOerr.append(0.0)
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

    def minus(self, b):
      c=YASPtwiss()
      for bpm in self.VNAME:
        aindx=self.Vindx[bpm]
        
        try:
          bindx=b.Vindx[bpm]
          c.VNAME.append(bpm)
          c.Vindx[bpm]=len(c.VNAME)-1
          c.VCO.append(  self.VCO[aindx] -  b.VCO[bindx] )
          c.VCOerr.append( self.VCOerr[aindx] + b.VCOerr[bindx])
        except:
          print "Doing YASP V minus:", bpm, "is not in all files"
      
      for bpm in self.HNAME:
        aindx=self.Hindx[bpm]
        try:
          bindx=b.Hindx[bpm]
          c.HNAME.append(bpm)
          c.Hindx[bpm]=len(c.HNAME)-1
          c.HCO.append(  self.HCO[aindx] -  b.HCO[bindx] )
          c.HCOerr.append( self.HCOerr[aindx] + b.HCOerr[bindx]) 
          
        except:                                                               
          print "Doing YASP H minus:", bpm, "is not in all files"

      return c



     

    def YASPwrite(self,ofile):
      f=open(ofile+'x',"w")
      for i in range(len(self.HNAME)):
        print >>f , self.HNAME[i], self.HCO[i], self.HCOerr[i]
      f.close()
      f=open(ofile+'y',"w")
      for i in range(len(self.VNAME)):
        print >>f , self.VNAME[i], self.VCO[i], self.VCOerr[i]
      f.close()


        




def  YASPaverage(list):
    "Average YASP classes"
    if len(list)==0:
        print "Nothing to average"
        return
    elif len(list)==1:
        return list[0]
    else:
        a=YASPtwiss()
        for bpm in list[0].VNAME:
            
            average=0
            rms=0
            count=0
            
            for files in list:
                flag=1
                
                try:
                    co=files.VCO[files.Vindx[bpm]]            
                    average=average+co
                    rms=rms+co**2
                    count=count+1
                except:
                    print "Doing V average:", bpm, "not found in all files"
                    flag=0

            if flag==1 and count > 0:
                
                a.VNAME.append(bpm)
                a.Vindx[bpm]=len(a.VNAME)-1
                a.VCO.append(average/count)
                a.VCOerr.append( sqrt(rms/count - average**2/count**2 + 1e-18) )

                
             
            count=0
            rms=0
            average=0
            flag=1
            for files in list:
                try:
                    co=files.HCO[files.Hindx[bpm]]            
                    average=average + co
                    rms=rms + co**2
                    count=count+1
                except:
                    print "Doing H average:", bpm, "not found in all files"
                    flag=0


            if flag==1 and count > 0:
                
                a.HNAME.append(bpm)
                a.Hindx[bpm]=len(a.HNAME)-1
                a.HCO.append(average/count)
                a.HCOerr.append( sqrt(rms/count - average**2/count**2  + 1e-18) )

    return a           




  
