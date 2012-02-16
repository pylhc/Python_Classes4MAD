                  #############################################
                  #                                           #
                  # @ author : Glenn Vanbavinckhove (02/07/09)#
                  #                                           #
                  #                                           #
                  #############################################

                  
import os,sys
from metaclass import twiss
from string import *
from LinearAlgebra import *
import linreg as linreg
from fit import FIT
try:
    from numpy import array
    
    print "Python will use numpy"
    
except ImportError:

     from Numeric import *

     print "Python wil use numeric"


####### info
#
# features =>
# - readtables(path,ext) => read twiss files with the given extension
# - writetables(name,info,names,specie,data) => writes twiss tables to specific name, info (@), names (column names), specie (%le,%s,...) and data

#### global functions
def intersect(ListOfFile):
    print "Find common BPMS in "+str(len(ListOfFile))+" files"
    if len(ListOfFile)==0:
       	print "Nothing to intersect!!!!"
       	sys.exit()
    z=ListOfFile[0].NAME
    for b in ListOfFile:
       	z=filter(lambda x: x in z   , b.NAME)
	#SORT by S
    result=[]
    x0=ListOfFile[0]
    for bpm in z:
       	result.append((x0.S[x0.indx[bpm]], bpm)) # first is location, second name

    print str(len(result))+" common BPMs have been found"
    
    result.sort()

    return result

def modelIntersect(expbpms, model):
	bpmsin=[]
	for bpm in expbpms:
		try:
			check=model.indx[bpm[1].upper()]
			bpmsin.append(bpm)
		except:
			print bpm, "Not in Model"
	if len(bpmsin)==0:
		print "Zero intersection of Exp and Model"
		print "Please, provide a good Dictionary"
		print "Now we better leave!"
		sys.exit()			
	return bpmsin

def buildString(write):
  stri =""
  for i in range(len(write)):
    if i==0:
        stri +=str(write[i])
    else:
        stri +='  '+str(write[i])
    
  return stri



#### class getData
class getData:

    ###### initialization
    def __init__(self):

        print "GETDATA =>initialising opticsclass"

    ###### handling data stuff
    def readtables(self,path,ext):

        self.twissfilesx=[]
        self.twissfilesy=[]
        
        
        files=os.listdir(path)

        for i in files:

            if ext in i:
                
                print "Twiss file found : "+path+"/"+i
                tempx=twiss(path+"/"+i)
                tempy=twiss(path+"/"+i)
        
                self.twissfilesx.append(tempx)
                self.twissfilesy.append(tempx)

        print "Number of twiss files found in directory "+str(len(self.twissfiles))

    ###### writing data to file
    def writetables(self,name,info,names,specie,data):

        file2write=open(name,'w')

        for i in info: # writing extra info to the table (@ ....  % ....)
            
            file2write.write(i+'\n')

        file2write.write(buildString(names)+'\n')
        file2write.write(buildString(specie)+'\n')
        bpms=data[0]
        for i in range(len(bpms)):
            bn=upper(bpms[i][1])

            appdata=[]
    
            for j in range(len(data)):

                if j==0:
                    appdata.append("\""+bpms[i][1]+"\"")
                    appdata.append(bpms[i][0])
                else:
                    appdata.append(data[j][bn])
            

            file2write.write(buildString(appdata)+'\n')

        file2write.close() 


    ###### linear stuff
class getLin:

    def getphases(self):

        print "METHOD => Starting to calculate phase advance"

    def getbetas(self):

        print "METHOD => Starting to calculate betas from phases"




#### class getNonLin
class getNonLin:
        
    #=> calculating fterms
    def getsextupole(self):

        print " Start to calculate the fterms for sextupoles"

class getOffMomentum:

    #=> calculating off-momentum beta-beat (06/07/2009)
    def findoffBB(self,t,plane,MADTwiss,betalist):

        print "METHOD => Starting to calculate off-momentum beta-beat"

        
	self.bpms=intersect(t)
	self.bpms=modelIntersect(self.bpms, MADTwiss)
	MADTwiss.chrombeat()

	self.slope={}
        self.slopeE={}
	self.slopeM={}

        dpplist=[]

        for i in range(len(betalist)): #### for sorting

            dpplist.append(betalist[i]['DPP'])

        
        dpplist.sort()
        
        #if dpplist[0]>dpplist[1]
        #dpplist.reverse()

        
	for i in range(len(self.bpms)):
            
            bn=upper(self.bpms[i][1])
            check=0
            slopei=0.0
            slopeiM=0.0
            beta=[]
            k=0
            

            for k in range(len(dpplist)):
                count=0
                for j in range(len(betalist)):
                    
                    if dpplist[k]==betalist[j]['DPP']:

                        if count==0:
                        
                            beta.append(betalist[j][bn][0])
                        count=1
                        
          
                  
            [slope, offset, R, slope_error, offset_error]=linreg.linreg(dpplist,beta)
           
            a=FIT(dpplist,beta,1,1e-15)
            #print ' list '+str(dpplist)
            coef,error=a.linreg()             
			
            if plane=='H':
                slopeiM=MADTwiss.dbx[MADTwiss.indx[upper(bn)]]
            else:
                slopeiM=MADTwiss.dby[MADTwiss.indx[upper(bn)]]
	
            #self.slope[bn]=coef[0]
            self.slope[bn]=coef[0]/100
            #self.slopeE[bn]=error[0]
            self.slopeE[bn]=error[0]/100
            self.slopeM[bn]=slopeiM
      
                

    #=> calculating chromatacity (03/07/2009)
    def findchrom(self,t,plane):

        print "METHOD => Starting to calculate the chromaticity "

        

        # t is collection of twiss files

        if len(t)<2:
            print " Not enough files \n ... I cannot work like this! => sys.exit()"
            sys.exit()


        bpms=intersect(t)

        slopem=0.0
        error=[]
        slopeme=0.0

        self.slope={}
        self.slope_error={}
        self.chrom=0.0
        self.chrome=0.0
        self.bpms=bpms
        
        for x in range(len(bpms)):
             bpm=upper(bpms[x][1])
             chrom=0.0
             chrom_error=0.0
             tune=[]
             dpp=[]

             for i in range(len(t)):

                 if plane=='H':
                     tune.append(t[i].TUNEX[t[i].indx[bpm]])
                     dpp.append(t[i].DPP)
                 else:
                     tune.append(t[i].TUNEY[t[i].indx[bpm]])
                     dpp.append(t[i].DPP)
                     
             
             [slope, offset, R, slope_error, offset_error]=linreg.linreg(dpp,tune)

             slopem+=slope
             slopeme+=slope_error
             self.slope[bpm]=float(slope)
             self.slope_error[bpm]=(float(slope_error))

        self.chrom=slopem/len(bpms)
        self.chrome=slopeme/len(bpms)
       
        #=> for finding off momentum phase
    
    def findoffPhase(self,t,plane,MADTwiss,phaselist):

        print "METHOD => Start calculating off-momentum phase"

        # t is collection of twiss files
        dpplist=[]
        for i in range(len(phaselist)): #### for sorting

            dpplist.append(phaselist[i]['DPP'])
        dpplist.sort()

        if len(t)<2:
            print " Not enough files \n ... I cannot work like this! => sys.exit()"
            sys.exit()

        self.bpms=intersect(t)
        self.bpms=modelIntersect(self.bpms, MADTwiss)
         
        self.slope={}
        self.slopeE={}
        self.slopeM={}

        print "The lenght of the bpms "+str(len(self.bpms))

        for i in range(len(self.bpms)):
            print i
            if i==len(self.bpms)-1:
                bn=upper(self.bpms[i][1])
                print bn
                print "hereeeeeee"
                bn2=upper(self.bpms[0][1])
            else:
                bn=upper(self.bpms[i][1])
                bn2=upper(self.bpms[i+1][1])
                
            phase=[]

            if plane=='H':
                slopeMp=MADTwiss.DMUX[MADTwiss.indx[bn2]]-MADTwiss.DMUX[MADTwiss.indx[bn]]	
            else:
                slopeMp=MADTwiss.DMUY[MADTwiss.indx[bn2]]-MADTwiss.DMUY[MADTwiss.indx[bn]]
       

            for k in range(len(dpplist)):
                count=0
                for j in range(len(phaselist)):

                    if dpplist[k]==phaselist[j]['DPP']:

                        if count==0:
                       
                            phase.append(phaselist[j][bn][0])
                            
                        count=1
                        
            a=FIT(dpplist,phase,1,1e-15)
            coef,error=a.linreg()

            self.slope[bn]=coef[0]
            self.slopeE[bn]=error[0]
            self.slopeM[bn]=slopeMp
