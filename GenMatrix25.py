#!/usr/bin/env pythonafs



# Just to make sure that the path to the libraires is defined 
import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy

import pickle
from numpy import *
#from numpy.oldnumeric.linear_algebra import generalized_inverse
from os import system
from metaclass25 import twiss
import random,re
from AllLists import *
#from LinearAlgebra import *


########################
def MakePairs(x):
###########################################
        t=[]
        for i in range(len(x.NAME)):
                t.append(x.NAME[i]+' '+x.NAME2[i])
        return t
						

##################################################
def writeparams(deltafamilie, variables, app=0):
########################################################
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    g = open ('changeparameters', mode)
    i=0
    print len(variables), len(deltafamilie)
    for var in variables:      
        g.write(var+' = '+ var+' + ( '+str(deltafamilie[i])+' );\n')
        i +=1
    g.close()
    return



#----calc. beta,phase,disp beat
def betabeat(a,b):
    rmsx=sqrt(sum(((a.BETX-b.BETX)/b.BETX)**2)/len(b.BETX))
    rmsy=sqrt(sum(((a.BETY-b.BETY)/b.BETY)**2)/len(b.BETY))
    peakx=max(abs((a.BETX-b.BETX)/b.BETX))
    peaky=max(abs((a.BETY-b.BETY)/b.BETY))
    rmsdx=sqrt(sum(((a.DX/sqrt(a.BETX)-b.DX/sqrt(b.BETX)))**2)/len(b.DX))
    peakdx=max(abs((a.DX/sqrt(a.BETX)-b.DX/sqrt(b.BETX))))
    dphixa=[];dphiya=[];dphixb=[];dphiyb=[]
    for i in range(1,len(a.MUX)):
        dphixa.append(a.MUX[i]-a.MUX[i-1])
        dphiya.append(a.MUY[i]-a.MUY[i-1])
        dphixb.append(b.MUX[i]-b.MUX[i-1])
        dphiyb.append(b.MUY[i]-b.MUY[i-1])
    dphixa=array(dphixa);dphiya=array(dphiya)
    dphixb=array(dphixb);dphiyb=array(dphiyb)
    rmsphix=sqrt(sum((dphixa-dphixb)**2)/len(dphixa))
    rmsphiy=sqrt(sum((dphiya-dphiyb)**2)/len(dphiya))
    peakphix=max(abs(dphixa-dphixb));
    peakphiy=max(abs(dphiya-dphiyb))

    #return array([rmsx,rmsy, peakx, peaky])
    return array([rmsx,rmsy, peakx,peaky, rmsphix,rmsphiy,
                  peakphix,peakphiy, rmsdx, peakdx])








def correctbeat(a,beat_input, cut=0.01, app=0):
    #global zerovector, R, sigmaphase, dx, sensitivity_matrix, normvector
    
    R=   transpose(beat_input.sensitivity_matrix)  #Transpose is dangerous, shares memory with initial array!!!
    vector=beat_input.computevector(a)
    #print len(vector), len(zerovector), len(generalized_inverse(R,cut)), len(R[0]) 

    #print sqrt(wg[0]), sqrt(wg[1]), sqrt(wg[2],sqrt(wg[3])),'\n'
    #---- weights for phix, phiy,betax, betay,  disp,tunes
    wg=beat_input.wg
    weisvec=array(concatenate([sqrt(wg[0])*ones(len(beat_input.phasexlist)),
                               sqrt(wg[1])*ones(len(beat_input.phaseylist)),
                               sqrt(wg[2])*ones(len(beat_input.betaxlist)),
                               sqrt(wg[3])*ones(len(beat_input.betaylist)),
                               sqrt(wg[4])*ones(len(beat_input.displist)),
                               sqrt(wg[5])*ones(2)]))
    #---- equiv to x=inv(A^T*W*A)*A^T*W*b, W is the weight matrix
    Rnew=transpose(transpose(R)*weisvec)
    #for lj in range(len(weisvec)):
    #    R[lj,:] = R[lj,:]*weisvec[lj]


    delta=-dot(linalg.pinv(Rnew,cut),(vector-beat_input.zerovector)/beat_input.normvector)

    writeparams(delta, beat_input.varslist, app)
    return




######################
class beat_input:
##########################################
	def __init__(self, varslist, phasexlist=[], phaseylist=[], betaxlist=[], betaylist=[], displist=[], wg=[1,1,1,1,1,10]):
			self.varslist=varslist
			self.phasexlist=phasexlist
			self.phaseylist=phaseylist
			self.betaxlist=betaxlist
			self.betaylist=betaylist
			self.displist=displist
			self.wg=wg


	#####################
	def computevector(self,x):
	##################################################################################
    		phix=[]
 		phiy=[]		
    		betx=[]
    		bety=[]
    		disp=[]
    		for ph in self.phasexlist:
        		[bpm1,bpm2]=ph.split()
        		phix.append(x.MUX[x.indx[bpm2]]-x.MUX[x.indx[bpm1]])
    		for ph in self.phaseylist:
        		[bpm1,bpm2]=ph.split()
        		phiy.append(x.MUY[x.indx[bpm2]]-x.MUY[x.indx[bpm1]])
    		for b in self.betaxlist:
        		betx.append(x.BETX[x.indx[b]])
    		for b in self.betaylist:
        		bety.append(x.BETY[x.indx[b]])
    		for b in self.displist:
        		disp.append(x.DX[x.indx[b]]/sqrt(x.BETX[x.indx[b]]))
    		return array(concatenate([phix,phiy,betx,bety,disp,[x.Q1],[x.Q2]]))


	###################################################################################
	def computeSensitivityMatrix(self,x):
	################################################################################
    		#global zerovector, normvector
    		self.zerovector=self.computevector(x['0'])
    		self.sensitivity_matrix=[]
    		incr=x['incr'][0]  # BUG! need to read it from FullResponse!

    		nph=len(self.phasexlist)+len(self.phaseylist)
    		nbet=len(self.betaxlist)+len(self.betaylist)
    		ndisp=len(self.displist)
    		self.normvector=array(concatenate([ones(len(self.phasexlist)), ones(len(self.phaseylist)), self.zerovector[nph:nph+nbet], ones(ndisp+2)]))*1.0  #To take dbeta/beta later
    		for var in self.varslist:
        		vector=self.computevector(x[var])
        		self.sensitivity_matrix.append((vector-self.zerovector)/self.normvector/incr)
    		self.sensitivity_matrix=array(self.sensitivity_matrix)
		return self.sensitivity_matrix
    

########### START ###############
'''



print "Starting loading Full Response optics"
FullResponse=pickle.load(open('/afs/cern.ch/user/r/rtomas/w1/MAD/LHC/Beta-Beat/NewRespMatrixGen/FullResponse','r'))
print "Loading ended"

varslist=quadvarsb1()
variables=varslist
phasexlist=[]
phaseylist=[]
betaxlist=IndependentQuadsb1()
betaylist=betaxlist
displist=BPMsb1()
wei=[1,1,1,1,1,10] # Weights of phasex phasey betax betay disp and tunes

beat_inp=beat_input(varslist, phasexlist, phaseylist, betaxlist, betaylist, displist, wei)



sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
    
y=twiss('twiss.base')

x=twiss('twiss.dat')

#--- evaluate beta-beat for output
print "Initial beatings:"
print "rmsx  rmsy  peakx  peaky  rmsphix  rmsphiy  peakphix  peakphiy   rmsdx  peakdx"
print betabeat(x,y)

#------ Apply svd correction
correctbeat(x, beat_inp, cut=0.04, app=0)
    
#----- compute twiss after correction
#system('madx < job.random.madx > scum  ')
#z=twiss('twiss.dat')
#bb1=betabeat(x,y);bb=betabeat(z,y)
    
'''
    
    

    
