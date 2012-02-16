#!/usr/bin/env pythonafs



# Just to make sure that the path to the libraires is defined 
import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')

 
#--- beta beat for store with numpy
try:
	from metaclass import *
except:
	from metaclass25 import *
try:
	from Numeric import *
	from LinearAlgebra import *
except:
	from numpy import *
	from numpy import linalg
	from numpy.linalg import pinv as generalized_inverse
	from numpy import dot as matrixmultiply

import pickle
#from Numeric import *
#from numpy.oldnumeric.linear_algebra import generalized_inverse
from os import system
#from metaclass import twiss
import random,re
#from AllLists import *
#from LinearAlgebra import *
import datetime
import time
import cmath




#############################################
def MakeList(x, m, modelcut, errorcut, CorD):
#############################################
        """Makes list in coupling correction """
        t=[]
        cou=0
	if x==[]:
		return []
	for i in range(len(x.NAME)):
		try:
			m.indx[x.NAME[i].upper()]
			bn=x.NAME[i].upper()
		except:
			print "Not in Response:", x.NAME[i].upper()
			cou=cou+1
		else:
			if (CorD=='C'):
				
				if ((abs(x.F1001W[x.indx[bn]]-m.F1001W[x.indx[bn]]) < modelcut) and (x.FWSTD1[x.indx[bn]] < errorcut)):
			
				
					t.append(x.NAME[i])
					
			elif(CorD=='D'):
				if ((abs(x.DY[x.indx[bn]]-m.DY[x.indx[bn]]) < modelcut) and (x.STDDY[x.indx[bn]] < errorcut)):
					t.append(x.NAME[i])
	if cou > 0:
		print "Warning: ", cou, "BPMs removed from data for not beeing in the model"

	
	return t	




##################################################
def writeparams(deltafamilie, variables, app=0, path="./"):
########################################################
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open (path+'changeparameters', mode)
    f = open (path+'changeparameters.tfs', mode)
    print >>f, "@", "APP", "%le", app
    print >>f, "@", "PATH","%s", path
    print >>f, "@", "DATE", "%s", a.ctime()
    print >>f, "*", "NAME", "DELTA"
    print >>f, "$", "%s", "%le"
    i=0
    print len(variables), len(deltafamilie)
    for var in variables:
        g.write(var+' = '+ var+' + ( '+str(deltafamilie[i])+' );\n')
	f.write(var+'   '+str(deltafamilie[i])+'\n')
        i +=1
    g.close();f.close()
    return




#################################################
def coupling(a,b):
#################################################
	# Applicable to both simulation-model and exp-model

	rmsf1001=sqrt(sum((a.F1001-b.F1001)**2)/len(b.F1001))
	rmsf1010=sqrt(sum((a.F1010-b.F1010)**2)/len(b.F1001))
	peakf1001=max(abs(a.F1001-b.F1001))
	peakf1010=max(abs(a.F1010-b.F1010))
	couplelist=MakePairs(a, b)
	return array([rmsf1001, rmsf1010, peakf1001, peakf1010])



##########################################
def correctcouple(a, dispy, couple_input, cut=0.01, app=0, path="./"):
################################################
	print "one"
	R=   transpose(couple_input.sensitivity_matrix) 
	vector=couple_input.computevector(a,dispy)
	wg=couple_input.wg
	print len(couple_input.couplelist), wg[2]
	weisvec=array(concatenate([sqrt(wg[0])*ones(len(couple_input.couplelist)),sqrt(wg[1])*ones(len(couple_input.couplelist)),sqrt(wg[2])*ones(len(couple_input.couplelist)),sqrt(wg[3])*ones(len(couple_input.couplelist)),sqrt(wg[4])*ones(len(couple_input.dispylist))]))
	Rnew=transpose(transpose(R)*weisvec)
	delta=-matrixmultiply(generalized_inverse(Rnew,cut),(vector-couple_input.zerovector)/couple_input.normvector)
	writeparams(delta, couple_input.varslist, app,  path=path)
	return [delta, couple_input.varslist]



##########################################
class couple_input:
##########################################
	def __init__(self, varslist, couplelist=[], dispylist=[], wg=[10,10,1,1,1]):
			self.varslist=varslist
			self.couplelist=couplelist
			self.dispylist=dispylist
			self.wg=wg


	##################################################
	def computevector(self,a,dispy):
	##################################################
		f1001r=[]
		f1001i=[]
		f1010r=[]
		f1010i=[]
	       	dy=[]
		for bpm in self.couplelist:
			f1001r.append(a.F1001R[a.indx[bpm]])
			f1001i.append(a.F1001I[a.indx[bpm]])
			f1010r.append(a.F1010R[a.indx[bpm]])
			f1010i.append(a.F1010I[a.indx[bpm]])
		for bpm in self.dispylist:
			dy.append(dispy.DY[dispy.indx[bpm]])

				
		return array(concatenate([f1001r,f1001i,f1010r,f1010i,dy]))



	###################################################################################
	def computeSensitivityMatrix(self,x):
	################################################################################
    		#global zerovector, normvector
    		self.zerovector=self.computevector(x['0'],x['0'])
    		self.sensitivity_matrix=[]
    		incr=x['incr'][0]  # BUG! need to read it from FullResponse!

    		#ncouple=4*len(self.couplelist)
    		self.normvector=array(concatenate([ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.dispylist)) ]))*1.0
    		for var in self.varslist:
        		vector=self.computevector(x[var], x[var])
        		self.sensitivity_matrix.append((vector-self.zerovector)/self.normvector/incr)
    		self.sensitivity_matrix=array(self.sensitivity_matrix)
		return self.sensitivity_matrix
    


##########################################
def correctcouple2(couple, couple2, dispy, dispy2, couple_input2, cut=0.01, app=0, path="./"):
################################################
	R=transpose(couple_input.sensitivity_matrix2) 
	vector=couple_input2.computevector2(couple, couple2, dispy, dispy2)
	wg=couple_input2.wg
	weisvec=array(concatenate([sqrt(wg[0])*ones(len(couple_input2.couplelist)),
				   sqrt(wg[0])*ones(len(couple_input2.couplelist2)),
				   sqrt(wg[1])*ones(len(couple_input2.couplelist)),
				   sqrt(wg[1])*ones(len(couple_input2.couplelist2)),
				   sqrt(wg[2])*ones(len(couple_input2.couplelist)),
				   sqrt(wg[2])*ones(len(couple_input2.couplelist2)),
				   sqrt(wg[3])*ones(len(couple_input2.couplelist)),
				   sqrt(wg[3])*ones(len(couple_input2.couplelist2)),
				   sqrt(wg[4])*ones(len(couple_input2.dispylist)),
				   sqrt(wg[4])*ones(len(couple_input2.dispylist2))]))
	Rnew=transpose(transpose(R)*weisvec)
	delta=-matrixmultiply(generalized_inverse(Rnew,cut),(vector-couple_input2.zerovector)/couple_input2.normvector)
	writeparams(delta, couple_input2.varslist, app,  path=path)
	return [delta, couple_input2.varslist]



##########################################
class couple_input2:
##########################################
	def __init__(self, varslist, couplelist=[], couplelist2=[], dispylist=[], dispylist2=[], wg=[1,1,1,1,1]):
			self.varslist=varslist
			self.couplelist=couplelist
			self.couplelist2=couplelist2
			self.dispylist=dispylist
			self.dispylist=dispylist2     	
			self.wg=wg


	##################################################
	def computevector2(self,couple,couple2,dispy,dispy2):
	##################################################
		f1001r=[]
		f1001i=[]
		f1010r=[]
		f1010i=[]
	       	dy=[]
		for bpm in self.couplelist:
			f1001r.append(couple.F1001R[couple.indx[bpm]])
			f1001i.append(couple.F1001I[couple.indx[bpm]])
			f1010r.append(couple.F1010R[couple.indx[bpm]])
			f1010i.append(couple.F1010I[couple.indx[bpm]])
		for bpm in self.couplelist2:
			f1001r.append(couple2.F1001R[couple.indx[bpm]])
			f1001i.append(couple2.F1001I[couple.indx[bpm]])
			f1010r.append(couple2.F1010R[couple.indx[bpm]])
			f1010i.append(couple2.F1010I[couple.indx[bpm]])
		
			
		for bpm in self.dispylist:
			dy.append(dispy.DY[dispy.indx[bpm]])
		for bpm in self.dispylist2:
			dy.append(dispy2.DY[dispy2.indx[bpm]])
			

		return array(concatenate([f1001r,f1001i,f1010r,f1010i,dy]))



	###################################################################################
	def computeSensitivityMatrix2(self,x,x2):
	################################################################################
    		#global zerovector, normvector
    		self.zerovector=self.computevector(x['0'],x2['0'],x['0'],x2['0'])
    		self.sensitivity_matrix=[]
    		incr=x['incr'][0]  # BUG! need to read it from FullResponse!
		incr2=x2['incr'][0]

    		#ncouple=4*len(self.couplelist)
    		self.normvector=array(concatenate([ones(len(self.couplelist)),
						   ones(len(self.couplelist2)),
						   ones(len(self.couplelist)),
						   ones(len(self.couplelist2)),
						   ones(len(self.couplelist)),
						   ones(len(self.couplelist2)),
						   ones(len(self.couplelist)),
						   ones(len(self.couplelist2)),
						   ones(len(self.dispylist)),
						   ones(len(self.dispylist2)),]))*1.0
    		for var in self.varslist:
        		vector=self.computevector2(x[var], x2[var], x[var], x2[var])
        		self.sensitivity_matrix.append((vector-self.zerovector)/self.normvector/incr) # Assuming the same incr...
    		self.sensitivity_matrix=array(self.sensitivity_matrix)
		return self.sensitivity_matrix






###################################################
def correctbeatEXP2(x,x2,y,y2,dx,dx2,beat_input2, cut=0.01, app=0, path="./"):
###########################################################
	R=   transpose(beat_input2.sensitivity_matrix)
        vector=beat_input2.computevectorEXP2(x,x2,y,y2,dx,dx2)
        wg=beat_input2.wg
	weisvec=array(concatenate([sqrt(wg[0])*ones(len(beat_input2.phasexlist)+len(beat_input2.phasexlist2)),
		sqrt(wg[1])*ones(len(beat_input2.phaseylist)+len(beat_input2.phaseylist2)),
                sqrt(wg[2])*ones(len(beat_input2.betaxlist)+len(beat_input2.betaxlist2)),
	        sqrt(wg[3])*ones(len(beat_input2.betaylist)+len(beat_input2.betaylist2)),
		sqrt(wg[4])*ones(len(beat_input2.displist)+len(beat_input2.displist2)),
		sqrt(wg[5])*ones(4)]))
	Rnew=transpose(transpose(R)*weisvec)
	#print (vector-beat_input.zerovector) # Good line to debug big correction strs
	delta=-matrixmultiply(generalized_inverse(Rnew,cut),(vector-beat_input2.zerovector)/beat_input2.normvector)
        writeparams(delta, beat_input2.varslist, app, path)
	return [delta, beat_input2.varslist]





######################
class beat_input2:
##########################################
	def __init__(self, varslist, phasexlist=[], phaseylist=[], betaxlist=[], betaylist=[], displist=[], phasexlist2=[], phaseylist2=[], betaxlist2=[], betaylist2=[], displist2=[], wg=[1,1,1,1,1,10]):
			self.varslist=varslist
			self.phasexlist=phasexlist
			self.phaseylist=phaseylist
			self.betaxlist=betaxlist
			self.betaylist=betaylist
			self.displist=displist
			self.phasexlist2=phasexlist2
			self.phaseylist2=phaseylist2
			self.betaxlist2=betaxlist2
			self.betaylist2=betaylist2
			self.displist2=displist2
			self.wg=wg


	################################
	def computevectorEXP2(self,x,x2,y,y2,dx,dx2):
	###############################################################################
		phix=[]
		phiy=[]
		betx=[]
		bety=[]
		disp=[]
		print "Exp tunes in beta_input class: ", x.Q1, y.Q2, x2.Q1, y2.Q2
		for ph in self.phasexlist:
			[bpm1,bpm2]=ph.split()
			phix.append(x.PHASEX[x.indx[bpm1]])
		for ph in self.phaseylist:
                        [bpm1,bpm2]=ph.split()
                        phiy.append(y.PHASEY[y.indx[bpm1]])
		for b in self.displist:
			disp.append(dx.NDX[dx.indx[b]])

		for ph in self.phasexlist2:
			[bpm1,bpm2]=ph.split()
			phix.append(x2.PHASEX[x2.indx[bpm1]])
		for ph in self.phaseylist2:
                        [bpm1,bpm2]=ph.split()
                        phiy.append(y2.PHASEY[y2.indx[bpm1]])
		for b in self.displist2:
			disp.append(dx2.NDX[dx2.indx[b]])
		
		return array(concatenate([phix,phiy,betx,bety,disp,[x.Q1],[y.Q2],[x2.Q1],[y2.Q2]]))



	#####################
	def computevector2(self,x,x2):
	##################################################################################
    		phix=[]
 		phiy=[]		
    		betx=[]
    		bety=[]
    		disp=[]
		#print "Model tunes in beta_input class: ", x.Q1, x.Q2
    		for ph in self.phasexlist:
        		[bpm1,bpm2]=ph.upper().split()
        		phix.append(x.MUX[x.indx[bpm2]]-x.MUX[x.indx[bpm1]])
    		for ph in self.phaseylist:
        		[bpm1,bpm2]=ph.upper().split()
        		phiy.append(x.MUY[x.indx[bpm2]]-x.MUY[x.indx[bpm1]])
    		for b in self.betaxlist:
        		betx.append(x.BETX[x.indx[b]])
    		for b in self.betaylist:
        		bety.append(x.BETY[x.indx[b]])
    		for b in self.displist:
        		disp.append(x.DX[x.indx[b]]/sqrt(x.BETX[x.indx[b]]))
		for ph in self.phasexlist2:
        		[bpm1,bpm2]=ph.upper().split()
        		phix.append(x2.MUX[x2.indx[bpm2]]-x2.MUX[x2.indx[bpm1]])
    		for ph in self.phaseylist2:
        		[bpm1,bpm2]=ph.upper().split()
        		phiy.append(x2.MUY[x2.indx[bpm2]]-x2.MUY[x2.indx[bpm1]])
    		for b in self.betaxlist2:
        		betx.append(x2.BETX[x2.indx[b]])
    		for b in self.betaylist2:
        		bety.append(x2.BETY[x2.indx[b]])
    		for b in self.displist2:
        		disp.append(x2.DX[x2.indx[b]]/sqrt(x2.BETX[x2.indx[b]]))

    		return array(concatenate([phix,phiy,betx,bety,disp,[x.Q1],[x.Q2],[x2.Q1],[x2.Q2]]))



	###################################################################################
	def computeSensitivityMatrix2(self,x,x2):
	################################################################################
    		#global zerovector, normvector
    		self.zerovector=self.computevector2(x['0'],x2['0'])
    		self.sensitivity_matrix=[]
    		incr=x['incr'][0]  # BUG! need to read it from FullResponse!
		incr2=x2['incr'][0]

    		nph=len(self.phasexlist)+len(self.phaseylist)+len(self.phasexlist2)+len(self.phaseylist2)
    		nbet=len(self.betaxlist)+len(self.betaylist)+len(self.betaxlist2)+len(self.betaylist2)
    		ndisp=len(self.displist)
		ndisp2=len(self.displist2)
    		self.normvector=array(concatenate([ones(len(self.phasexlist)), ones(len(self.phasexlist2))*incr/incr2, ones(len(self.phaseylist)), ones(len(self.phaseylist2))*incr/incr2, self.zerovector[nph:nph+nbet], ones(ndisp), ones(ndisp2)*incr/incr2,ones(4)]))*1.0  #To take dbeta/beta later
    		for var in self.varslist:
			try:
				xx1=x[var]
			except:
				xx1=x['0']
			try:
				xx2=x2[var]
			except:
				xx2=x2['0']
        		vector=self.computevector2(xx1,xx2)
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
    
    

    
