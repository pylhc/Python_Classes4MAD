#!/usr/bin/env pythonafs



# Just to make sure that the path to the libraires is defined
import sys
sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy

import pickle
from metaclass import *
from numpy import *
from numpy import dot as matrixmultiply
from numpy.linalg import pinv as generalized_inverse

from os import system
import random,re
#from AllLists import *
#from LinearAlgebra import *
import datetime
import time
import cmath




#############################################
def MakeListCD(x, m, modelcut, errorcut, CorD):
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
    g.write('! minus sign to construct error model (+ sign in changeparmeters.tfs to correct)')
    for var in variables:
        g.write(var+' = '+ var+' - ( '+str(deltafamilie[i])+' );\n')
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
	R=transpose(couple_input2.sensitivity_matrix)
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
			self.dispylist2=dispylist2
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
			f1001r.append(couple2.F1001R[couple2.indx[bpm]])
			f1001i.append(couple2.F1001I[couple2.indx[bpm]])
			f1010r.append(couple2.F1010R[couple2.indx[bpm]])
			f1010i.append(couple2.F1010I[couple2.indx[bpm]])


		for bpm in self.dispylist:
			dy.append(dispy.DY[dispy.indx[bpm]])
		for bpm in self.dispylist2:
			dy.append(dispy2.DY[dispy2.indx[bpm]])


		return array(concatenate([f1001r,f1001i,f1010r,f1010i,dy]))



	###################################################################################
	def computeSensitivityMatrix2(self,x,x2):
	################################################################################
    		#global zerovector, normvector
    		self.zerovector=self.computevector2(x['0'],x2['0'],x['0'],x2['0'])
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
			try:
				xx1=x[var]
			except:
				xx1=x['0']
				xx1.Cmatrix()
			try:
				xx2=x2[var]
			except:
				xx2=x2['0']
				xx2.Cmatrix()
        		vector=self.computevector2(xx1, xx2, xx1, xx2)
        		self.sensitivity_matrix.append((vector-self.zerovector)/self.normvector/incr) # Assuming the same incr...
    		self.sensitivity_matrix=array(self.sensitivity_matrix)
		return self.sensitivity_matrix

