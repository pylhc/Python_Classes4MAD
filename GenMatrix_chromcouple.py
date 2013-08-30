



# Just to make sure that the path to the libraires is defined
import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy
from metaclass import *
#try:
#	from Numeric import *
#	from LinearAlgebra import *
#except:
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
#from AllLists_couple import *
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


			if ((abs(x.Cf1001r[x.indx[bn]]-m.Cf1001r[x.indx[bn]]) < modelcut)
			    and (abs(x.Cf1001i[x.indx[bn]]-m.Cf1001i[x.indx[bn]]) < modelcut)
			    and (x.Cf1001iERR[x.indx[bn]] < errorcut)
			    and (x.Cf1001rERR[x.indx[bn]] < errorcut)):


				t.append(x.NAME[i])


	if cou > 0:
		print "Warning: ", cou, "BPMs removed from data for not beeing in the model"


	return t



########################
def MakePairs(x,m, modelcut=0.1, errorcut=0.027):
###########################################
# not used actually for coupling and Dy correction
	t=[]
        cou=0
	keys=x.__dict__.keys()
	if "PHYMDL" in keys:
		phmdl="PHYMDL"
	else:
		phmdl="PHXMDL"
	for i in range(len(x.NAME)):
		phm=x.__dict__[phmdl][i]

		if (x.STDPH[i] < errorcut and abs(x.PHASE[i]-phm) < modelcut):
			#print x.STDPH[i], errorcut, abs(x.PHASE[i]-phm)
			try:
				m.indx[x.NAME[i].upper()]
				m.indx[x.NAME2[i].upper()]
			except:
				print "Not in Response:", x.NAME[i].upper(), x.NAME2[i].upper()
				cou=cou+1
			else:
				t.append(x.NAME[i]+' '+x.NAME2[i])
			#print "Good:", x.NAME[i].upper(), x.NAME2[i].upper()
		else:
			cou=cou+1
	if cou > 0:
		print "Warning: ", cou, "BPM pairs removed from data for not beeing in the model or having too large error deviations: ", phmdl, modelcut, "STDPH",errorcut, "LEN", len(t)
	return t





##################################################
def writeparams(deltafamilie, variables, app=0, path="./"):
########################################################
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open (path+'/changeparameters_couple', mode)
    f = open (path+'/changeparameters_couple.tfs', mode)
    print "Until here "+path
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








##########################################
def correctcouple(a, chromcouple_input, cut=0.01, app=0, path="./"):
################################################



	R=   transpose(chromcouple_input.sensitivity_matrix)
	vector=chromcouple_input.computevector(a)
	wg=chromcouple_input.wg
#        print wg
	#print couple_input.couplelist


	lc=len(chromcouple_input.couplelist)
	weisvec=array(concatenate(
          [sqrt(wg[0])*ones(lc),
            sqrt(wg[1])*ones(lc),
            sqrt(wg[2])*ones(lc),
            sqrt(wg[3])*ones(lc)]))

	#print "here here",str(weisvec)

	Rnew=transpose(transpose(R)*weisvec)



	delta=-matrixmultiply(generalized_inverse(Rnew,cut),(vector-chromcouple_input.zerovector)/chromcouple_input.normvector)


	writeparams(delta, chromcouple_input.varslist, app,  path=path)

	return [delta, chromcouple_input.varslist]



##########################################
class chromcouple_input:
##########################################
	def __init__(self, varslist, couplelist=[], wg=[1,1,1,1,1]):
			self.varslist=varslist
			self.couplelist=couplelist

#			print couplelist.F1001R
#			sys.exit


			self.wg=wg


	##################################################
	def computevector(self,a):
	##################################################
		Cf1001r=[]
		Cf1001i=[]
		Cf1010r=[]
		Cf1010i=[]

		for bpm in self.couplelist:
			Cf1001r.append(a.Cf1001r[a.indx[bpm]])
			Cf1001i.append(a.Cf1001i[a.indx[bpm]])
			Cf1010r.append(a.Cf1010r[a.indx[bpm]])
			Cf1010i.append(a.Cf1010i[a.indx[bpm]])



		#print concatenate([f1001r,f1001i,f1010r,f1010i,dy])
		vector=array(concatenate([Cf1001r,Cf1001i,Cf1010r,Cf1010i]))

		return array(concatenate([Cf1001r,Cf1001i,Cf1010r,Cf1010i]))



	###################################################################################
	def computeSensitivityMatrix(self,x):
	################################################################################
    		#global zerovector, normvector
    		self.zerovector=self.computevector(x['0'])
    		self.sensitivity_matrix=[]
    		incr=x['incr'][0]  # BUG! need to read it from FullResponse!

    		#ncouple=4*len(self.couplelist)
    		self.normvector=array(concatenate([ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.couplelist)),ones(len(self.couplelist))]))*1.0
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




