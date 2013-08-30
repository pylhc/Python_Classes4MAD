#!/usr/bin/env pythonafs



# Just to make sure that the path to the libraires is defined
import sys
#sys.path.append('/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/')


#--- beta beat for store with numpy
from metaclass import *
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




################################
def MakeList(x,m):
################################
	t=[]
        cou=0
	if x==[]:
		return []
	for i in range(len(x.NAME)):
		try:
			m.indx[x.NAME[i].upper()]
		except:
			print "Not in Response:", x.NAME[i].upper()
			cou=cou+1
		else:
			t.append(x.NAME[i])
	if cou > 0:
		print "Warning: ", cou, "BPMs removed from data for not beeing in the model"
	return t


########################
def MakePairs(x,m, modelcut=0.1, errorcut=0.027):
###########################################
	t=[]
        cou=0
	keys=x.__dict__.keys()
	if "PHYMDL" in keys:
		phmdl="PHYMDL"
		STD=x.STDPHY
		PHASE=x.PHASEY

	else:
		phmdl="PHXMDL"
		STD=x.STDPHX
		PHASE=x.PHASEX

	for i in range(len(x.NAME)-1):    # The -1 comes from the fact that the last BPM is again the first BPM, model not ready for this
		phm=x.__dict__[phmdl][i]

		if (STD[i] < errorcut and abs(PHASE[i]-phm) < modelcut):
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


#######################################################
def MakePairsBeta(x,m, modelcut=5.0, errorcut=5.0):
######################################################
	t=[]
        cou=0
	keys=x.__dict__.keys()
	if "BETYMDL" in keys:
		betmdl="BETYMDL"
		plane="V"
	else:
		betmdl="BETXMDL"
		plane="H"
	for i in range(len(x.NAME)):
		betm=x.__dict__[betmdl][i]
		if (plane=="H"):
			if (x.STDBETX[i] < errorcut and abs(x.BETX[i]-betm) < modelcut):
				try:
					m.indx[x.NAME[i].upper()]
					t.append(x.NAME[i])
				except:
					print "Not in Response:", x.NAME[i].upper()
					cou=cou+1
			else:
				cou=cou+1
		if (plane=="V"):
			if (x.STDBETX[i] < errorcut and abs(x.BETY[i]-betm) < modelcut):
				try:
					m.indx[x.NAME[i].upper()]
					t.append(x.NAME[i])
				except:
					print "Not in Response:", x.NAME[i].upper()
					cou=cou+1
			else:
				cou=cou+1

	if cou > 0:
		print "Warning: ", cou, "BPM removed from data for not beeing in the model or having too large error deviations: ", betmdl, modelcut, "STDBET",errorcut, "LEN", len(t)
	return t



##################################################
def writeparams(deltafamilie, variables, app=0, path="./"):
########################################################
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open (path+'/changeparameters', mode)
    f = open (path+'/changeparameters.tfs', mode)
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
def betabeat(a,b):
#############################################################
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


#################################################
def betabeatEXP(x,y,b):
####################################################
    #rmsx=sqrt(sum(((a.BETX-b.BETX)/b.BETX)**2)/len(b.BETX))
    #rmsy=sqrt(sum(((a.BETY-b.BETY)/b.BETY)**2)/len(b.BETY))
    #peakx=max(abs((a.BETX-b.BETX)/b.BETX))
    #peaky=max(abs((a.BETY-b.BETY)/b.BETY))
    #rmsdx=sqrt(sum(((a.DX/sqrt(a.BETX)-b.DX/sqrt(b.BETX)))**2)/len(b.DX))
    #peakdx=max(abs((a.DX/sqrt(a.BETX)-b.DX/sqrt(b.BETX))))
    rmsx=0
    rmsy=0
    peakx=0
    peaky=0
    rmsdx=0
    peakdx=0
    phasexlist=MakePairs(x, b)
    phaseylist=MakePairs(y, b)

    dphixa=[];dphiya=[];dphixb=[];dphiyb=[]
    for ph in phasexlist:
        [bpm1,bpm2]=ph.upper().split()
        dphixa.append(x.PHASE[x.indx[bpm1]])
	dphixb.append(b.MUX[b.indx[bpm2]]-b.MUX[b.indx[bpm1]])

    for ph in phaseylist:
        [bpm1,bpm2]=ph.upper().split()
        dphiya.append(y.PHASE[y.indx[bpm1]])
	dphiyb.append(b.MUY[b.indx[bpm2]]-b.MUY[b.indx[bpm1]])
    dphixa=array(dphixa);dphiya=array(dphiya)
    dphixb=array(dphixb);dphiyb=array(dphiyb)
    if len(dphixa)>0:
	    rmsphix=sqrt(sum((dphixa-dphixb)**2)/len(dphixa))
	    peakphix=max(abs(dphixa-dphixb));
    else:
	    rmsphix=0
	    peakphix=0
	    print "Warning: No horizontal BPMs matching when computing beta-beating"
    if len(dphiya)>0:
	    rmsphiy=sqrt(sum((dphiya-dphiyb)**2)/len(dphiya))
	    peakphiy=max(abs(dphiya-dphiyb))
    else:
	    rmsphiy=0
	    peakphiy=0
	    print "Warning: No vertical BPMs matching when computing beta-beating"



    #return array([rmsx,rmsy, peakx, peaky])
    return array([rmsx,rmsy, peakx,peaky, rmsphix,rmsphiy,
                  peakphix,peakphiy, rmsdx, peakdx])




##########################################
def correctbeat(a,beat_input, cut=0.01, app=0, path="./"):
################################################
	R=   transpose(beat_input.sensitivity_matrix)
	vector=beat_input.computevector(a)
	wg=beat_input.wg
	weisvec=array(concatenate([sqrt(wg[0])*ones(len(beat_input.phasexlist)),
		sqrt(wg[1])*ones(len(beat_input.phaseylist)),
                sqrt(wg[2])*ones(len(beat_input.betaxlist)),
                sqrt(wg[3])*ones(len(beat_input.betaylist)),
		sqrt(wg[4])*ones(len(beat_input.displist)),
                sqrt(wg[5])*ones(2)]))
	Rnew=transpose(transpose(R)*weisvec)
	delta=-matrixmultiply(generalized_inverse(Rnew,cut),(vector-beat_input.zerovector)/beat_input.normvector)
	writeparams(delta, beat_input.varslist, app,  path=path)
	return


###################################################
def correctbeatEXP(x,y,dx,beat_input, cut=0.01, app=0, path="./"):
###########################################################
	R=   transpose(beat_input.sensitivity_matrix)
        vector=beat_input.computevectorEXP(x,y,dx)
        wg=beat_input.wg
	weisvec=array(concatenate([sqrt(wg[0])*ones(len(beat_input.phasexlist)),
		sqrt(wg[1])*ones(len(beat_input.phaseylist)),
                sqrt(wg[2])*ones(len(beat_input.betaxlist)),
	        sqrt(wg[3])*ones(len(beat_input.betaylist)),
		sqrt(wg[4])*ones(len(beat_input.displist)),
		sqrt(wg[5])*ones(2)]))
	Rnew=transpose(transpose(R)*weisvec)
	#print (vector-beat_input.zerovector) # Good line to debug big correction strs
	delta=-matrixmultiply(generalized_inverse(Rnew,cut),(vector-beat_input.zerovector)/beat_input.normvector)
        writeparams(delta, beat_input.varslist, app, path)
	return [delta, beat_input.varslist]





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


	################################
	def computevectorEXP(self,x,y,dx):
	###############################################################################
		phix=[]
		phiy=[]
		betx=[]
		bety=[]
		disp=[]
		print "Exp tunes in beta_input class: ", x.Q1, y.Q2
		for ph in self.phasexlist:
			[bpm1,bpm2]=ph.split()
			phix.append(x.PHASEX[x.indx[bpm1]])

		for ph in self.phaseylist:
                        [bpm1,bpm2]=ph.split()
                        phiy.append(y.PHASEY[y.indx[bpm1]])
		for b in self.displist:
			disp.append(dx.NDX[dx.indx[b]])
		return array(concatenate([phix,phiy,betx,bety,disp,[x.Q1],[y.Q2]]))



	#####################
	def computevector(self,x):
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




