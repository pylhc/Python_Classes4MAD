#!/afs/cern.ch/eng/sl/lintrack/Python-2.5/bin//python2.5

from string import *
from numpy import *
from sys import *
import os
import re
from Sussix import *




###########################
def gettune(filename):
###########################
    x=[]
    y=[]
    for line in open(filename ):
        x.append(float(split(line)[1]))
        y.append(float(split(line)[3]))

    x=array(x)
    y=array(y)



    sussix_inp(ir=1, turns=512, tunex=0.127, tuney=0.177, istun=0.03, idam=2, narm=300)
    a=sussix(x,y)
    

    return a



#print "your file is: ",argv[1]

a=gettune(argv[1])

#Frequency and amplitude of lines:

outf=open(argv[1]+'.sussix','w')

for i in range(len(a.ox)):
    outf.write( str(a.ox[i])+' '+str(a.ax[i])+' '+str(a.oy[i])+' '+str(a.ay[i])+'\n')

outf.close()
#Tunes are here
#print a.tunexy
