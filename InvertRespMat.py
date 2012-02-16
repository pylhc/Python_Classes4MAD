#!/usr/bin/env python

from Numeric import *
from string import *
from MLab import mean
#from LinearAlgebra import singular_value_decomposition
from LinearAlgebra import *
import sys
from random import gauss, seed
from FFT import real_fft
import pickle
from random import random
def initMatrix(filename):
    global matrix
    matrix = []
    s=[]
    for line in open(filename):
        s=[]
        dat=split(line)
        for j in range(0,len(dat)):
            s.append(float(dat[j])+random()*1e-11*0.0)
        matrix.append(s)
    matrix=array(matrix)


initMatrix('f1001re')
matrix_f1001re=matrix
initMatrix('f1001im')
matrix_f1001im=matrix
initMatrix('f1010re')
matrix_f1010re=matrix
initMatrix('f1010im')
matrix_f1010im=matrix
initMatrix('dispy')
matrix_dispy=matrix

initMatrix('errarray')
matrix_err=ravel(matrix) 

#print matrix_err

#print matrix_f1010im

matrix_all=concatenate([matrix_f1001re,matrix_f1001im,matrix_f1010re,matrix_f1010im,matrix_dispy])

vals1=singular_value_decomposition(matrix_f1001re)[1]
vals2=singular_value_decomposition(matrix_f1001im)[1]
vals3=singular_value_decomposition(matrix_f1010re)[1]
#vals4=singular_value_decomposition(matrix_f1010im)[1]
vals5=singular_value_decomposition(matrix_dispy)[1]
vals6=singular_value_decomposition(matrix_all)[1]

vals4=zeros(len(vals3))

#generalized_inverse(matrix)

#inv_x=generalized_inverse(matrix_x, 0.01)
#inv_y=generalized_inverse(matrix_y, 0.002)

# Dumping the inverse matrix in pickle format
# HIGHEST PROTOL is used with -1
#pickle.dump(inv_x,open('Inv_x','w'),-1)
#pickle.dump(inv_y,open('Inv_y','w'),-1)
#pickle.dump(matrix_x,open('Mat_x','w'),-1)
#pickle.dump(matrix_y,open('Mat_y','w'),-1)


f = open('singvals.dat','w')
for i in range(len(vals1)):
    f.write(str(vals1[i])+' '+str(vals2[i])+' '+\
            str(vals3[i])+' '+str(vals4[i])+' '+\
            str(vals5[i])+' '+str(vals6[i])+'\n')
f.close()

stre=-matrixmultiply(generalized_inverse(matrix_all,0.04),matrix_err)

f = open('SkewsStr.madx','w')
for i in range(len(stre)):
    f.write('str'+str(i+1)+'= str'+str(i+1)+' + ('+str(stre[i])+');'+'\n')
f.close()


#print matrixmultiply(generalized_inverse(matrix_all,0.001),matrix_err)

