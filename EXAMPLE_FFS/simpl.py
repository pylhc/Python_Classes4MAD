import math
import copy
import sys
from os import system
from simplex import Simplex
from mapclass import *

#########################
def writeparams(deltafamilie):
#########################
    "Writting a change of parameters of MADX variables in relative units"
    global variables
    g = open ('changeparameters', 'w')
    i=0
    for var in variables: 
        g.write(var+' = '+var+'*(1+('+str(deltafamilie[i])+'));\n')
        i +=1
    g.close()
    return



#########################
def chi2(deltafamilies):
#########################
    "Provide a chi2"
    global sigmaFFS
    print deltafamilies
    writeparams(deltafamilies)
    system('madxdev < ff.madx > scum')
    map=Map()
    chi=map.sigma('x',sigmaFFS)*weis[0]+map.sigma('y',sigmaFFS)*weis[1]
    print chi
    return  chi



##############
def main():
##############
    global deltaall,weis, variables, sigmaFFS
    #Definition of lattice parameters and sigmas at FFS start
    betx=64.9999
    bety=17.9997
    gamma=3e6
    ex=68e-8
    ey=1e-8    
    sigmaFFS=[sqrt(ex*betx/gamma), sqrt(ex/betx/gamma), sqrt(ey*bety/gamma), sqrt(ey/bety/gamma), 0.01]
    
    print 'Sigmas at FFS start:', sigmaFFSstart
    #Definition of MADX variables to vary (sextupole str)
    variables=['ksf6','ksf5', 'ksd4', 'ksf1', 'ksd0']
    #Definition of weights:
    weis=[10,10000]
    #Starting point for the MADX variables
    deltafamilies1=[-0.0,0.0,-0.00,-0.00,0.00]
    #Starting increment of the MADX variables.
    #Changes are relative in this example, see writeparams
    incr=[0.05,0.05,-0.05,-0.01,0.01]
    s = Simplex(chi2, deltafamilies1, incr)
    values, err, iter = s.minimize(epsilon =6e-18)
    print 'args = ', values
    print 'error = ', err
    print 'iterations = ', iter

if __name__ == '__main__':
    main()
