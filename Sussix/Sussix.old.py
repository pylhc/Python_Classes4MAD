import sussix4drivexxNoO
from Numeric import *
from random import gauss
from string import *

#######################################################
def sussix_inp(ir=1, turns=500, tunex=0.1, tuney=0.2, istun=0.05, idam=2):
#######################################################
    
    si=open("sussix_v4.inp", "w")
    si.write('C\nC INPUT FOR SUSSIX_V4 ---17/09/1997---\n')
    si.write('C DETAILS ARE IN THE MAIN PROGRAM SUSSIX_V4.F\nC\n\n')
    si.write('ISIX  = 0\nNTOT  = 1\nIANA  = 1\nICONV = 0\n')
    si.write('TURNS = 1 '+str(turns)+'\n')
    si.write('NARM  = 2\nISTUN = 1 '+str(istun)+'  '+str(istun)+'\n')
    si.write('TUNES = '+str(tunex)+'  '+str(tuney)+' .07\n')
    si.write('NSUS  = 0\nIDAM  = '+ str(idam)+ '\n')
    si.write('NTWIX = 1\nIR    = '+str(ir) + '\nIMETH = 2\nNRC   = 6\nEPS   = 1D-3\n')
    si.write('NLINE = 0\nL,M,K = \nIDAMX = 1\nNFIN  = 500\nISME  = 1\n')
    si.write('IUSME = 200\nINV   = 0\nIINV  = 250\nICF   = 0\nIICF  = 350\n')
    si.close()     
    return



#####################
class sussixout:
#####################
    "Structure for sussix output"
    def __init__(self, z):
        self.tunexy=z[0]
        self.amplitude=z[1]
        self.phase=z[2]
        self.ox=z[3]
        self.ax=z[4]
        self.oy=z[5]
        self.ay=z[6]

        


#####################
def fillorcut(vec, n):
####################
    lv=len(vec)
    if (lv>n):
        return vec[0:n]
    elif   (lv<n):
        return concatenate([vec, zeros(n-lv)])
    elif   (lv==n):
        return vec


################################################
def sussix(x0,x1):
################################################
    dim=(4404, 2, 14,14,300,300,300,300)
    y0=fillorcut(x0, 1100)
    y1=fillorcut(x1, 1100)
    y3=fillorcut(concatenate([y0,y1]), dim[0])
    z=sussixout(sussix4drivexxNoO.sussix4drivenoise(y3))
    return z



#x=array(map(cos ,  0.35*2*pi*array(range(0,100))))
#y=array(map(cos ,  0.14*2*pi*array(range(0,100))))


#sussix_inp(ir=1, turns=100, tunex=0.35, tuney=0.14, istun=0.05, idam=2)
#a=sussix(x,y)
#print a.tunexy
#print a.amplitude
#for i in range(100):
#    print a.ox[i], a.ax[i],a.oy[i], a.ay[i]


