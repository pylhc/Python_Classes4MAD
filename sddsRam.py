#-- usage:
# a=sddsRead('filname')
# print a.data -- to see a dict with all data from sdds file

try: import numpy as npy
except ImportError:
    try: import Numeric as npy
    except ImportError: sys.exit()
    
from string import split,replace
import sys,os
import gzip


class sddsRead:   
    def __init__(self, filename):
        name=[];self.data={}
        iterPar=0;iterArr=0
        
        if '.gz' in filename:
            f=gzip.open(filename, 'rb')
        else:
            f=open(filename, 'r')
                                            
        
        line=f.readlines()
        for j in range(1,len(line)):
            if "&parameter" in line[j]:
                iterPar+=1
                name.append(replz(split(line[j])[1]))
            if "&array" in line[j]:
                iterArr+=1
                name.append(replz(split(line[j])[1]))
            if "!page number" in line[j]: break
        line=line[j+1:]
        for k in range(iterPar):
            
            self.data[name[k]]=replz(line[k])
        line=line[k+1:];name=name[k+1:];j=0
        for k in range(iterArr):
            arrSize=int(line[j]);j=j+1
            #DATA=readLine(line[j]);j=j+1
            DATA=split(line[j]);j=j+1

            #--- need this if newlines in the middle of arrays
            while len(DATA)<arrSize:
                #DATA=npy.concatenate((DATA,readLine(line[j])));j=j+1
                DATA=DATA+split(line[j]); j=j+1
            self.data[name[k]]=DATA

def readLine(lin):
    return npy.fromstring(lin,dtype=int,sep=' ')[:-1]
    
def replz(par):
    par=replace(par,'name=','')
    par=replace(par,',','')
    par=replace(par,'"','')
    return par

    
def writeCO2TFS(ab,dirTFS='./', corr='None', Q1=0.0,Q2=0.0):
    HPOS=split(ab.data['HMON-POS_B1'])
    VPOS=split(ab.data['VMON-POS_B1'])
    HPOSRMS=split(ab.data['HMON-POS-RMS_B1'])
    VPOSRMS=split(ab.data['VMON-POS-RMS_B1'])
    HNAMES=split(ab.data['HMON-NAME_B1'])
    VNAMES=split(ab.data['VMON-NAME_B1'])
    
    if len(HNAMES)!=len(HPOS):
        print 'unequal arrays, continuing w/o writing TFS file'
    else:
        f=open(dirTFS+'CO.TFSX','w')
        f.write('@ DATE %s '+ab.data['DATE'])
        f.write('@ TIME %s '+ab.data['TIME'])
        f.write('@ Q1 %le '+str(Q1)+'\n')
        f.write('@ EXCITER %le '+corr+'\n')
        if 'H' in corr: f.write('@ DIR %le H\n')
        if 'V' in corr: f.write('@ DIR %le V\n')
        f.write("*"+'%10s %10s %10s' % ("NAME", "CO", "CORMS")+"\n")
        f.write("$"+'%10s %10s %10s' % ("%s", "%le", "%le")+"\n")
        for j in range(len(HNAMES)):
            f.write('%10s %10s %10s' %(HNAMES[j], HPOS[j], HPOSRMS[j])+'\n')
        f.close()
    if len(VNAMES)!=len(VPOS):
        print 'unequal arrays, continuing w/o writing TFS file'
    else:
        f=open(dirTFS+'CO.TFSY','w')
        f.write('@ DATE %s '+ab.data['DATE'])
        f.write('@ TIME %s '+ab.data['TIME'])
        f.write('@ Q2 %le '+str(Q2)+'\n')
        f.write('@ EXCITER %le '+corr+'\n')
        if 'H' in corr: f.write('@ DIR %le H\n')
        if 'V' in corr: f.write('@ DIR %le V\n')

        f.write("*"+'%10s %10s %10s' % ("NAME", "CO", "CORMS")+"\n")
        f.write("$"+'%10s %10s %10s' % ("%s", "%le", "%le")+"\n")
        
        for j in range(len(VNAMES)):
            f.write('%10s %10s %10s' %(VNAMES[j], VPOS[j], VPOSRMS[j])+'\n')
        f.close()
    return


if __name__ == "__main__":
    dir='/afs/cern.ch/user/g/grumolo/scratch1/beta.beat.jun2.08/'
    file='tttt.sdds'
    a=sddsRead(dir+file)
    writeCO2TFS(a,dir,Q1=0.0, Q2=0.0)
    
#    import pylab
#    fileName=sys.argv[1]
#    a=sddsRead(fileName)
#    print a.data.keys()
    
#    ##--- multiq example
#    fac=1.0 #-- to avoid overflow error in pylab
#    for ii in a.data:
#        if 'rawData' in ii:
#            #indx=int(replace(ii,'fftData-H-',''))
#            #xdata=npy.arange(len(a.data[ii]))/1024.
#           #ydata=npy.ones(len(a.data[ii]))*a.data['acqTimes'][indx-1]
#            pylab.plot(a.data[ii])
#            #pylab.scatter(xdata[:512],ydata[:512],c=a.data[ii][:512]*fac,\
#            #              faceted=False,cmap=pylab.cm.jet,s=10)
#            pylab.show()
            
