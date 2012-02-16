#--- set of python routines:
#--- 1. best corrector for beta-beat (July 28, 2008)
#--- 2. Best N Corrector ~ like micado
#--- 3. changeparams contains all corrs, even if zero value
#--- Both Numeric + Numpy should work

import sys,os,time,datetime,string
ver=sys.version; print "Using Python",string.split(ver)[0]
#try:
#    import numpy as npy; default=0
#    from GenMatrix25 import *
#    from metaclass25 import twiss
#except ImportError:
#    try:
#        import Numeric as npy; default=1
#        from GenMatrix import *
#        from metaclass import twiss
#        import MLab
#    except ImportError: sys.exit()


#import Numeric as npy
try:
	from metaclass import *
except:
	from metaclass25 import *
try:
	from Numeric import *
except:
	from numpy import *


default=1
from GenMatrix import *
#from metaclass import twiss
from AllLists import *


def wrtparOLD(dfam,app=0,path="./"):
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    a = datetime.datetime.fromtimestamp(time.time())
    #dfam=[(k,v) for (k,v) in dfam.items()]
    g = open (path+'changeparameters', mode);i=0
    f = open (path+'changeparameters.tfs', mode)
    print >>f, "@", "APP", "%le", app
    print >>f, "@", "PATH","%s", path
    print >>f, "@", "DATE", "%s", a.ctime()
    print >>f, "*", "NAME    ", "    DELTA"
    print >>f, "$", "%s    ",       "    %le"
    if default==0: items=sorted(dfam.iteritems(),key=lambda (k,v):\
                                (npy.abs(v),k),reverse=True)
    if default==1: items=[(k,v) for (k, v) in dfam.items()]
    for k,v in items:
        g.write('%12s = %12s'%(k,k)+' +  '+\
                '( %e )'%(v)+';\n')
        f.write('%12s   %e'%(k, v)+'\n'); i+=1
    g.close();f.close()
    return

def wrtpar(corr, dfam,app=0,path="./"):
    if (app == 0):mode='w'
    if (app == 1):mode='a'
    a = datetime.datetime.fromtimestamp(time.time())
    g = open (path+'changeparameters', mode);i=0
    f = open (path+'changeparameters.tfs', mode)
    #--- note minus sign for change.tfs to do corr
    print >>f, "@", "APP", "%le", app
    print >>f, "@", "PATH","%s", path
    print >>f, "@", "DATE", "%s", a.ctime()
    print >>f, "*", "NAME    ", "    DELTA"
    print >>f, "$", "%s    ",       "    %le"
    for k in corr:
        if k in dfam.keys():
            g.write('%12s = %12s'%(k,k)+' +  '+\
                    '( %e )'%(dfam[k])+';\n')
            f.write('%12s   %e'%(k, -dfam[k])+'\n'); i+=1
        else:
            #g.write('%12s = %12s'%(k,k)+' +  '+\
            #        '( %e )'%(0.0)+';\n')
            f.write('%12s   %e'%(k, 0.0)+'\n'); i+=1
    g.close();f.close()
    return

def calcRMS(R):
    rms=npy.sqrt(npy.sum(R**2)/len(R))
    ptp=npy.max(R)-npy.min(R)
    return rms, ptp

def calcRMSNumeric(R):
    rms=sqrt(sum(R**2)/len(R))
    ptp=max(R)-min(R)
    return rms, ptp

def sortDict(adict):
    list=sorted(adict.items(), key=lambda (k,v): (v,k))
    for item in list: print item[0],":",item[1]
    
def bCorr(X,Y,DX, beat_input, cut=0.001,app=0,path="./"):
    R=transpose(beat_input.sensitivity_matrix)
    b=beat_input.computevectorEXP(X,Y,DX)-beat_input.zerovector
    corr=beat_input.varslist
    m,n=npy.shape(R)
    if len(b)==m and len(corr)==n:
        rms,ptop=calcRMS(b); inva=npy.linalg.pinv(R)
        print "initial {RMS, Peak}: {", '%e' %(rms),',', \
              '%e' %(ptop), "} mm"
        print "finding best over",n,"correctors"
        for i in range(n):
            dStren=npy.dot(inva[i,:],b)
            bvec=b-npy.dot(R[:,i],dStren)
            rm,ptp=calcRMS(bvec)
            if rm < rms: rms=rm; rbest=i; rStren=dStren
            if ptp < ptop: ptop=ptp; pbest=i; pStren=dStren
        print "final {RMS, Peak}: {", '%e' %(rms),',', \
              '%e' %(ptop), "}"
        if rbest==pbest:
            print "best corr:", corr[rbest], '%e'%(rStren), "1/m"
        else:
            print "--- warning: best corr for rms & peak are not same"
            print "RMS best corr:", corr[rbest], '%e'%(rStren), "1/m"
            print "Peak best corr:", corr[pbest], '%e'%(pStren), "1/m"
    else: print "dimensional mismatch in input variables"
    return

def bCorrNumeric(X,Y,DX, beat_input, cut=0.001,app=0,path="./"):
    R=npy.transpose(beat_input.sensitivity_matrix)
    b=beat_input.computevectorEXP(X,Y,DX)-beat_input.zerovector
    corr=beat_input.varslist; m,n=npy.shape(R)
    if len(b)==m and len(corr)==n:
        rms,ptop=calcRMSNumeric(b); inva=generalized_inverse(R,cut)
        print "initial {RMS, Peak}: {", '%e' %(rms),',', \
              '%e' %(ptop), "} mm"
        print "finding best over",n,"correctors"
        for i in range(n):
            dStren=npy.matrixmultiply(inva[i,:],b)
            bvec=b-npy.matrixmultiply(R[:,i],dStren)
            rm,ptp=calcRMSNumeric(bvec)
            if rm < rms: rms=rm; rbest=i; rStren=dStren
            if ptp < ptop: ptop=ptp; pbest=i; pStren=dStren
        print "final {RMS, Peak}: {", '%e' %(rms),',', \
              '%e' %(ptop), "}"
        if rbest==pbest:
            print "best corr:", corr[rbest], '%e'%(rStren), "1/m"
        else:
            print "--- warning: best corr for rms & peak are not same"
            print "RMS best corr:", corr[rbest], '%e'%(rStren), "1/m"
            print "Peak best corr:", corr[pbest], '%e'%(pStren), "1/m"
    else: print "dimensional mismatch in input variables"
    return

def bNCorr(X,Y,DX, beat_input, cut=0.001,ncorr=3, app=0, tol=1e-9,path="./"):
    R=transpose(beat_input.sensitivity_matrix);m,n=npy.shape(R)
    b=beat_input.computevectorEXP(X,Y,DX)-beat_input.zerovector
    corr=beat_input.varslist;
    inva=npy.linalg.pinv(R,cut);RHO2={};rmss,ptopp=calcRMS(b)
    for ITER in range(ncorr):
        for j in range(n):
            if j not in RHO2:
                RHO=[k for k in RHO2.keys()];RHO.append(j)
                RR=npy.take(R,RHO,1);invaa=npy.take(inva,RHO,0)
                dStren=npy.dot(invaa,b)
                bvec=b-npy.dot(RR,dStren)
                rm,ptp=calcRMS(bvec)
                if rm < rmss:
                    rmss=rm; ptopp=ptp; rbest=j; rStren=dStren
        print 'ITER:', ITER+1, '  RMS,PK2PK:',rmss,ptopp
        RHO2[rbest]=(rmss,ptopp)
        if (rm < tol):
            print "RMS converged with",ITER, "correctors";break
        if (rm < rmss):
            print "stopped after",ITER, "correctors";break
    itr=0;RHO3={} #-- make dict RHO3={corr:strength,...}
    for j in RHO2: RHO3[corr[j]]=rStren[itr];itr+=1
    print '\n',sortDict(RHO3);wrtpar(corr,RHO3,app,path)
    return RHO3

def itrSVD(X,Y,DX,beat_input, cut=0.001, iter=1, app=0, tol=1e-9,path="./"):
    R=transpose(beat_input.sensitivity_matrix);m,n=npy.shape(R)
    b=beat_input.computevectorEXP(X,Y,DX)-beat_input.zerovector
    corr=beat_input.varslist;inva= generalized_inverse(R,cut);
    rmss,ptopp=calcRMSNumeric(b);RHO2={};dStren=zeros(len(corr))
    print 'Initial Phase-Beat:', '{RMS,PK2PK}',rmss,ptopp
    for ITER in range(iter):
        dStren=dStren+npy.matrixmultiply(inva,b)
        bvec=b-npy.matrixmultiply(R,dStren)
        rm,ptp=calcRMSNumeric(bvec); print 'ITER',iter, '{RMS,PK2PK}',rm,ptp
    for j in range(len(corr)): RHO2[corr[j]]=dStren[j]
    wrtpar(corr,RHO2,app,path)

def bNCorrNumeric(X,Y,DX, beat_input, cut=0.001, ncorr=3, app=0, tol=1e-9,path="./"):
    R=transpose(beat_input.sensitivity_matrix);m,n=npy.shape(R)
    b=beat_input.computevectorEXP(X,Y,DX)-beat_input.zerovector
    corr=beat_input.varslist;
    # "-" sign added, to be verified by Ram (1-Aug-2008)
    #-- sign removed by Ram, look comment below in loop(Nov 27, 2009)
    inva= generalized_inverse(R,cut);  
    RHO2={};rmss,ptopp=calcRMSNumeric(b)
    print 'Initial Phase-beat {RMS,PK2PK}:', rmss,ptopp
    for ITER in range(ncorr):
        for j in range(n):
            if j not in RHO2:
                RHO=[k for k in RHO2.keys()];RHO.append(j)
                RR=npy.take(R,RHO,1);invaa=npy.take(inva,RHO,0)
                dStren=npy.matrixmultiply(invaa,b)
                #--- calculate residual due to 1..nth corrector
                bvec=b-npy.matrixmultiply(RR,dStren)
                rm,ptp=calcRMSNumeric(bvec)
                if rm < rmss:
                    rmss=rm; ptopp=ptp; rbest=j; rStren=dStren
        print 'ITER:', ITER+1, '  RMS,PK2PK:',rmss,ptopp
        RHO2[rbest]=(rmss,ptopp)
        if (rm < tol):
            print "RMS converged with",ITER, "correctors";break
        if (rm < rmss):
            print "stopped after",ITER, "correctors";break
    itr=0;RHO3={} #-- make dict RHO3={corr:strength,...}
    for j in RHO2: RHO3[corr[j]]=rStren[itr];itr+=1
    print '\n',sortDict(RHO3); wrtpar(corr,RHO3,app,path)
    return RHO3

def bNCorrNumericSim(a, beat_input, cut=0.1,ncorr=3, app=0, tol=1e-9,path="./"):
    R=transpose(beat_input.sensitivity_matrix);m,n=npy.shape(R)
    b=beat_input.computevector(a)-beat_input.zerovector
    corr=beat_input.varslist;
    inva=generalized_inverse(R,cut);RHO2={};rmss,ptopp=calcRMSNumeric(b)
    for ITER in range(ncorr):
        for j in range(n):
            if j not in RHO2:
                RHO=[k for k in RHO2.keys()];RHO.append(j)
                RR=npy.take(R,RHO,1);invaa=npy.take(inva,RHO,0)
                dStren=npy.matrixmultiply(invaa,b)
                bvec=b-npy.matrixmultiply(RR,dStren)
                rm,ptp=calcRMSNumeric(bvec)
                if rm < rmss:
                    rmss=rm; ptopp=ptp; rbest=j; rStren=dStren
        print 'ITER:', ITER+1, '  RMS,PK2PK:',rmss,ptopp
        RHO2[rbest]=(rmss,ptopp)
        if (rm < tol):
            print "RMS converged with",ITER, "correctors";break
        if (rm < rmss):
            print "stopped after",ITER, "correctors";break
    itr=0;RHO3={} #-- make dict RHO3={corr:strength,...}
    for j in RHO2: RHO3[corr[j]]=rStren[itr];itr+=1
    print RHO3, wrtpar(corr,RHO3,app,path)
    return RHO3


'''
if __name__ == "__main__":
    mpath='/afs/cern.ch/user/r/rtomas/w1/MAD/LHC/Beta-Beat/NewRespMatrixGen/'
    mpath='/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/MODEL/LHCB1/nominal.opt/'
    print "loading Full Response optics"
    if default==0:
        FullResponse=pickle.load(open(mpath+'/FullResponse','r'))
    else: FullResponse=pickle.load(open(mpath+'/FullResponse.Numeric','r'))

    #x=twiss('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Correction/getphasex.out')
    #y=twiss('/afs/cern.ch/eng/sl/lintrack/Beta-Beat.src/Correction/getphasey.out')
    #--- LHC B1 data Nov 26 overnight
    x=twiss('./1:21:20/getphasex.out')
    y=twiss('./1:21:20/getphasey.out')

    #varslist=["b2", "b3", "b4", "b5", "b6", "b8", "b10", "b11", "b12", \
    #          "b13", "b14", "b15", "b16", "b17", "b18", "b19", "b20", \
    #          "b21", "b22", "b23", "b24", "b25", "b26", "b27", "b28", "\
    #          b29", "b30", "b31", "b32", "b33", "b34", "b35", "b36",\
    #          "b37", "b38", "b39", "b40", "b41", "b42", "b43", "b44", "\
    #          b45", "b46", "b47", "b48", "b49", "b50", "b51", "b52", "b53",\
    #          "b54", "b55", "b56", "b57", "b58", "b59", "b60", "b61", "b62",\
    #          "b63", "b64", "b65", "b66", "b67", "b68","b69", "b70", "b71",\
    #          "b72", "b73", "b74", "b75", "b76", "b77", "b78", "b79", "b80",\
    #          "b81", "b82", "b83", "b84", "b85", "b86", "b87", "b88", "b89",\
    #          "b90", "b91", "b92", "b93","b94", "b95", "b96", "b97", "b99", \
    #          "b101", "b102", "b103", "b104", "b105", "b106", "b107", "b108",\
    #          "b109"]

    varslist=quadvarsb1();variables=varslist
    
    phasexlist=MakePairs(x, FullResponse['0'],modelcut=0.02, errorcut=0.013)
    phaseylist=MakePairs(y, FullResponse['0'],modelcut=0.02, errorcut=0.013)
    betaxlist=[]; betaylist=betaxlist; dx=[]
    displist=MakeList(dx, FullResponse['0'])

    #x=twiss('twiss.dat'); y=twiss('twiss.base')
    #phasexlist=[];phaseylist=[];betaxlist=BPMsb1()
    #betaylist=betaxlist; displist=BPMsb1();
    wei=[1,1,1,1,1,10]
    beat_inp=beat_input(varslist, phasexlist, phaseylist, \
                        betaxlist, betaylist, displist, wei)
    sensitivity_matrix=beat_inp.computeSensitivityMatrix(FullResponse)
    #os.system('madx < job.lhc.madx > scum')
    #y=twiss('twiss.base');x=twiss('twiss.dat')
    #print "Initial Beta-beat:"
    #print betabeat(x,y)


    #------ Apply Correction
    if default==0:
        #correctbeat(x, beat_inp, cut=0.04, app=0)
        #bCorr(x,y,dx beat_inp, app=0,path)
        quads=bNCorr(x,y,dx,beat_inp, cut=0.1,ncorr=5, app=0, tol=1e-4,path='./')
    elif default==1:
        #bCorrNumeric(x,y,dx,beat_inp, cut=0.01,app=0,path='./')
       # bNCorrNumeric(x,y,dx,beat_inp, cut=0.1,ncorr=10, app=0,tol=1e-4,path='./')
       itrSVD(x,y,dx,beat_inp, cut=0.1,iter=1, app=0,tol=1e-4,path='./')
    else: print "No numpy or Numeric modules"
    
    sys.exit()
'''
