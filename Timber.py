import os, sys, time
from numpy import *

def parseout(file,sep=','):
  data={}
  dataon=False
  header=True
  log=[]
  datavars=[]
  fh=open(file,"r")
  for l in fh.readlines():
    if l.startswith('VARIABLE'):
      vname=l.split()[1]
      print 'Found var %s' % vname
      if vname in data:
        t,v=data[vname]
      else:
        t,v=[],[]
        data[vname]=(t,v)
        datavars.append(vname)
      dataon=False
      header=False
    elif l.startswith('Timestamp'):
      dataon=True
    elif l=='\n':
      dataon=False
    elif 'null' in l:
      dataon=False
    elif dataon:
      ll=l.strip().split(sep)
      t.append(ll[0])
      v.append(ll[1:])    # This is always a vector
    elif header:
      log.append(l)
  if len(log)>0:
    data['log']=log
  data['datavars']=datavars
  fh.close()
  return data



def toseconds(s):
  try:
    miliseconds=float(s.split(".")[1])
  except:
    print "Error with ",s, "in toseconds()"
    sys.exit()
  tup=time.strptime(s.split(".")[0], "%Y-%m-%d %H:%M:%S")
  return time.mktime(tup)+miliseconds/1000.


def extract(t,d,t1,t2):
    
    t1s=toseconds(t1)
    t2s=toseconds(t2)
    if t2s<t1s:
        print "Bad range in extract(), t2<t1"
        sys.exit()
    rt,rv=[],[]
    for i in range(len(d)):
        ts=toseconds(t[i])
        if ts>=t1s and ts<=t2s:
            rt.append(t[i])
            try:
              rv.append(float(d[i][0]))   # [0] is because data is always a vector, but for us it will be scalar
            except:
              print "Error with ",d[i][0], "in extract"
        if ts>t2 and len(rv)>0:
            return average(rv), std(rv)
    return average(rv), std(rv)





def extractnew(t,d,t1,t2):    
    t1s=toseconds(t1)
    t2s=toseconds(t2)
    if t2s<t1s:
        print "Bad range in extract(), t2<t1"
        sys.exit()
    rt,rv=[],[]
    for i in range(len(d)):
        ts=toseconds(t[i])
        if ts>=t1s and ts<=t2s:
            rt.append(t[i])
            try:
              rv.append(float(d[i][0]))   # [0] is because data is always a vector, but for us it will be scalar
            except:
              print "Error with ",d[i][0], "in extract"
        if ts>t2 and len(rv)>0:
            return average(rv), std(rv), rv
    return average(rv), std(rv), rv




        
def avepoints(data,t1,t2,file):
    vars= data['datavars']
    N=len(vars)
    f=open(file,"w")
    print >>f, "* NAME   AVE  ERR"
    print >>f, "$ %s   %le %le" 
    for i in range(N):
        t=data[vars[i]][0]
        d=data[vars[i]][1]
        if len(t)> 0 and len(d)>0:
          ave, rms= extract(t, d, t1,t2)
          print >>f,vars[i],  ave, rms
        else:
          print "**** No data in ", vars[i]
    return ave, rms
                   
