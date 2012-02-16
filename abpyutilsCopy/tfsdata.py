import os,gzip
from numpy import array
from numpy.core.records import fromarrays

class tfsdata(object):
  def __init__(self,param={},data=[],labels=[]):
    self.param=param
    self.data=data
    self.labels=labels
  def __getitem__(self,k):
#    if hasattr(k,'__iter__'):
    if type(k) is tuple:
      if len(k)==1:
        return self.data[k[0]]
      else:
        return self.data[k[0]][self.labels[k[1]]]
    else:
      return self.data[k]
  def __setitem__(self,k,v):
#    if hasattr(k,'__iter__'):
    if type(k) is tuple:
      if len(k)==1:
        self.data[k[0]]=v
      else:
        self.data[k[0]][k[1]]=v
    else:
      self.data[k]=v
  def __len__(self):
    return len(self.data)

  def fromfile(self,filename):
    self.data={}
    self.param={}
    if hasattr(filename,'__iter__'):
      file=filename
    else:
      if filename[-3:]=='.gz':
        file=gzip.open(filename)
      else:
        try:
          file=open(filename)
        except IOError:
          file=gzip.open(filename+'.gz')
    for line in file:
      if line.strip():
        f=line.split()
        if (f[0] == '@'):  # descriptor lines
          try:
            self.param[f[1]]=conv(f[2],f[3])
          except:
            print "bad descriptor"," ".join(f)
        elif ( f[0] == '*'): # self.labels lines
          f.pop(0)
          f=[pythonname(l) for l in f]
          self.labels=f
          for l in self.labels: self.data[l]=[]
        elif (f[0] == '$'):  # type lines
          f.pop(0) ; self.types=f
        elif (f[0].startswith('#')):  # comment lines
          pass
        else :   # data lines
          f=map(conv,self.types,f)
          for l in self.labels:
            d=f.pop(0)
            self.data[l].append(d)
    data=[self.data[n] for n in self.labels]
    self.data=fromarrays(data,names=','.join(self.labels))
    return self

  def fromstring(self,string):
    self.fromfile(string.split('\n'))
    return self

def conv(t,i) :
  if ('e' in t): i=float(i)
  if ( ('s' in t) and (i[0]=='"') ): i=i[1:-1]
  if ('d' in t): i=int(i)
  return i

def pythonname(string):
  string=string.replace('[','')
  string=string.replace(']','')
  string=string.replace('.','_')
  string=string.replace('$','_')
  return string.lower()


def test():
  self=tfsdata()
#  self.fromfile('../twiss.lhcb1.data.gz')
  return True




