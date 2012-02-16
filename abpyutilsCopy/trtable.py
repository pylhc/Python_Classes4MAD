try:
  import transport as _transport
except ImportError:
  import os
  os.system('(cd %s; rm transport.so; make transport.so)' %  os.path.dirname(__file__))
  import transport as _transport


_cols=_transport.colname().split()
_sindex=_cols.index('s')
_lindex=len(_cols)-_sindex

from table import table
from plotfunc import plotfunc
from numpy import zeros

class trtable(table,plotfunc):
  def __init__(self,elems,cols=_cols,rowattr=True):
    table.__init__(self,elems,_cols,rowattr=rowattr)
    (self.br,self.gr,self.betx,self.bety)=(1,7.E6/938,1,1)
  def track(self,stop=None):
    if stop:
      stop=stop._row_index[0]
    else:
      start=self._row_index[0]
    _transport.loop(self._data[start:stop].T,_sindex,_lindex)
#    if hasattr(self,'currplot'):
#      self.currplot.update()
    return self
  def getindex(self,elems,cols):
    indices=[]
    for i in elems:
      for j in cols:
        indices.append(eval('self.%s._addr[%d]' % (i,j)))
    return indices
  def put(self,indices,x):
    self._data[zip(*indices)]=x
  def take(self,indices):
    return self._data[zip(*indices)]

  def __mul__(self,n):
    olddata=self._data
    oldnames=self._row_names
    r,c=olddata.shape
    newdata=zeros((r*n,c),dtype=float)
    newnames=zeros(r*n,dtype=oldnames.dtype)
    for i in range(n):
      newdata[i::n]=olddata
      newnames[i::n]=oldnames
    t=trtable.from_row_names(newnames,self._col_names)
    t.copydata(newdata)
    t.l=t.l/n
    t.kn1l=t.kn1l/n
    t.ks1l=t.ks1l/n
    t.kn0l=t.kn0l/n
    t.ks0l=t.ks0l/n
    return t

def from_tfstable(t,rowattr=True,kicker=False):
  names=t._row_names.tolist() + ['end']
  l=trtable(names,rowattr=rowattr)
  if kicker:
    begcolm ='l k1l angle hkick vkick'.split()
    begcol  ='l kn1l kn0l kn0l ks0l'.split()
  else:
    begcolm ='l k1l angle'.split()
    begcol  ='l kn1l kn0l'.split()
  endcol='s betx bety alfx alfy dx dpx dy dpx x px y py mux muy'.split()
  for n,o in zip(begcol,begcolm):
    if hasattr(t,o):
      getattr(l,n)[:-1]+=getattr(t,o)
      getattr(l,n)[0]+=getattr(t,o)[0]
  for n in endcol:
    if hasattr(t,n):
      getattr(l,n)[1:]=getattr(t,n)
      getattr(l,n)[0]=getattr(t,n)[0]
  return l

def reverse_table(t):
  names=t._elems[:]
  names.pop()
  names.reverse()
  names.append('stop')
  tn=trtable(names)
  for i in range(len(tn)-1):
    tn[i,_sindex:]=t[::-1][i,_sindex:]
  for i in range(len(tn)-1):
    tn[i,:_sindex]=t[::-1][i+1,:_sindex]
  tn.alfx=-tn.alfx
  tn.alfy=-tn.alfy
  tn.dpx=-tn.dpx
  tn.dpy=-tn.dpy
  tn.px=-tn.px
  tn.py=-tn.py
  return tn


_sindex=7
if __name__=='__main__':
  from scipy.optimize import fsolve
  s=trtable('d0 q1 d1 q2 d2 q3 d3 q4 d4 end')
  s.d0.l, s.d1.l, s.d2.l, s.d3.l, s.d4.l = 23, 2.7, 1, 2.9, 10
  s.q1.l, s.q2.l = 9.2, 7.8
  s.q1.kn1l = 122./23350*s.q1.l
  def cons():
    s.q3.l, s.q4.l = s.q2.l, s.q1.l
    s.q2.kn1l = -s.q1.kn1l/s.q1.l*s.q2.l
    s.q3.kn1l, s.q4.kn1l = s.q2.kn1l, s.q1.kn1l

  s.d0.betx,s.d0.bety=0.25,0.25
  cons()
  s.track()
#  s.plot('betx bety')
  sm=(s*10)
  sm.track()
#  sm.plot('betx bety')
  from tfstable import tfstable
  from utils import *
  o=tfstable('twiss.ir5b1.data')
  c=from_tfstable(o)
  t=c.copy()
  t.track()
#  o.show(o.s<1,'^s$ ^l$ betx dx')
#  c.show(c.s<1,'^s$ ^l$ betx dx')
  print 'trtable.py: test OK'
  o=tfstable('twiss.lhcb1.data',rowattr=False)
  c=from_tfstable(o,rowattr=False,kicker=True)
  t=c.copy(rowattr=False)
  t.track()


