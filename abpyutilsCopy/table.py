from numpy import array,zeros,where
from view import view
from utils import numtostr, opticsutils

import re

def _count(d,k,i):
  d.setdefault(k, []).append(i)

def _mystr(d,nd):
  if hasattr(d,'__coerce__'):
    return numtostr(d)
  else:
    return ('%%%ds' % nd) % d

def _mknames(elems):
  row_names=[]
  for n in elems:
    if hasattr(n,'__iter__'):
      na,tbl=n
      for tbln in tbl._row_names:
        ntbln='%s.%s' % (na,tbln)
        row_names.append(ntbln)
    else:
      row_names.append(n)
  return row_names

def _mkind(elems,row_index):
  i=-1
  naind={}
  elind={}
  tblelems={}
  bknaind={}
  bktblelems={}
  for n in elems:
    if hasattr(n,'__iter__'):
      na,tbl=n
      _count(elind,tbl,row_index[i]+1) # for copying data if needed
      for j in tbl._row_names:
        i+=1
        _count(naind,na,row_index[i]) # for complex child
      for j in tbl._elems:
        _count(tblelems,na,j) # for complex child
    else:
      i+=1
      _count(naind,n,row_index[i]) # for making simple child
  return naind,elind,tblelems,bknaind,bktblelems

def _mkdata(row_names,col_names):
  data=zeros((len(row_names),len(col_names)))
  return data

def _mkcopydata(data,elind,col_names):
  for tbl in elind:
    ncol=[col_names.index(k) for k in tbl._col_names if k in col_names]
    vcol=[tbl._col_names.index(k) for k in tbl._col_names if k in col_names]
    for i in elind[tbl]:
      data[i:i+len(tbl._row_names),ncol]=tbl.all[:,vcol]

def _mkrowindex(row_names,col_names):
  row_index=range(len(row_names))
  addr=dict([j,(slice(None),i)] for i,j in enumerate(col_names))
  addr['all']=slice(None)
  return row_index,addr

def _mkaddr(row_index,col_names):
  addr=dict([j,(row_index,i)] for i,j in enumerate(col_names))
  addr['all']=(row_index,slice(None))
  return addr

def _rowattr(self,elems,tblelems,col_names,data,naind):
  for n in elems:
    if hasattr(n,'__iter__'):
      na,tbl=n
      ctable=table(tblelems[na],col_names,
          data=data,row_index=naind[na],rowattr=True)
      setattr(self,na,ctable)
    else:
      ind=naind[n]
      ctable=table([n]*len(ind),col_names,
          data=data,row_index=ind,rowattr=False)
      setattr(self,n,ctable)


class table(view,opticsutils):
  _show_char=12
  _show_len=5
  def __init__(self,elems,col_names,data=None,row_index=None,rowattr=True):
    if not hasattr(elems,'__iter__'):  elems=elems.split()
    if not hasattr(col_names,'__iter__'): col_names=col_names.split()
    row_names=_mknames(elems)
    if data is None:
      data=_mkdata(row_names,col_names)
    if row_index is None:
      row_index,addr=_mkrowindex(row_names,col_names)
    else:
      addr=_mkaddr(row_index,col_names)
    naind,elind,tblelems,bknaind,bktblelems=_mkind(elems,row_index)
    _mkcopydata(data,elind,col_names)
    view.__init__(self,data,**addr) #must exist before setting address
    if rowattr:
      _rowattr(self,elems,tblelems,col_names,data,naind)
    self._row_names=array(row_names)
    self._row_index=array(row_index)
    self._row_iref=dict(zip(row_index,row_names))
    self._row_ref=dict(zip(row_names,row_index))
    self._col_names=col_names
    self._elems=elems
    self._col_ref=dict([(n,i) for i,n in enumerate(col_names) ])
    self._col_show={'_name':15}

  def __getitem__(self,index):
    return self._data[index]

  def __setitem__(self,index,value):
    self._data[index]=value

  def __len__(self):
    return len(self._row_names)

  def __str__(self):
    out=[]
    colsn=self._col_names
    cols=[self._col_names.index(n) for n in colsn]
    row=[ ('%%-%ds' % self._col_show.get(j,self._show_char))%j for j in ['_name']+colsn]
    out.append(' '.join(row))
    c1,c2=0,len(self._row_index)
    for i,m in zip(self._row_index,self._row_names):
      c1+=1
      if c1<self._show_len or c2-c1<self._show_len:
        row=[('%%-%ds' % self._col_show['_name']) % m]
        for j,n in zip(cols,colsn):
          d=self._data[i,j]
          nd=self._col_show.get(n,self._show_char)
          row.append(_mystr(d,nd))
        out.append(' '.join(row))
      elif c1==self._show_len: out.append('   ...')
    return '\n'.join(out)

  __repr__=__str__

  def cpattern(self,regexp):
    regexp=regexp.split()
    out=[]
    for r in regexp:
       out+=[ n for n in self._col_names if n==r]
    return out

  def pattern(self,regexp):
    out=zeros(len(self._data),dtype=bool)
    c=re.compile(regexp)
    for n,i in zip(self._row_names,self._row_index):
      if c.search(n): out[i]=True
    return out
  __floordiv__=pattern

  def where(self,ind):
    return where(ind)[0]

  def show(self,rows='.*',cols='.*'):
    if isinstance(rows,str):
      rows=self.pattern(rows)
    if isinstance(cols,str):
      colsn=self.cpattern(cols)
      cols=[getattr(self,n) for n in colsn]
    else:
      colsn=['']*len(cols)
    irows=where(rows)[0]
    out=[]
    row=[ ('%%-%ds' % self._col_show.get(j,self._show_char))%j for j in ['_name']+colsn]
    out.append(' '.join(row))
    for i in irows:
      row=[ ('%%-%ds' % self._col_show['_name'])% self._row_iref[i] ]
      for n,c in zip(colsn,cols):
        d=c[i]
        nd=self._col_show.get(n,self._show_char)
        row.append(_mystr(d,nd))
      out.append(' '.join(row))
    return '\n'.join(out)


  def copydata(self,newdata):
    self._data[:]=newdata

  def copy(self,rowattr=True):
    t=self.__class__.from_row_names(self._row_names,self._col_names,rowattr=rowattr)
    t.copydata(self._data)
    return t

  def from_row_names(cls,row_names,col_names,rowattr=True):
    tmp={}
    out=[]
    for i in row_names:
      if '.' in i:
        n,r=i.split('.',1)
        if n in tmp:
          tmp[n][1].append(r)
        else:
          l=[n,[r]]
          out.append(l)
          tmp[n]=l
      else:
        out.append(i)
    for i in tmp:
      tmp[i][1]=from_row_names(cls,tmp[i][1][:],col_names,rowattr=rowattr)
    t=cls(out,col_names,rowattr=rowattr)
    return t
  from_row_names=classmethod(from_row_names)


if __name__=='__main__':
  print "Test table"
  out=[]
  t=table('qd d qf d','x y z')
  t.x=4
  t.x[2]=5
  t.d.y=4
  t.qf.z=3
  out.append(t._data[0,0]==4)
  out.append(t._data[2,0]==5)
  out.append(t._data[1,1]==4)
  out.append(t._data[3,1]==4)
  out.append(t._data[2,2]==3)
  tt=table(['start',('f1',t),('f2',t),('f1',t),'end'],'x z')
  tt.f1.z=1
  out.append(tt._data[1,1]==1)
  ttt=table(['start',('f1',tt),('f2',t),('f1',t),'end'],'x z a')
  ttt.f1.f1.a=.5
  out.append(ttt._data[2,2]==0.5)
  ttt.show('\.d','x z')
  print ttt.f2
  res=reduce(lambda x,y: x and y,out)
  if res:
    print 'table.py: test OK'
  else:
    print 'table.py: test FAILED'

