import re
from numpy import *
from  gzip import GzipFile
from utils import *

class madtable:
  """
  from pymad import *
  t=madtable('twiss.lhcb1.x1y5.data')
  """
  lookup={}
  def __init__(self, filename=None, model=None, gen=None):
    self.rows=0
    self.cols=[]
    self.labels=[]
    self.descs={}
    self.idx={}
    self.data={}

    if filename:
      self.fromfile(filename)
    elif model:
      self.frommodel()
      self.model=model
      self.gen=gen

  def fromfile(self,filename):

    # dispatch convert filename in list of lines
    self.filename=filename
    if type(filename) is str:
      if len(filename)<200:
        try:
          file=open(filename)
          if filename[-3:]=='.gz':
            file=GzipFile(fileobj=file)
        except IOError:
          try:
            file=open(filename+'.gz')
            file=GzipFile(fileobj=file)
          except IOError:
            file=[]
      else:
        file=filename.split('\n')
    elif type(filename) is list:
      file=filename
    else:
      file=[]

    if file:
      file=[line.split() for line in file]

      while 1:
        f=file.pop(0)
        if (f[0] == '@'):
          self.descs[f[1]]=typeconv(f[2]," ".join(f[3:]))
        elif ( f[0] == '*'): # self.labels lines
          f.pop(0)
          f=map(pythonname,f)
          self.labels=f
        elif (f[0] == '$'):  # type lines
          f.pop(0)
          self.types=f
          break
        elif (f[0].startswith('#')):  # comment lines
          pass

      file=[ tuple(map(typeconv,self.types,f)) for f in file ]
      data=rec.fromrecords(file,names=self.labels)
      self.idx=dict(enumerate(data['NAME']))
      for l in self.labels:
        self.data[l]=data[l]
        setattr(self,l,self.data[l])
        setattr(self,l.lower(),self.data[l])
        self.__class__.lookup[id(self.data[l])]=l

  def add_column(self,d,label="USER"):
    self.data[label]=d
    temp=d[0]
    try :
      d[0]='ds'
      self.types.append('%s')
      d[0]=temp
    except TypeError:
      self.types.append('%le')
    self.labels.append(label)
    self.data[label]=array(self.data[label])
    setattr(self,label.lower(),self.data[label])
    setattr(self,label,self.data[label])
    self.__class__.lookup[id(self.data[label])]=label

  def update(self):
    t=madtable(self.filename)
    self=t


  def pattern(self,regexp,*a) :
    """
    t.name[t.pattern("MQ")]
    t.name[t.pattern("^IP1$")]
    t.name[t.pattern(1.387,t.s)]
    """
    if len(a)==1 :
      col=a[0]
    else :
      col=self.data[self.labels[0]]

    regexp=str(regexp)
    c=re.compile(regexp, re.IGNORECASE)
    idx=slist()
    for i in xrange(0,len(col)) :
      c.search(str(col[i])) and idx.append(i)
    return idx

  def name2idx(self,name) :
    """
    t.name[t.name2idx(["MQXC.3L5.B2","MQXB.B2L5.B2"])]
    """
    idx=[ self.idx[str(s)] for s in name]
    return slist(idx)

  def elem(self,ra,*a) :
    """
    t.s[t.elem("MQXC.3L5.B2")]
    t.name[t.elem(483.937625784,t.s)]
    """
    if len(a)==1 :
      col=a[0]
    else :
      col=self.data[self.labels[0]]

    ra=str(ra).lower()

    idx=slist()
    for i in xrange(0,len(col)) :
      if (str(col[i]).lower() == ra) :
        idx.append(i)
    return idx


  def range(self,ra,rb,*a) :
    """
    t.name[t.range("MQT.13L5.B2","MS.13L5.B2")]
    """
    if len(a)==1 :
      col=a[0]
    else :
      col=self.data[self.labels[0]]

    ra=str(ra).lower()
    rb=str(rb).lower()

    c=re.compile(ra, re.IGNORECASE)
    d=re.compile(rb, re.IGNORECASE)
    try:
      ida=self.pattern(ra)[0]
      idb=self.pattern(rb)[0]
      idx=range(ida,idb+1)
    except IndexError:
      idx=slist()
    return slist(idx)

  def nrange(self,sa,sb,*a) :
    """
    t.name[t.nrange(1,7)]
    t.name[t.nrange(100,200,t.betx)]
    """
    if len(a)==1 :
      col=a[0]
    else :
      try :
        col=self.s
      except KeyError:
        col=self.data[self.labels[2]]

    idx=slist()
    for i in xrange(0,len(col)) :
      if (col[i] >= sa) and (col[i] <= sb) :
        idx.append(i)
    return idx


  def cpattern(self,regexp) :
    """
    t.cpattern("BET")
    array(t.cpattern("BET"))[:,t.nrange(0,7)]
    """
    c=re.compile(regexp, re.IGNORECASE)
    col = [getattr(self,l) for l in self.labels if  c.search(l) ]
    return slist(col)

  def lpattern(self,regexp) :
    """
    t.lpattern("BET")
    """
    c=re.compile(regexp, re.IGNORECASE)
    lab = [l for l in self.labels if  c.search(l) ]
    return slist(lab)

  def lab2col(self,lab) :
    """
    t.lab2col(t.lpattern("BET"))
    """
    col=[getattr(self,l) for l in lab]
    return slist(col)

  def select(self,row=".*",col=[0]) :
    """
    t.select("MQ","N BE ")
    t.select(["MQT.13L5.B2", "MCBH.13L5.B2"],["betx", "bety"])
    t.select([1,3],[t.betx, t.bety])
    """

    idx=slist()
    try :
      if ( type(row) is type("") ) :
        for a in row.split() :
          idx.extend(self.pattern(a))
      elif isinstance(row,list):
        if row :
          if type(row[0]) is type("") :
            idx=self.name2idx(row)
          elif type(row[0]) is type(0) :
            idx=row
    except :
      print "Error in row selection:", sys.exc_info()[0]
      idx=slist()

    cols=slist()
    try :
      if   (type(col) is type("")) :
        for a in col.split() :
          cols.extend(self.cpattern(a))
      elif isinstance(col,list) :
        for i in col:
          if isinstance(i,str):
            cols.extend(self.lab2col([i]))
          elif type(i) is type( array([]) ):
            cols.append(i)
          elif isinstance(i,int):
            cols.extend(self.lab2col([self.labels[i]]))
    except :
      print "Error in col selection:", sys.exc_info()[0]
      cols=slist()

    labels=[]
    for i in cols:
      try :
         labels.append(self.__class__.lookup[id(i)])
      except KeyError:
        labels.append("USER")
    return idx,cols,labels

  def notzero(self,row=".*",col=[0]) :
    """
    t1.notzero("[LR]5\.","k1l")
    t1.show(t1.notzero(t1.range("s.ds.l5.b1","s.ds.r5.b1"),"k1l"),"NAME ^L$ K1L")

    """
    a=(row,col)
    ind,cols,lbl=self.select(*a)
#    print ind,cols,lbl
    a=array(cols).transpose()
    idx=slist()
    for i in ind:
      if all(abs(a[i,:])>1E-15) :
        idx.append(i)
    return idx

  def extract(self,row=".*",col=[0]) :
    """
    t.extract("MQ","BET")
    t.extract("MQ",["NAME","S"])
    t.extract("MQ",[t.name,t.s])
    """
    a=(row,col)
    ind,cols,lbl=self.select(*a)
    a=array(cols).transpose()
    a=a[ind,:]
    return a

  def show(self,row=".*",col=[0],a3=None) :
    """
    t.show("MQ","NAME ^S^ BET")
    t.show("MQ","NAME ^S^ BET","out.dat")
    """
    if a3:
      f=open(a3,"w")
    else:
      f=None
    a=(row,col)
    ind,cols,labels=self.select(*a)

    if labels:
      s="%-"+str(21)+"s "
      s=s*len(cols)
      a=array(cols,dtype='|S20')
      print >> f, s % tuple(labels)
      for i in ind:
        print >> f, s % tuple(a[:,i])
    else :
      print "No column chosen"

    f and f.close()
    return

  def export(self,filename):
    f=open(filename,'w')
    keys=self.descs.keys()
    keys.remove('NAME')
    keys.insert(0,'NAME')
    for i in keys:
      print >>f, '@ %-16s %%%is "%s"' % (i,len(self.descs[i]),self.descs[i])
    fmt="%-16s "* len(self.labels)
    print >>f, "*", fmt % tuple(self.labels)
    print >>f, "$", fmt % tuple(self.types)
    for i in range(self.lines+1):
      print >>f, fmt % tuple(map(lambda n:self.data[n][i],self.labels))



class slist(list):

  def __or__(self,other):
    u={}
    for i in self+other: u[i]=True
    out=u.keys()
    out.sort()
    return slist(out)

  def __and__(self,other):
    u={}
    for i in self :
      try :
        other.index(i)
        u[i]=True
      except ValueError:
        pass
    out=u.keys()
    out.sort()
    return slist(out)

  def __sub__(self,other):
    u={}
    for i in self :
      try :
        other.index(i)
      except ValueError:
        u[i]=True
    out=u.keys()
    out.sort()
    return slist(out)


