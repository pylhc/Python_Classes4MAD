import re
import sys
import os
import pexpect
from math import acos, asin, atan, atan2, ceil, cos, cosh, degrees, e, exp, fabs, floor, fmod, frexp, hypot, ldexp, log, log10, modf, pi, pow, radians, sin, sinh, sqrt, tan, tanh
from pylab import *
from numpy import *
import readline
from user import home
import gzip
#from string import split

class madtable:
  """
  from pymad import *
  t=madtable('twiss.lhcb1.x1y5.data')
  """
  dict={}
  def __init__(self, filename, regexp=None,istring=None):
    self.data={}
    self.descs={}
    self.idx={}
    self.lines=-1
    if (regexp): c=re.compile(regexp, re.IGNORECASE)
    if istring:
      file=filename.split('\n')
    elif filename[-3:]=='.gz':
      file=gzip.open(filename)
      self.filename=filename
    else:
      file=open(filename)
      self.filename=filename
    for line in file:
      if line.strip():
        f=line.split()
        if (f[0] == '@'):  # descriptor lines
          try:
            self.descs[f[1]]=self.conv(f[2],f[3])
          except:
            print "bad descriptor"," ".join(f)
        elif ( f[0] == '*'): # self.labels lines
          f.pop(0)
          f=[l.replace('[','') for l in f]
          f=[l.replace(']','') for l in f]
          f=map(pythonname,f)
          self.labels=f
          for l in self.labels: self.data[l]=[]
        elif (f[0] == '$'):  # type lines
          f.pop(0) ; self.types=f
        elif (f[0].startswith('#')):  # comment lines
          pass
        else :   # data lines
          if (regexp):
            if (c.search(line)): save=True
            else : save=False
          else : save=True
          if (save):
            self.lines=self.lines+1
            f=map(self.conv,self.types,f)
    #        print self.labels
            for l in self.labels:
              d=f.pop(0)
              self.data[l].append(d)
            try:
              self.idx[self.data['NAME'][-1]]=self.lines
            except KeyError:
              self.idx[self.lines]=self.lines
    # End of file
#    for l in self.descs.keys():
#      setattr(self,l.lower(),self.descs[l])
#      setattr(self,l,self.descs[l])
    for l in self.labels:
      self.data[l]=array(self.data[l])
      setattr(self,l.lower(),self.data[l])
      setattr(self,l,self.data[l])
      self.__class__.dict[id(self.data[l])]=l

  def addcol(self,d,label="USER"):
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
    self.__class__.dict[id(self.data[label])]=label

  def conv(self,t,i) :
    if ('e' in t): i=float(i)
    if ( ('s' in t) and (i[0]=='"') ): i=i[1:-1]
    if ('d' in t): i=int(i)
    return i


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
      idbb=self.pattern(rb)
      idb=0
      for i in idbb:
        if i>=ida:
          idb=i
          break
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
         labels.append(self.__class__.dict[id(i)])
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

class envelope:
  """
  find aperture given a survey and a twiss table, example:
  import pymad as m
  s1=m.madtable('survey.lhcb1.data')
  s2=m.madtable('survey.lhcb2.data')
  t1=m.madtable('twiss.lhcb1.data')
  t2=m.madtable('twiss.lhcb2.data')
  ap1=m.envelope(s1,t1)
  ap2=m.envelope(s1,t1)

  hold(False)
  figure(figsize=(6,6))
  hold(True)
  plot(ap1.co[:,2],ap1.co[:,0])
  plot(ap1.xp[:,2],ap1.xp[:,0],'g',linewidth=.1)
  plot(ap1.yp2D[:,2],ap1.yp2D[:,0],'r',linewidth=.1)
  plot(ap1.xm[:,2],ap1.xm[:,0],'g',linewidth=.1)
  plot(ap1.ym2D[:,2],ap1.ym2D[:,0],'r',linewidth=.1)
  axis([-10000,10000,-15000,5000])

  t1.select(t1.pattern("IP1"),t1.name,t1.s,ap1.xsize,ap1.ysize)
  savefig('ring.eps',dpi=600)
  """

#    kbeta=1.1,      # beta beating
#    nsigma=9.5,      # 9.5
#    emit=3.75E-6,   #3.75E-6
#    delta=1.129E-4, # RMS energy spread
#    tol=4.6E-3,  # CO=3mm + dtol=1.6mm
#    deltamax=8E-4,   # for chromaticity measurment
#    betamaxarc=180,  # maximum beta in the arcs
#    dxmaxarc=2,  # maximum beta in the arc
  def __init__(self,s,t, kbeta=1.1, nsigma=9.5, nemit=3.75E-6, delta=1.129E-4, tol=4.6E-3, deltamax=8E-4, betamaxarc=180, dxmaxarc=2, gamma=7000):
    self.co=zeros([len(s.x),3],float)
    self.xp=zeros([len(s.x),3],float)
    self.xm=zeros([len(s.x),3],float)
    self.yp=zeros([len(s.x),3],float)
    self.ym=zeros([len(s.x),3],float)
    self.yp2D=zeros([len(s.x),3],float)
    self.ym2D=zeros([len(s.x),3],float)
    self.xsize=zeros(len(s.x),float)
    self.ysize=zeros(len(s.x),float)

    for i in xrange(len(s.x)):
      vro = array([s.x[i],s.y[i],s.z[i]])
      theta,phi,psi = s.theta[i],s.phi[i],s.psi[i]
      betx,bety,dx,dy,x,y= t.betx[i],t.bety[i],t.dx[i],t.dy[i],t.x[i],t.y[i]
      thetam=array([[cos(theta) ,           0,sin(theta)],
           [          0,           1,         0],
           [-sin(theta),           0,cos(theta)]])
      phim=  array([[          1,          0,          0],
          [          0,cos(phi)   ,   sin(phi)],
          [          0,-sin(phi)  ,   cos(phi)]])
      psim=  array([[   cos(psi),  -sin(psi),          0],
          [   sin(psi),   cos(psi),          0],
          [          0,          0,          1]])
      wm=dot(thetam,dot(phim,psim))
      ex=dot(wm,array([1,0,0]))
      ey=dot(wm,array([0,1,0]))
      self.co[i]=vro+x * ex + y * ey
      emit=nemit/gamma
      dx+= dxmaxarc*sqrt(betx/betamaxarc)
      dy+= dxmaxarc*sqrt(bety/betamaxarc)
      xsize=kbeta* (nsigma*sqrt(betx*emit + (dx*delta)**2) + deltamax*dx)+ tol
      ysize=kbeta* (nsigma*sqrt(bety*emit + (dy*delta)**2) + deltamax*dx)+ tol
      self.xp[i]=self.co[i] + xsize * ex
      self.xm[i]=self.co[i] - xsize * ex
      self.yp[i]=self.co[i] + ysize * ey
      self.ym[i]=self.co[i] - ysize * ey
      self.yp2D[i]=self.co[i] + ysize * ex
      self.ym2D[i]=self.co[i] - ysize * ex
      self.xsize[i]=xsize
      self.ysize[i]=ysize

def madname(string):
  string=string.replace('_','.')
  string=string.replace('underscore','_')
  return string

def pythonname(string):
  string=string.replace('_','underscore')
  string=string.replace('.','_')
  return string


class madinput:
  """
  import pymad as m
  reload(m); i=m.madinput("dipfirst.thin.seq")
  """

  def __init__(self, filename,istring=None):
    self._v={}
    self._e={}
    self.varlist=[]
    self.elemlist=[]
    self.seqlist=[]
    if istring:
      s=filename
    elif filename[-3:0]=='.gz':
      f=gzip.open(filename)
      s=f.read()
    else:
      f=open(filename,"rU")
      s=f.read()
    self.decodestring(s)

  def decodestring(self,s):
    s=s.lower()
    s=s.replace("const ","")
    s=s.replace("shared ","")
    r=re.compile("!.*")
    s=r.sub("",s)
    s=s.replace('\n','')
    s=s.replace(' ','')
    r=re.compile("[^;]*;")
    #rvar=re.compile("([\w\.]+):?=([\w\.]+);")
    rvar=re.compile("(^[\w\.]+)(:?=)([^;:,]+)")
    relem=re.compile("(^[\w\.]+):([\w\.]+)(,.+)?;")
    relem2=re.compile("(^[\w\.]+)(,.+);")
    rdot1=re.compile("([a-z][a-z0-9]+)\.([a-z0-9])")
    rdot2=re.compile("([a-z])\.([a-z])")


    currseq=None
    for i in r.finditer(s):
      st=i.group()
      st=st.replace("_","underscore")
      st=rdot1.sub(r"\1_\2",st)
      st=rdot2.sub(r"\1_\2",st)
      st=st.replace("..","__")
      st=st.replace("^","**")
      # variable parsing
      if (rvar.match(st)):
        tok=rvar.match(st).groups()
        # store var in the dictionary
        v=var(tok[2],self._v)
        v.name=tok[0]
        if not ':=' in st:
          v.actvalue=v.value()
        self._v[tok[0]]=v
        self.varlist.append(v)
      # elements parsing
      elif (relem.match(st) or relem2.match(st)):
#        print st
#        print currseq
        try :
          tok=relem.match(st).groups()
          el=elem(tok[0])
          try :
            el.parent=self._e[tok[1]]
          except KeyError:
            el.parent=tok[1]
          attr=tok[2]
        except AttributeError:
          tok=relem2.match(st).groups()
          try:
            el=self._e[tok[0]]
            attr=tok[1]
          except KeyError:
            print "No parent:" + tok[0]
            el=elem(tok[0])
        # store elem in the dictionary
        el._env=self._v
        self._e[ el.name ]=el
#        el.show()
        # attribute parsing
        if attr :
          attr=list(attr)
          ch=False
          for i in range(len(attr)):
            if attr[i]=='{':
              ch=True
            if attr[i]=='}':
              ch=False
            if (ch) and (attr[i]==',') :
              attr[i]=' '
          attrl=(''.join(attr[1:])).split(',')
          for a in attrl:
            # a attribute statement
            try :
              # match attribute definition
              stok=rvar.match(a).groups()
              av=stok[2]
              an=stok[0]
              at=stok[1]
              # list case
              if av[0]=='{':
                av=av[1:-1].split()
              # trasnform string in variables
              if not an in ['refpos','from','bv']:
                if type(av) is str:
                  av=var(av,self._v)
                  if at=='=':
                    av.actvalue=av.value()
                else:
                  # list case
                  for i in range(len(av)):
                    av[i]=var(av[i],self._v)
                    if at=='=':
                      av[i].actvalue=av[i].value()
              # manage sequence information
              if an=='at':
                el._at[currseq]=av
                el._vartype['_at']=at
              elif an=='from':
                el._from[currseq]=av
                el._vartype['_from']=at
              else:
                setattr(el,an,av)
              el._vartype[an]=at

        
            except AttributeError :
              pass
#        el.show()
        # sequence management
        if (el.parent=='sequence'):
#          raise
#        if (st.startswith('bpmsw')): raise
          currseq=el.name
          el.elements=[]
          self.seqlist.append(el)
        elif (currseq) :
          try:
            el.sequence.append(currseq)
          except AttributeError:
            el.sequence=[ currseq ]
          self._e[currseq].elements.append(el)
        else:
          self.elemlist.append(el)
#        el.show()
#        print
      # sequence parsing
      elif (st=='endsequence;'):
        currseq=None

    self.__dict__.update(self._e)
    self.__dict__.update(self._v)
#    print 'globals().update(i.e)'
#    print 'globals().update(i.v)'



  def writevar(self, filename):
    try:
      f=open(filename,"w")
    except:
      return
    init="""@ NAME             %06s "VARIABLES"
@ TYPE             %06s "VARIABLES"
@ TITLE            %16s "VARIABLES"
@ DATE  %08s    14/10/05
@ TIME  %08s    14:15:06
* NAME VALUE EXPR
$ %s           %le        %s
"""
    f.write(init)
    ns = [v for v in self.__dict__.keys() if isinstance(self.__dict__[v],var)]
    ns.sort()
    for n in ns:
      v=self.__dict__[n]
      f.write( "%-20s %22.15e \"%s\"\n" % (n,v.value(),v.expr) )
    f.close()

  def select(self,a='.*',env=None):
    if not (env): env=self.__dict__
    n=a.split()
    names=[]
    for i in n:
      c=re.compile(i, re.IGNORECASE)
      names.extend(filter(c.search,env.keys()))
    names.sort()
    return names

  def show(self,a,env=None):
    if not (env): env=self.__dict__
    for i in self.select(a,env):
      env[i].show(env=self._v)

  def vars(self,a='.*'):
    return [i for i in self.select(a,env=self._v)]

  def elems(self,a='.*'):
    return [i for i in self.select(a,env=self._e)]

  def tomad(self,file=None):
    try:
      for i in self.varlist:
        print i.tomad()
      for i in self.elemlist:
        print i.tomad()
      for i in self.seqlist:
        print i.tomad()
        for j in i.elements:
          print j.tomad(sequence=i.name)
        print 'endsequence;'
    except :
      print "Error in row selection:", sys.exc_info()[0]
      print i




class slist(list):

  def __or__(self,other):
    u={}
    for i in list(self)+list(other): u[i]=True
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

  def __div__(self,other):
    u={}
    for i in self :
      try :
        other.index(i)
      except ValueError:
        u[i]=True
    out=u.keys()
    out.sort()
    return slist(out)

  def __add__(self,n):
    out=slist(self)
    for i,j in enumerate(self) :
      out[i]+=n
    return out

  def __sub__(self,n):
    out=slist(self)
    for i,j in enumerate(self) :
      out[i]-=n
    return out


class var:
  def __init__(self,e,v):
    self.expr=e.replace(' ','')
    self._env=v
    self.actvalue=None
  def value(self,env=None):
    if self.actvalue:
      return  self.actvalue
    else:
      if not (env): env=self._env
      try:
        return eval(self.expr,globals(),env)
      except:
        return self.expr
  def __coerce__(self,other): return float(self),float(other)
  def __float__(self): return float(self.value())
  def __pos__(self): return float(self)
  def __neg__(self): return -float(self)
  def __str__(self):
    try:
      return "%e" % self.value()
    except NameError:
      (t,v,b)=sys.exc_info()
      return "expr='%s' Error: %s" % (self.name,self.expr,v.args[0])
  def __repr__(self):
    try:
      return "%e" % self.value()
    except NameError:
      (t,v,b)=sys.exc_info()
      return "expr='%s' Error: %s" % (self.name,self.expr,v.args[0])
  def show(self,env=None):
    try:
      v="%e" % self.value(env=env)
    except NameError:
      v=None
    if v:
      print "%s=%s value:%s" % (self.name,self.expr,v)

  def tomad(self):
    if self.actvalue==None:
      return madname("%s:=%s;" % (self.name,self.expr))
    else:
      return madname("%s=%s;" % (self.name,self.expr))

class elem:
  """
  mad elements
  """
  def __init__(self, name,**kwds):
    self.name=name
    self._vartype={'name': '=', 'parent': '=', 'sequence': '='}
    self._at={}
    self._from={}
    # store attribute via dictionary
    self.__dict__.update(kwds)

  def show(self,env=None):
    if not (env): env=self._env
    print "--- %s ---" % self.name
    for k in self.__dict__:
      if not k.startswith('_') and not k in ['elements']:
        a=getattr(self, k)
        if isinstance(a,str):
          print "%s%s%s" % (k,self._vartype[k],a)
        elif isinstance(a,var):
          print "%s%s%s value:%e" % (k,self._vartype[k],a.expr,a.value())
        elif isinstance(a,list):
          print "%s%s[" % (k,self._vartype[k]),
          for i in a:
            if isinstance(i,str):
              print "%s, " % (i),
            elif isinstance(i,var):
              print "%s value:%e, " % (i.expr,i.value()),
          print "]"
    if isinstance(self.parent,elem):
      self.parent.show()

  def __repr__(self):
    return  "<%s>" % self.name

  def madname(self):
    name=self.name.replace('_','.')
    name=name.replace('underscore','_')
    return name

  def tomad(self,sequence=None):
    if isinstance(self.parent,elem):
      out= "%s:%s" % (self.name,self.parent.name)
    elif isinstance(self.parent,str):
      out= "%s:%s" % (self.name,self.parent)
    att=self.__dict__.copy()
    del att['name']
    del att['parent']
    if 'sequence' in att: del att['sequence']
    if 'elements' in att: del att['elements']
    if sequence in self._at:
      out+=',%s%s%s' %('at',self._vartype['_at'],self._at[sequence])
    if sequence in self._from:
      out+=',%s%s%s' %('from',self._vartype['_from'],self._from[sequence])
    for k in att:
      if not k.startswith('_'):
        a=getattr(self, k)
        if isinstance(a,var):
          out+=',%s%s%s' %(k,self._vartype[k],a.expr)
        elif isinstance(a,str):
          out+=',%s%s%s' %(k,self._vartype[k],a)
        elif isinstance(a,list):
          out+=',%s%s{' % (k,self._vartype[k])
          for i in a:
            out+='%s,' %(i.expr)
          out+='}'
    out+=';'
    return madname(out)

class madrun:
  """
  r=m.madrun();
  r.s('a=4;value,a^2;')
  r.p('a=4;value,a^2;')
  """
  prompt=['\r\nX: ==>\r\n', pexpect.TIMEOUT]
  def __init__(self,prg='madx'):
    self.mad=pexpect.spawn(prg)
    self.mad.expect(self.prompt)
    print self.mad.before;
    self.mad.sendline('option,-echo;')
    self.mad.expect(self.prompt)


  def s(self,string):
    """
    Send command to madx and return string
    print r.s('call,file=job.sample.madx;')
    """
    self.mad.sendline(string)
    self.mad.expect(self.prompt)
    return self.mad.before

  def p(self,string):
    """
    Send command to madx and print string
    r.p('call,file=job.sample.madx;')
    """
    print self.s(string)

  def i(self):
    """
    Enter in madx command line
    r.i()
    """
    readline.parse_and_bind("tab: complete")
    fn1=home+'/.history_madx'
    fn2=home+'/.ipython/history'
    os.system("touch " + fn1)
    os.system("touch " + fn2)
    print "Type Ctrl-D to enter ipython"
    cmd=''
    readline.clear_history()
    readline.read_history_file(fn1)
    print 'X: ==> '
    try:
      while( 1):
        cmd=cmd+re.sub('\!.*','',raw_input())
        if cmd.endswith(';'):
          self.p(cmd)
          cmd=''
          print 'X: ==> '
    except EOFError:
      pass
    readline.write_history_file(fn1)
    readline.clear_history()
    readline.read_history_file(fn2)
    print "Type madx() to enter MadX"

def printmatch(r):
  """
  m.printmatch(t.select('MQY.4R5.B1','NAME BET[XY] ALF[XY] ^DX ^DPX'))
  """
  idx=r[0][0]
  n=r[1][0][idx]
  for i in range(1,len(r[2])):
    l=r[2][i]
    v=r[1][i][idx]
    print "constraint,sequence=lhcb1,range=%s,%s=%s;" % (n,l,v)

def plotenvelope(ip,t1,ap1,sele,eele,t2,ap2,sele2,eele2):
  """
m.plotenvelope("IP5",t1,ap1,"MQY.4L5.B1","MQY.4R5.B1",t2,ap2,"MQY.4L5.B2","MQY.4R5.B2")
  """

  yip5=ap1.xp[t1.elem(ip)[0],0]

  idxs1=t1.elem(sele)[0]
  idxe1=t1.elem(eele)[0]
  idxs2=t2.elem(sele2)[0]
  idxe2=t2.elem(eele2)[0]

  hold(True)
  title("Beam envelope")
  xlabel(r"$x [\rm{m}]$")
  ylabel(r"$y [\rm{m}]$")
  grid(True)

  # closed orbit
  x1=ap1.co[idxs1:idxe1,2]
  y1=ap1.co[idxs1:idxe1,0]-yip5
  plot(x1,y1,color=[0,0,1])
  x1=ap2.co[idxs2:idxe2,2]
  y1=ap2.co[idxs2:idxe2,0]-yip5

  # beam1
  plot(x1,y1,color=[1,0,0])
  x1=ap1.xp[idxs1:idxe1,2]
  y1=ap1.xp[idxs1:idxe1,0]-yip5
  x2=ap1.xm[idxs1:idxe1,2]
  y2=ap1.xm[idxs1:idxe1,0]-yip5
  x = concatenate( (x1,x2[::-1]) )
  y = concatenate( (y1,y2[::-1]) )
  fill(x, y, facecolor='b',alpha=0.2)

  # beam2
  x1=ap2.xp[idxs2:idxe2,2]
  y1=ap2.xp[idxs2:idxe2,0]-yip5
  x2=ap2.xm[idxs2:idxe2,2]
  y2=ap2.xm[idxs2:idxe2,0]-yip5
  x = concatenate( (x1,x2[::-1]) )
  y = concatenate( (y1,y2[::-1]) )
  fill(x, y, facecolor='r',alpha=0.2)



def plotenvelopey2D(ip,t1,ap1,sele,eele,t2,ap2,sele2,eele2):
  """
m.plotenvelopey2D("IP5",t1,ap1,"MQY.4L5.B1","MQY.4R5.B1",t2,ap2,"MQY.4L5.B2","MQY.4R5.B2")
  """

  yip5=ap1.xp[t1.elem(ip)[0],0]

  idxs1=t1.elem(sele)[0]
  idxe1=t1.elem(eele)[0]
  idxs2=t2.elem(sele2)[0]
  idxe2=t2.elem(eele2)[0]

  hold(True)
  title("Beam envelope")
  xlabel(r"$x [\rm{m}]$")
  ylabel(r"$y [\rm{m}]$")
  grid(True)

  # closed orbit
  x1=ap1.co[idxs1:idxe1,2]
  y1=ap1.co[idxs1:idxe1,0]-yip5
  plot(x1,y1,color=[0,0,1])
  x1=ap2.co[idxs2:idxe2,2]
  y1=ap2.co[idxs2:idxe2,0]-yip5

  # beam1
  plot(x1,y1,color=[1,0,0])
  x1=ap1.yp2D[idxs1:idxe1,2]
  y1=ap1.yp2D[idxs1:idxe1,0]-yip5
  x2=ap1.ym2D[idxs1:idxe1,2]
  y2=ap1.ym2D[idxs1:idxe1,0]-yip5
  x = concatenate( (x1,x2[::-1]) )
  y = concatenate( (y1,y2[::-1]) )
  fill(x, y, facecolor='b',alpha=0.2)

  # beam2
  x1=ap2.yp2D[idxs2:idxe2,2]
  y1=ap2.yp2D[idxs2:idxe2,0]-yip5
  x2=ap2.ym2D[idxs2:idxe2,2]
  y2=ap2.ym2D[idxs2:idxe2,0]-yip5
  x = concatenate( (x1,x2[::-1]) )
  y = concatenate( (y1,y2[::-1]) )
  fill(x, y, facecolor='r',alpha=0.2)


def pymadhelp():
  print "to start madx:            r=madrun()         "
  print "to start madxp:           r=madrun('madxp')  "
  print "to send command to mad:   print r.s('twiss;')"
  print "to enter in madx console: r.i()              "
  print "to exit:                  Ctrl-D             "

def plotbeta(t,thin=False):
  """
use after:
select,flag=twiss,column=name,s,l,betx,bety,dx,k0l,k1l;
use,sequence=lhcb1,range=IP5/e.DS.R5.B1;
twiss,betx=0.25,bety=0.25,file=twiss.data;

t=madtable('twiss.data')
plotbeta(t)
  """
  f=figure()
  elemAxes=subplot(111)
  subplots_adjust(right=0.75)
  if thin:
    mask=t.nrange(1E-5,1,abs(t.angle))
    tk0=t.extract(mask,[t.s,1E-2,t.angle,t.k1l])
    mask=t.nrange(1E-5,1,abs(t.k1l))
    tk1=t.extract(mask,[t.s,1E-2,t.k1l,t.k1l])
  else:
    mask=t.nrange(1E-5,1,abs(t.angle))
    tk0=t.extract(mask,[t.s-t.l,t.l,t.angle/t.l,t.k1l/t.l])
    mask=t.nrange(1E-5,1,abs(t.k1l))
    tk1=t.extract(mask,[t.s-t.l,t.l,t.k1l/t.l,t.k1l/t.l])

  k0p=bar(tk0[:,0],tk0[:,2],tk0[:,1],color='#aaffaa')
  setp(k0p,facecolor="#e0ffe0",edgecolor="#666666")
  k1p=bar(tk1[:,0],tk1[:,2],tk1[:,1],color='#aaaaff')
  setp(k1p,facecolor="#e0e0ff",edgecolor="#666666")
  setp( elemAxes.yaxis, visible=False)

  betaAxes=gcf().add_axes(elemAxes.get_position(), sharex=elemAxes, frameon=False, label="beta")
  bxp,=plot(t.s,sqrt(t.betx),'k',label=r'$\beta_x$')
  byp,=plot(t.s,sqrt(t.bety),'r',label=r'$\beta_y$')
  ylim(ymax=200)
  betaAxes.yaxis.set_label_position('left')
  betaAxes.yaxis.set_ticks_position('left')
  xlabel(r'$s [m]$')
  ylabel(r'$\sqrt\beta [m]$')

  dispAxes=gcf().add_axes(elemAxes.get_position(), sharex=betaAxes, frameon=False, label="disp")
  dispAxes.yaxis.tick_right()
  dispAxes.yaxis.set_label_position('right')
  dispAxes.yaxis.set_ticks_position('right')
  dxp,=plot(t.s,t.dx,'g',label=r'$D_x$')
  ylabel(r"$D [m]$")


  figlegend([bxp,byp,dxp,k0p[0],k1p[0]],[r'$\beta_x$',r'$\beta_y$',r'$D_x$','Bend','Quad'],'upper right')
  i=t.elem('IP5')
  draw()
  bs=text(t.s[i][0],1.8,'$\\beta^*=%6.4f$' % t.betx[i][0])
  return (t,bxp,byp,dxp,k0p,k1p,bs)


def plotbetaupdate(plot, cmd=';',thin=False,newfile=None):
  (t,bxp,byp,dxp,k0p,k1p,bs)=plot
#  r.p(cmd)
#  r.s('twiss,betx=0.25,bety=0.25,file=twiss.data;')
  if newfile:
    t=madtable(newfile)
  else:
    t=madtable(t.filename)
  bxp.set_ydata(sqrt(t.betx))
  byp.set_ydata(sqrt(t.bety))
  dxp.set_ydata(t.dx)
  bxp.set_xdata(t.s)
  byp.set_xdata(t.s)
  dxp.set_xdata(t.s)

  if thin:
    mask=t.nrange(1E-5,1,abs(t.angle))
    tk0=t.extract(mask,[t.s,1E-2,t.angle,t.k1l])
    mask=t.nrange(1E-5,1,abs(t.k1l))
    tk1=t.extract(mask,[t.s,1E-2,t.k1l,t.k1l])
  else:
    mask=t.nrange(1E-5,1,abs(t.angle))
    tk0=t.extract(mask,[t.s-t.l,t.l,t.angle/t.l,t.k1l/t.l])
    mask=t.nrange(1E-5,1,abs(t.k1l))
    tk1=t.extract(mask,[t.s-t.l,t.l,t.k1l/t.l,t.k1l/t.l])

  for i in range(len(k0p)):
    k0p[i].set_x(float(tk0[i,0]))
    k0p[i].set_width(float(tk0[i,1]))
    k0p[i].set_y(0)
    k0p[i].set_height(float(tk0[i,2]))
  for i in range(len(k1p)):
    k1p[i].set_x(float(tk1[i,0]))
    k1p[i].set_width(float(tk1[i,1]))
    k1p[i].set_y(0)
    k1p[i].set_height(float(tk1[i,2]))
  i=t.elem('IP5')
  setp(bs,'text','$\\beta^*=%g$' % t.betx[i][0])
  draw()
  return (t,bxp,byp,dxp,k0p,k1p,bs)




def plottune(s,lbl=''):
  title(r"$\rm{Tune versus} \delta$")
  xlabel("$\delta$")
  xlabel("Fractional tune")
  tt=r'$%s \rm{%s}$'
  plot(s.deltap,s.q1-s.q1.round(),label=tt %('Q_x',lbl))
  plot(s.deltap,s.q2-s.q2.round(),label=tt %('Q_y',lbl))
  text(0.0,0.31,r"$Q_x$")
  text(0.0,0.32,r"$Q_y$")
  grid(True)
  legend()


def plotbetabeat(t,t1,dp='0.0003'):
  title(r"$\rm{Beta beat: 1 - \beta(\delta=%s)/\beta(\delta=0)}$" % dp)
  ylabel(r"$\Delta\beta/\beta$")
  xlabel(r"$s$")
  plot(t.s,1-t1.betx/t.betx,label=r'$\Delta\beta_x/\beta_x$')
  plot(t.s,1-t1.bety/t.bety,label=r'$\Delta\beta_y/\beta_y$')
  grid(True)
  legend()

def plotw(t,lbl=''):
  title(r"Chromatic function: %s"%lbl)
#  ylabel(r"$w=(\Delta\beta/\beta)/\delta$")
  ylabel(r"$w$")
  xlabel(r"$s$")
  plot(t.s,t.wx,label=r'$w_x$')
  plot(t.s,t.wy,label=r'$w_y$')
  grid(True)
  legend()

