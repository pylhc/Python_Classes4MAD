try :
  from numpy import *
  isnumeric=False
except :
  from Numeric import *
  isnumeric=True



import re
import sys
from math import acos, asin, atan, atan2, ceil, cos, cosh, degrees, e, exp, fabs, floor, fmod, frexp, hypot, ldexp, log, log10, modf, pi, pow, radians, sin, sinh, sqrt, tan, tanh
#from string import split

class madtable:
  """
  import a mad table, examples:
  from pymadtable import *
  t=madtable('twiss.lhcb1.x1y5.data')
  t.select(t.pattern("IP"),t.name,t.s)
  t.select(t.range("IP"),t.name,t.s)
  t.select(t.rmul(t.range('IP5','MQY.*R5'), t.radd(t.pattern('^MQ'),[i-1 for i in t.pattern('^MQ')]) ),t.name,t.s-t.s[t.elem('IP5')])
  """
  dict={}
  def __init__(self, filename, regexp=None):
    self.data={}
    self.descs={}
    self.idx={}
    self.lines=-1
    if (regexp): c=re.compile(regexp, re.IGNORECASE)
    for line in open(filename):
      f=line.split()
      if (f[0] == '@'):  # descriptor lines
        self.descs[f[1]]=self.conv(f[2],f[3])
      elif ( f[0] == '*'): # self.labels lines
        f.pop(0) ; self.labels=f
        for l in self.labels: self.data[l]=[]
      elif (f[0] == '$'):  # type lines
        f.pop(0) ; self.types=f
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
      if (not (isnumeric and (l == 'NAME') )) :
        self.data[l]=array(self.data[l])
      setattr(self,l.lower(),self.data[l])
      setattr(self,l,self.data[l])
      self.__class__.dict[id(self.data[l])]=l


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
    idx=[]
    for i in xrange(0,len(col)) :
      c.search(str(col[i])) and idx.append(i)
    return idx

  def name2idx(self,name) :
    """
    t.name[t.name2idx(["MQXC.3L5.B2","MQXB.B2L5.B2"])]
    """
    idx=[ self.idx[str(s)] for s in name]
    return idx

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

    idx=[]
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
      return range(ida,idb+1)
    except IndexError:
      return []

  def nrange(self,sa,sb,*a) :
    """
    t.name[t.nrange(1,7)]
    t.name[t.nrange(100,200,t.betx)]
    """
    if len(a)==1 :
      col=a[0]
    else :
      col=self.data[self.labels[1]]

    idx=[]
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
    return col

  def lpattern(self,regexp) :
    """
    t.lpattern("BET")
    """
    c=re.compile(regexp, re.IGNORECASE)
    lab = [l for l in self.labels if  c.search(l) ]
    return lab

  def lab2col(self,lab) :
    """
    t.lab2col(t.lpattern("BET"))
    """
    col=[getattr(self,l) for l in lab]
    return col

  def select(self,*a) :
    """
    t.select("MQ","N BE ")
    t.select(["MQT.13L5.B2", "MCBH.13L5.B2"],["betx", "bety"])
    t.select([1,3],[t.betx, t.bety])
    """

    if len(a)==0 :
      a1=".*"
      a2=[ self.data[self.labels[0]] ]
    elif len(a)==1 :
      a1=a[0]
      a2=[ self.data[self.labels[0]] ]
    elif len(a)>1 :
      a1=a[0]
      a2=a[1]

    idx=[]
    try :
      if   (type(a1) is type("")) :
        for a in a1.split() :
          idx.extend(self.pattern(a))
      elif type(a1) is type([]):
        if a1 :
          if   type(a1[0]) is type("") :
            idx=self.name2idx(a1)
          elif type(a1[0]) is type(0) :
            idx=a1
    except :
      print "Unexpected error:", sys.exc_info()[0]
      idx=[]

    col=[]
    try :
      if   (type(a2) is type("")) :
        for a in a2.split() :
          col.extend(self.cpattern(a))
      elif type(a2) is type([]):
        if a2 :
          if   ( type(a2[0]) is type("") ) :
            col=self.lab2col(a2)
          elif ( type(a2[0]) is type(array([])) ) :
            col=a2
    except :
      print "Unexpected error:", sys.exc_info()[0]

    return idx,col


  def extract(self,*a) :
    """
    t.extract("MQ","BET")
    t.extract("MQ",["NAME","S"])
    t.extract("MQ",[t.name,t.s])
    """
    ind,col=self.select(*a)
    a=array(col).transpose()
    a=a[ind,:]
    return a

  def show(self,*a) :
    """
    t.show("MQ","BET")
    t.show("MQ","BET","out.dat")
    """
    if (len(a)==3) :
      a=list(a)
      file=a.pop()
      f=open(file,"w")
      a=tuple(a)
    else:
      f=None
    ind,col=self.select(*a)
    labels=[]
    for i in col:
      try :
         labels.append(self.__class__.dict[id(i)])
      except KeyError:
        labels.append("USER")

    s="%-"+str(15)+"s "
    s=s*len(col)
    a=array(col)
    if f :
      s+="\n"
      f.write(s % tuple(labels))
      for i in ind:
        f.write(s % tuple(a[:,i]))
      f.close()
    else :
      print s % tuple(labels)
      for i in ind:
        print s % tuple(a[:,i])
    return

  def rsort(self,a):
    """
    """
    u={}
    for i in a: u[i]=True
    a=u.keys()
    a.sort()
    return a

  def radd(self,a,b) :
    """
    i.show( i.radd(i.pattern("E"),i.pattern("I")) )
    """
    u=a+b
    u=self.rsort(u)
    return u

  def rmul(self,a,b) :
    """
    i.show( i.radd(i.pattern("E"),i.pattern("I")) )
    """
    u=[]
    for i in a :
      try :
        b.index(i)
        u.append(i)
      except ValueError:
        pass
    return u

class envelope:
  """
  find aperture given a survey and a twiss table, example:
  s1=madtable('survey.lhcb1.data')
  s2=madtable('survey.lhcb2.data')
  t1=madtable('twiss.lhcb1.data')
  t2=madtable('twiss.lhcb2.data')
  ap1=aperture(s1,t1)
  ap2=aperture(s1,t1)

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
      wm=matrixmultiply(thetam,matrixmultiply(phim,psim))
      ex=matrixmultiply(wm,array([1,0,0]))
      ey=matrixmultiply(wm,array([0,1,0]))
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



class var:
  def __init__(self,e):
    self.expr=e

  def __float__(self):
    return eval(self.expr,self._env)

  def __neg__(self):
    return -1*eval(self.expr[1:],self._env)

  def __coerce__(self,other):
    return eval(self.expr,self._env),float(other)

  def __repr__(self):
    rep= 'expr=%s, value=%e' % (self.expr,self.value())
    return  rep

  def value(self):
    try:
      value=float(eval(self.expr,self._env))
    except :
      value=0
    return  value


class elem:
  """
  mad elements
  """
  def __init__(self, name,**kwds):
    self.name=name
    # store attribute via dictionary
    self.__dict__.update(kwds)

  def show(self):
    globals().update(self._env)
    for k in self.__dict__:
      if not k.startswith('_'):
        a=getattr(self, k)
  #        if a is type(""):
        try :
          if a is type([]):
            value=map(eval,a)
          else :
            value=eval(a)
          print '%s=%s value: %s' % (k,a,value)
        except:
          print '%s=%s' % (k,a)


  def __str__(self):
    return  self.name

  def __repr__(self):
    return  self.name

class madinput:
  """
  import pymad as m
  reload(m); i=m.madinput("dipfirst.thin.seq")
  """
  def __init__(self, filename):
    self.v={}
    self.e={}
    f=open(filename,"rU")
    s=f.read()
    s=s.lower()
    r=re.compile("!.*")
    s=r.sub("",s)
    s=s.replace('\n','')
    s=s.replace(' ','')
    r=re.compile("[^;]*;")
    #rvar=re.compile("([\w\.]+):?=([\w\.]+);")
    rvar=re.compile("(^[\w\.]+):?=([^;:,]+)")
    relem=re.compile("(^[\w\.]+):([\w\.]+)(,.+)?;")
    relem2=re.compile("(^[\w\.]+)(,.+);")
    rdot1=re.compile("([a-z])\.([a-z0-9])")
    rdot2=re.compile("([a-z0-9])\.([a-z])")
    rdot2=re.compile("([a-z0-9])\.([a-z])")


    currseq=None
    for i in r.finditer(s):
      st=i.group()
      st=st.replace("const","")
      st=rdot1.sub(r"\1_\2",st)
      st=rdot2.sub(r"\1_\2",st)
      st=st.replace("..","__")
      # variable parsing
      if (rvar.match(st)):
        tok=rvar.match(st).groups()
        # store var in the dictionary
#        self.v[ tok[0] ]=tok[1]
        v=var(tok[1])
        v._env=self.v
        self.v[ tok[0] ] = v
      # elements parsing
      elif (relem.match(st) or relem2.match(st)):
#        print st
#        print currseq
        try :
          tok=relem.match(st).groups()
          el=elem(tok[0])
          try :
            el.parent=self.e[tok[1]]
          except KeyError:
            el.parent=tok[1]
          attr=tok[2]
        except AttributeError:
          tok=relem2.match(st).groups()
          try:
            el=self.e[tok[0]]
            attr=tok[1]
          except KeyError:
            print "No parent:" + tok[0]
            el=elem(tok[0])
        # store elem in the dictionary
        el._env=self.v
        self.e[ el.name ]=el
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
            try :
              stok=rvar.match(a).groups()
              av=stok[1]
              if av[0]=='{':
                av=av[1:-1].split()
              setattr(el,stok[0],av)
            except AttributeError :
              pass
#        el.show()
        # sequence management
        if (el.parent=='sequence'):
#          raise
#        if (st.startswith('bpmsw')): raise
          currseq=el.name
          el.elements=[]
        elif (currseq) :
          try:
            el.sequence.append(currseq)
          except AttributeError:
            el.sequence=[ currseq ]
          self.e[currseq].elements.append(el.name)
#        el.show()
#        print
      # sequence parsing
      elif (st=='endsequence;'):
        currseq=None

    self.__dict__.update(self.e)
    self.__dict__.update(self.v)
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
    ns=self.v.keys()
    ns.sort()
    for n in ns:
      v=self.v[n]
      f.write( "%-20s %e %s\n" % (n,v.value(),v.expr) )
    f.close()

