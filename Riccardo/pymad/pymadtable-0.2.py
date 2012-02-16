try :
  from numpy import *
  isnumeric=False
except :
  from Numeric import *
  isnumeric=True



import re
import sys
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
        f.pop(0) ; types=f
      else :   # data lines
        if (regexp):
          if (c.search(line)): save=True
          else : save=False
        else : save=True
        if (save):
          self.lines=self.lines+1
          f=map(self.conv,types,f)
  #        print self.labels
          for l in self.labels:
            d=f.pop(0)
            self.data[l].append(d)
            if (l == 'NAME'):
              self.idx[d]=self.lines
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
      col=self.name

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
      col=self.name

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
      col=self.name

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
      col=self.s

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
    t.select("MQ","BET")
    t.select(["MQT.13L5.B2", "MCBH.13L5.B2"],["betx", "bety"])
    t.select([1,3],[t.betx, t.bety])
    """

    if len(a)==0 :
      a1=""
      a2=self.name
    elif len(a)==1 :
      a1=a[0]
      a2=self.name
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
    """
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
    print s % tuple(labels)
    for i in ind:
      print s % tuple(a[:,i])
    return

  def write(self,*a) :
    """
    """
    file=a.pop()
    ind,col=self.select(*a)
    a=array( col )
    a=a[ind,:]
    return a

  def rsort(self,a):
    u={}
    for i in a: u[i]=True
    a=u.keys()
    a.sort()
    return a

  def radd(self,a,b) :
    u=a+b
    u=self.rsort(u)
    return u

  def rmul(self,a,b) :
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

if __name__ == '__main__':

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

  savefig('ring.eps',dpi=600)

  plot(s2.z,s2.x)


if __name__ == '__main__':
  import sys
  usage = 'Usage: %s madtable' % sys.argv[0]
#  try:
#       infilename = sys.argv.pop()
#  except:
#       print usage; sys.exit(1)
  x=madtable('twiss.lhcb1.data')
  y=madtable('twiss.lhcb1.data',r'IP.*')
  x.select(x.name,'IP')
  x.select([t.name,t.s],'mq')
#  print x.descs
#  print x.data

class elem:
  pass

class madinput:
  """
  """
  def __init__(self, filename)::
    self.v={}
    self.e={}
    filename="dipfirst.thin.seq"
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
    relem=re.compile("(^[\w\.]+):([\w\.]+),(.+)?;")


    for i in r.finditer(s):
      st=i.group()
      # variable parsing
      if (rvar.match(st)):
        tok=rvar.match(st).groups()
        self.v[ tok[0] ]=tok[1]
      # elements parsing
      elif (relem.match(st)):
        tok=relem.match(st).groups()
        el=elem()
        el.name=tok[0]
        el.parent=tok[1]
        el.sequence=[]
        # attribute parsing
        attr=list(tok[2])
        ch=False
        for i in range(len(attr)):
          if attr[i]=='{':
            ch=True
          if attr[i]=='}':
            ch=False
          if (ch) and (attr[i]==',') :
            attr[i]=' '
        attr=(''.join(attr)).split(',')
        for a in attr:
          stok=rvar.match(a).groups()
          setattr(el,stok[0],stok[1])
        # store elem in the dictionary
        self.e[ tok[0] ]=el
        # sequence management
        if (el.parent=='sequence'):
          currseq=el.name
          el.elements=[]
        elif (currseq) :
          el.sequence.append(currseq)
          self.e[currseq].elements.append(el.name)
      # sequence parsing
      elif (st=='endsequence;'):
        currseq=None

#   globals().update(self.e)




    
