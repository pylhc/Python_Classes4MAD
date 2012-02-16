try :
  from numpy import *
  isnumeric=False
except :
  from Numeric import *
  isnumeric=True



import re
#import sys
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

  def pattern(self,regexp) :
    """
    create an index list based on the selection
    t.select(t.pattern("IP"),t.name,t.s)
    """
    c=re.compile(regexp, re.IGNORECASE)
    indexes=[ self.idx[str(s)] for s in self.name if  c.search(str(s)) ]
    return indexes


  def cpattern(self,regexp) :
    """
    create a labels list  based on the selection
    """
    c=re.compile(regexp, re.IGNORECASE)
    lab = [l for l in self.labels if  c.search(l) ]
    return lab

  def extract(self,ind,lab) :
    """
    create an array based on row col regexp
    t.extract("MQXA","^[ab]")
    """
    if (type(ind)!=list) :
      ind=self.pattern(ind)
    if (type(lab)!=list) :
      lab=self.cpattern(lab)

    a=array([ self.data[l] for l in lab ])
    a=a[:,ind]
    return a.transpose()


  def elem(self,ra) :
    """
    return an index list of index for name=ra
    t.select(t.pattern("IP5"),t.name,t.s)
    """
    ra=ra.lower()
    ids=[ self.idx[str(s)] for s in self.name if  str(s).lower() == ra ]
    return ids

  def range(self,ra,rb,find=None) :
    """
    create index from element ra to element rb
    t.select(t.range("IP","IR11"),t.name,t.s)
    t.select(t.range("IP1","IR11$END",find=t.elem),t.name,t.s)
    """
    c=re.compile(ra, re.IGNORECASE)
    d=re.compile(rb, re.IGNORECASE)
    if not find: find=self.pattern
    try:
      ida=find(ra)[0]
      idb=find(rb)[0]
      return range(ida,idb+1)
    except IndexError:
      return []

  def srange(self,sa,sb) :
    """
    create index from s_a to s_b
    t.select(t.srange(0,100),t.name,t.s)
    """
    for ida in xrange(self.lines):
      if (self.s[ida] >= sa): break
    for idb in xrange(self.lines,0,-1):
      if (self.s[idb] <= sb): break
    return range(ida,idb+1)


  def select(self,index,*a) :
    """
    print a list based on the a list of index"
    t.select(t.pattern("MQ"),t.name,t.s)
    t.select(t.range("IP","IR11"),t.name,t.s)
    """
    if (len(a)==0):
      a=[getattr(self,i) for i in self.labels]

    a=list(a)
    if (type(a[-1])==str) :
      file=a.pop()
    else :
      file=None

    d=array([ i.take(index) for i in a ])
    d=d.transpose()
    slen=d.itemsize
    s="%-" + str(slen+1) + "s "
    s=s*d.shape[1]
#    print s
    labels=[]
    for i in a:
      try :
        labels.append(self.__class__.dict[id(i)])
      except KeyError:
        labels.append("USER")

    if (file) :
      out=open(file,'w')
    print s % tuple(labels)
    if (file) :
      out.write( (s+'\n') % tuple(labels))
#    print s % tuple([id(i) for i in a])
#    print self.__class__.dict
    for i in d:
      print s % tuple(i)
      if (file) :
        out.write( (s+'\n')  % tuple(i))
    return

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

  def pr(self,a,fn='out') :
    out=open(fn,'w')
    for i in a:
      i=map(str,i)
      l=len(i)
      fmt='%-20s ' * l
      line=fmt % tuple(i)
      print line
      out.write(line + '\n')
    out.close()


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

