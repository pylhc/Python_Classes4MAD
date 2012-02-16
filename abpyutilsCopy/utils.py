from numpy import log10

def gt(a,b):
  if a<b:
    return (b-a)**2
  else:
    return 0.

def tr(a,b,c):
  if a<b:
    return (b-a)**2
  elif a>c:
    return (a-c)**2
  else:
    return 0.

def eq(a,b):
  return (a-b)**2

def rng(a,b,c):
  return (a>b) & (a<c)

class container(object):
  def __init__(self,**args):
    self.__dict__.update(args)
  def __repr__(self):
    out=['%s(' % self.__class__.__name__]
    out+=['%s=%s,' % (k,v) for k,v in  self.__dict__.items() if not k.startswith('__')]
    out+=[')']
    return ''.join(out)

def numtostr(n):
  if abs(n)>0:
    l=log10(abs(n))
    if l<-3:
      l=int(l)
      o=(l-1)//3*3
      fmt='%%%d.3fe%+03d' % (8,o)
      n=n/10**o
    elif -3<=l and l< 0: fmt='%11.6f '
    elif  0<=l and l< 3: fmt='%8.3f    '
    elif  3<=l:
      l=int(l)
      o=(l)//3*3
      fmt='%%%d.3fe%+03d' % (8,o)
      n=n/10**o
  else:
    fmt='%4.0f.       '
  return fmt % n

#print numtostr(0)
#for i in range(-10,10):
#  print '"%s"' % numtostr(-1.23456789*10**i)
#print numtostr(0.009000),'a'
#print numtostr(0.09000),'a'


class opticsutils(object):
  def maxbetx(f):
    return f.betx+f.alfx**2/f.betx/abs(f.kn1l/f.l)
  def maxbety(f):
    return f.bety+f.alfy**2/f.bety/abs(f.kn1l/f.l)
  def chromx(f):
    return -sum(f.kn1l*f.betx)/4/pi
  def chromy(f):
    return sum(f.kn1l*f.bety)/4/pi

def minbetapeak(bs,ls,k):
  return (2+ 1.11*exp(-ls*sqrt(k))  + ls*sqrt(k))**2/(k*bs)



