from numpy import array,zeros,ones,hstack

vars=2
order=3
class pol(object):
  def __init__(self,c=ones(1,dtype=float),e=ones((1,vars),dtype=int)):
    self.c=c
    self.e=e
  def __repr__(self):
    out=[]
    for i in range(len(self.e)):
      out.append(' '.join(map(str,[self.c[i] ]+self.e[i].tolist())))
    return '\n'.join(out)

  def __add__(a,b):
    c=hstack(a.c,b.c)
    e=hstack(a.e,b.e)
    ind=lexsort(e.T,axis=0)


  def merge




self=pol()
self
