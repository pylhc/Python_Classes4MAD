from math import sqrt


class interval:
  def __init__(self,d,u):
    self.u=u
    self.d=d

  def __str__(self):
    return "(%e,%e)" %(self.u,self.d)

  def __add__(x,y):
    return interval(x.d+y.d,x.u+y.u)

  def __sub__(x,y):
    return interval(x.d-y.u,x.u-y.d)

  def __neg__(x):
    return interval(-x.u,x.d)

  def __mul__(x,y):
    m1=x.u*y.u
    m2=x.d*y.u
    m3=x.u*y.d
    m4=x.d*y.d
    u=max(m1,m2,m3,m4)
    d=min(m1,m2,m3,m4)
    return interval(u,d)

class der:
  def __init__(self,x,dx):
    self.x=x
    self.dx=dx

  def __str__(self):
    return str(self.x)+','+str(self.dx)

  def __mul__(a,b):
    x=a.x*b.x
    dx=a.x*b.dx+a.dx*b.x
    return der(x,dx)

  def __add__(a,b):
    x=a.x+b.x
    dx=a.dx+b.dx
    return der(x,dx)



def test():
  a=interval(-1,1)
  b=interval(-1,1)
  print a*b
  c=der(a,b)
  print 

test()
