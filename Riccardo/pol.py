class map:
  def __init__(self,**vars):
    self.__dict__.update(vars)

  def __call__(self,map=None,**vars):
    print map


from decimal import Decimal

class pol:
  """ Class for multivariate polinomials
from pol import pol as p
x=p('x')
y=p('y')
l=p('l')
driftx=x+l*y
drifty=y
beta=0.25
alpha=0
gamma=(1-alpha**2)/beta
twiss=beta*y**2+2*alpha*x*y+gamma*x**2
twiss(x=driftx,y=drifty)
drift={'x':x+l*y,'y':y}
twiss(**drift)
"""
#  _float=Decimal
  _float=float
  _zero=_float(0)
  _one=_float(1)
  def __init__(self,expr=''):
    self.coef={}
    self.vars=[]
    if expr:
      expr=expr.split(':')
      if len(expr)==1:
        self.vars=expr[0].split()
        self.coef[(1,)*len(self.vars)]=pol._one
      elif len(expr)==2:
        (vars,poly)=expr
        self.vars=vars.split()
        for mon in poly.split('+'):
          m=mon.split()
          exp=tuple(map(int,m[1:]))
          self.coef[exp]=self.coef.get(exp,pol._zero)+pol._float(m[0])

  def copy(self):
    new=pol()
    new.vars=self.vars[:]
    new.coef.update(self.coef)
    return new

  def addfloat(self,other):
    new=self.copy()
    i=(0,)*len(new.vars)
    new.coef[i]=new.coef.get(i,pol._zero)+pol._float(other)
    return new

  def mulfloat(self,other):
    new=self.copy()
    for i in new.coef:
      new.coef[i]*=pol._float(other)
    return new

  def addpol(self,other):
    new=pol()
    new.vars=list(set(other.vars+self.vars))
    for exp in self.coef:
      expd=dict(zip(self.vars,exp))
      newexp=tuple([expd.get(j,0) for j in new.vars])
      new.coef[newexp]=new.coef.get(newexp,pol._zero)+self.coef[exp]
    for exp in other.coef:
      expd=dict(zip(other.vars,exp))
      newexp=tuple([expd.get(j,0) for j in new.vars])
      new.coef[newexp]=new.coef.get(newexp,pol._zero)+other.coef[exp]
    return new

  def mulpol(self,other):
    new=pol()
    new.vars=list(set(other.vars+self.vars))
    for i in self.coef:
      for j in other.coef:
        c=self.coef[i]*other.coef[j]
        expi=dict(zip(self.vars,i))
        expj=dict(zip(other.vars,j))
        newexp=tuple([expi.get(k,pol._zero)+expj.get(k,pol._zero) for k in new.vars])
        new.coef[newexp]=new.coef.get(newexp,pol._zero)+c
    return new

  def __pow__(self,n):
    new=self.copy()
    if n==0:
      return 1
    elif n <0:
      return None
    else:
      for i in range(n-1):
        new*=self
    return new

  def __sub__(self,other):
    if type(other) in (pol._float,int):
      return self.addfloat(-1*other)
    else:
      return self.addpol(-1*other)

  def __add__(self,other):
    if type(other) in (pol._float,int):
      return self.addfloat(other)
    else:
      return self.addpol(other)

  def __mul__(self,other):
    if type(other) in (pol._float,int):
      return self.mulfloat(other)
    else:
      return self.mulpol(other)

  __radd__=addfloat
  __rsub__=__sub__
  __rmul__=mulfloat

  def __neg__(self):
    return self*-self._one

  def __pos__(self):
    return self

  def __repr__(self):
    s='%12s   ' %''
    for i in self.vars:
      s+=' %2s' % (i)
    s+=' :\n'
    for i in sorted(self.coef):
      s+="%12.5e   " % (self.coef[i])
      for j in i:
        s+=' %2d' %(j)
      s+=' +\n'
    return s[:-1]

  def __call__(self,map=None,**loc):
    for i in self.vars:
      loc.setdefault(i,pol(i))
    return eval(self.toexpr(),{},loc)

  def toexpr(self,fmt='%g'):
    m=[]
    for i in self.coef:
      c=[self.coef[i]]
      for j in range(len(i)):
        if i[j]==1:
          c+=['%s' %(self.vars[j])]
        elif i[j]!=0:
          c+=['%s**%d' %(self.vars[j],i[j])]
      if len(c)>1 and c[0]==self._one:
#        print c
        c.pop(0)
      else:
        c[0]=fmt % c[0]
      m+=['*'.join(c)]
    return ' + '.join(m)


if __name__== '__main__':
  x=pol('x')
  y=pol('y')
  a=x**2+y*x+1
  print a.toexpr()
  x=1+x+x**2
  print x
  x=x(x=x); print x
  x=x(x=x); print x
  x=x(x=x); print x
  print x(x=1.)
