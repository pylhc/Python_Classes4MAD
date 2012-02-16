from parser import expr
from symbol import sym_name
from token import tok_name
from py2texdata import texfun,texmod,texsym


d=tok_name.copy()
d.update(sym_name)

texop={
    'LPAR':r'\left(',
    'RPAR':r'\right)',
    'PLUS':r'+',
    'MINUS':r'-',
    'GREATER':r'>',
    'SMALLER':r'<',
    'TILDE':r'\~',
    'GREATEREQUAL':r'\ge',
    'SMALLEREQUAL':r'\le',
    'EQEQUAL':r'=',
    'NOTEQUAL':r'\ne',
    'in':r'\in',
    'is':r'\equiv',
    'COMMA':r',',
    }

class tex:
  def __init__(self,e):
    self.out=[]
    try:
      self.t=expr(e).tolist()[1]
      self.cc(self.t)
#    print self.t
      self.out=self.totex(self.t)
    except SyntaxError:
      self.out=e

  def __str__(self):
    return self.out


  def cc(self,t):
    for i,l in enumerate(t):
      if type(l) is int:
        t[i]=d[l]
      elif type(l) is list:
        self.cc(l)

  def totex(self,t):
    f=tex.__dict__
    if type(t) is list:
#      print "***parsing:",t[0],
      n='mk_'+t[0]
      if (n in f):
#        print "execute:",n
        return f[n](self,*t[1:])
      elif (t[0] in texop):
#        print "picking:",t[0]
        return texop[t[0]]
      else:
#        print "skip"
        return "".join([self.totex(i)  for i in t[1:] ])

  def mk_mymod(self,s):
    for (a,b) in texmod.iteritems():
      if s.endswith(a):
        s=b % s[:-len(a)]
    for (a,b) in texsym.iteritems():
      s=s.replace(a,b)
    return s
#    l=[(len(i[0]),i) for i in texsym.items()]
#    l.sort()
#    l.reverse()
#    l=[ i[1] for i in l]
#    j=0
#    ss=s
#    for (i,c) in enumerate(s):
#      for (a,b) in l:
#        if s[j:i]==a:
#          ss=ss.replace(a,b)
#          j=i
#    return ss


  def mk_NAME(self,*args):
    s=args[0]
    self.lastname=s
    s=[ self.mk_mymod(i) for i in  s.split('_')]
    try:
      return r'{%s}_{\roman %s}' %  (s[0],s[1])
    except IndexError:
      return '{%s}' % (s[0])

  def mk_NUMBER(self,*args):
    s=args[0]
    s=s.split('E')
    try:
      return r'{%s \cdot 10^{%s}' %  (s[0],s[1])
    except IndexError:
      return '{%s}' % (s[0])

  def mk_trailer(self,*args):
    p=args[0]
    t=self.totex(args[1])
    if p[1]=='(':
      p1=self.totex(args[0])
      p2=self.totex(args[2])
      return '%s%s%s' % (p1,t,p2)
    elif p[1]=='[':
      return r'_{%s}' % t.replace(',','')
    elif p[1]=='.':
      return r'_{\rm %s}' % t


#  def mk_atom(self,*args):
#    print len(args)
#    if len(args)==1:
#      return self.totex(args[0])
#    elif len(args)==3:
#      self.plevel+=1
#      t=self.totex(args[1])
#      p1=self.totex(args[0])
#      p2=self.totex(args[2])
#      return "".join([p1,t,p2 ])


  def mk_power(self,*args):
#    print args
#    print len(args)
    if len(args)==1:
      return self.totex(args[0])
    if len(args)==2:
      a1=self.totex(args[0])
      a2=self.totex(args[1])
      if args[1][1][0]=='LPAR':
        a1=self.totex(args[0])
        a2=self.totex(args[1])
        n=a1[1:-1]
        a=a2[6:-7]
        a=a.split(',')
        if n in texfun:
          s=texfun[n]
          na=s.count('%s')
          if len(a)<na:
            a.extend(['{}']*(na-len(a)))
          elif len(a)>na:
            a=a[:na]
          s=s % tuple(a)
          return '{%s}' % s
        else:
          return '{%s%s}' % (a1,a2)
      else:
        a1=self.totex(args[0])
        a2=self.totex(args[1])
        return '{%s%s}' % (a1,a2)
    elif len(args)==3:
      a1=self.totex(args[0])
      a2=self.totex(args[2])
      return '{{%s}^{%s}}' % (a1,a2)

  def mk_factor(self,*args):
    if len(args)==1:
      return self.totex(args[0])
    elif len(args)==2:
      a1=self.totex(args[0])
      a2=self.totex(args[1])
      return '{%s%s}' % (a1,a2)



  def mk_term(self,*args):
    out=[]
    if len(args)==1:
      return self.totex(args[0])
    elif len(args)>2:
      out.append(self.totex(args[0]))
      for i in range(len(args)/2):
        op=args[2*i+1]
        ar=args[2*i+2]
        ar=self.totex(ar)
        if op[1]=='*':
          if ar[1].isdigit():
             op=r'\cdot'
          else:
             op=r'\;'
          out.append(op)
          out.append(ar)
        elif op[1]=='/':
          ar1=out.pop()
          isover=False
          if ar1.endswith('right)'):
            isover=True
            ar1=ar1[6:-7]
          if ar.startswith(r'\left('):
            isover=True
            ar=ar[6:-7]
          if isover:
            out.append('{%s \over %s}'% (ar1,ar))
          else:
            out.append('{%s/%s}'% (ar1,ar))
        elif op[1]=='%':
          op=r'\%'
          out.append(op)
          out.append(ar)
        elif op[1]=='//':
          op=op[1]
          out.append(op)
          out.append(ar)
    return "".join(out)

  def mk_arith_expr(self,*args):
#    print args
#    print len(args)
    out=[]
    if len(args)==1:
      return self.totex(args[0])
    elif len(args)>2:
      out.append(self.totex(args[0]))
      for i in range(len(args)/2):
        op=args[2*i+1]
        ar=args[2*i+2]
        ar=self.totex(ar)
        op=self.totex(op)
        out.append(op)
        out.append(ar)
    return "".join(out)

#n=tex('sa_rew**2')
#print n.out
#n=tex('sa_rew*2/a')
#print n.out
#n=tex('sa_rew*2%a')
#print n.out
#n=tex('+sa_rew/2')
#print n.out
#n=tex('-2+2-a')
#print n.out
#n=tex('a/(b)')
#print n.out
#n=tex('a in (b)')
#print n.out
#n=tex('a + (b*c/4)')
#print n.out
#n=tex('a.x')
#print n.out
#n=tex('(x,a,c)')
#print n.out
#n=tex('a[x]')
#print n.out
#n=tex('sqrt(x,a,c)')
#print n.out
#n=tex('alpha_avec')
#print n.out
