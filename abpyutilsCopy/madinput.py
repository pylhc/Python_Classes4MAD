#import re
#import sys
#import os
#import pexpect
#from math import acos, asin, atan, atan2, ceil, cos, cosh, degrees, e, exp, fabs, floor, fmod, frexp, hypot, ldexp, log, log10, modf, pi, pow, radians, sin, sinh, sqrt, tan, tanh
#from pylab import *
#from numpy import *
#import readline
#from user import home
#import gzip
#from string import split

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

