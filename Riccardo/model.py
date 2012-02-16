from os import system

class model(object):

  def __init__(self,s=None):
    self._formula={}
    self._compiled={}
    self._fdep={}
    self._bdep={}
    self._computed={}
    self._value={}
    self._comments={}
    self._units={}
    self._tools={}

    self._chain=[]
    self._dbg=0
    self._visited=[]
    if s:
      self.update(s)

  def setdebug(self,i):
    self._dbg=i

  def __delitem__(self, key):
    self._debug('del','%s' % (key),0)
    self._bdep[key]={}
    if key not in self._fdep:
      self._fdep[key]={}

    for i in self._bdep[key]:
      self._bdep[key].pop(i,0)
      self._fdep[i].pop(key,0)

    del self._formula[key]
    del self._compiled[key]
    if key in self._computed:
      del self._computed[key]
    del self._value[key]

  __delattr__=__delitem__

  def __setitem__(self, key, formula):
    self._debug('set','%s -> %s' % (key,formula),0)
    self._bdep[key]={}
    if key not in self._fdep:
      self._fdep[key]={}

    for i in self._bdep[key]:
      self._bdep[key].pop(i,0)
      self._fdep[i].pop(key,0)

    self._formula[key]=formula
    self._compiled[key]=compile(formula,'key','eval')
    if key in self._computed:
       self._computed.pop(key,0)
    self._value[key]=self[key]
    self._last=None

    self._propfdep(key)

  def __getitem__(self,key):
    if key in self._tools:
      raise KeyError

    self._chain.append(key)
    self._debug('chain',self._chain,1)
    if len(self._chain)==2:
      self._fdep[key][self._chain[0]]=1
      self._bdep[self._chain[0]][key]=1

    if key in self._computed:
      value=self._value[key]
      self._debug('recall','%s -> %s' % (key,value),0)
    else:
      self._debug('eval','%s' % (key),0)
      value=eval(self._compiled[key],self._tools,self)
      self._debug('value','%s -> %s' % (key,value),0)
      self._computed[key]=1
      self._value[key]=value

    self._chain.pop()

    return value

  def _propfdep(self,key):
    for i in self._fdep[key]:
      self._flag(i)
      self._propfdep(i)

  def _flag(self,key):
    self._debug('flag','%s' % (key),0)
    if key in self._computed:
      self._computed.pop(key,0)

  def __setattr__(self, key,formula):
    if key.startswith('_'):
      self.__dict__[key]=formula
    else:
      if type(formula) is not str:
        formula=repr(formula)
      f=formula.split('#')
      self[key]=f[0].strip()
      self._comments[key]=""
      self._units[key]=""
      if len(f)>2:
        self._comments[key]=f[2].strip()
      if len(f)>1:
        self._units[key]=f[1].strip()

  def __getattr__(self, key):
    self._debug('getattr','%s' % (key),1)
    if key in self._formula:
      return self[key]
    elif key=='__members__':
      return self._formula.keys()
    else:
      raise AttributeError

  def _debug(self,msg1,msg2,level):
    if level<self._dbg:
      print "%-7s: %s" % (msg1,msg2)

  def _visit(self,key):
    if key not in self._visited:
      for k in self._bdep[key]:
        self._visit(k)
      self._visited.append(key)

  def print_dep(self):
    for key in self._fdep:
      for i in self._fdep[key]:
        print i,'depends on',key


  def update(self,other):
    other._visited=[]
    for k in other._fdep:
      other._visit(k)
    for i in "_comments _units _tools".split():
      getattr(self,i).update(getattr(other,i))
    for i in other._visited:
      self[i]=other._formula[i]


  def graph(self,fn="out.ps",formula=1,value=1,comment=1):
    out=["""digraph W {
   node [shape=Mrecord];
    rotate=90;"""]
    for i in self._formula:
      out+=[' "%s" [label="' % i]
      if comment and (i in self._comments):
        out+=['%s\\n' % self._comments[i]]
      out+=['%s' %i]
      if formula:
        out+=['=%s' % self._formula[i]]
      if value:
        out+=['\\n=%s %s' % (getattr(self,i),self._units[i])]
      out+=['"];\n']

    for key in self._bdep:
      if self._bdep[key]:
        for i in self._bdep[key]:
          out+=['"%s" -> "%s";\n' % (i,key)]
    out+=['}']
    out=''.join(out)
    dfn=fn.split('.')[0]+".dot"
    open(dfn,"w").write(out)
    system("dot -Tps %s >%s"%(dfn,fn))
    print "Print %s" % fn
    return

  def getformula(self,i):
    return self._formula[i]

  def getcomment(self,i):
    return self._comments[i]

  def getunit(self,i):
    return self._units[i]

  def gettex(self,i):
#    from eng_notation import num_to_str
    d={}
    d['beta']=r'{\beta}'
    d['sigma']=r'{\sigma}'
    d['Delta']=r'{\Delta}'
    d['emit']=r'{\varepsilon}'
    d['phi']=r'{\phi}'
    d['star']=r'{^*}'
    n=i.split('_')
    n0=n[0]
    for i in d:
      n0=n0.replace(i,d[i])
    out=[n0]
    if len(n)==2:
      out.append('_{\\rm{%s}}' % n[1])
#    (v,s)=num_to_str(self[i],d=4,tex=1)
#    out.append(' &= %g \\rm{%s %s}\\\\' % (float(v),s,self.getunit(i)))
    return "".join(out)





if __name__== '__main__':
  m=model()
  m.a='4'
  m.b='4*a'
  del m.b
#  m.d='5'
#  m.c='4*b+d'
#  m.print_dep()
#  m.a='5'
#  m.b='44'
#  m.print_dep()

  #print 'prova'
#  m.x='0'
#  m.p='x**2'
#  m.p2='p**2'

#  print [m.p2 for m.x in [1,2,4]]
#  m.x='numpy.array([1,2,3])'
#  m.p2
  #print [m['p'] for m.x in [1]]

