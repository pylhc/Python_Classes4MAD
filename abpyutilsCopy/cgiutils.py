import cgi
import cgitb; cgitb.enable()
import time
import Cookie
import os


from urllib import urlencode



class tag:
  """Usage:
  t1=tag('x','y',z='b')
  print t1 => <x z="b">y</x>
  print tag('a','b',t1,'c',d='e') => <a d="e">b<x z="b">y</x>c</a>
  t1.add('f')
  print t1 => <x z="b">y<f/></x>
  t1.cat('f')
  print t1 => <x z="b">y<f/>f</x>
    """
  def __init__(self,*text,**attr):
    """Example:
    tag('a','my link',href='index.html')
    """
    text=list(text)
    self.name=str(text.pop(0))
    self.text=text
    self.attr=attr

  def __str__(self):
    if self.name.startswith('<'):
      return self.name
    s=["<%s" % self.name]
    for k,v in self.attr.items():
      s+=[' %s="%s"' % (k,htmlfixtext(str(v)))]
    if self.text:
      s+=[">"]
      for i in self.text:
        if hasattr(i,'text'):
          s+=[str(i)]
        else:
          s+=[htmlfixtext(str(i))]
      s+=["</%s>" % (self.name)]
    else:
      s+=["/>"]
    return "".join(s)

  def add(self,*text,**attr):
    """Define and add a new tag to self"""
    text=list(text)
    name=text.pop(0)
    t=tag(str(name),*text,**attr)
    self.text.append(t)
    return t

  def cat(self,*text):
    """Add new texts or tags to self"""
    self.text.extend(text)
    return self


def cmd(text,base,query={},**args):
  if query:
    base+='?'+urlencode(query)
  return tag('a',text,href=base,**args)


def htmlfixtext(text):
  l=list(text)
  d={'"' :'&quot;', "'" :'&#39;', "&" :'&amp;',
     '<' :'&lt;', '>' :'&gt;', }
  l=[d.get(i,i) for i in l]
  return ''.join(l)


def simplepage():
  page=tag('html')
  head=page.add('head')
  body=page.add('body')
  return page,head,body

def parsedate(s,fmt='%Y-%m-%d %H:%M:%S'):
  cmd="date -d '%s' " % s
  cmd+="+'%s'" % fmt
  date=os.popen(cmd).readline()
  return date.strip()



class cgipage:
  def __init__(self):
    self.page,self.head,self.body=simplepage()
    self.environ=os.environ
    self.formob=cgi.FieldStorage(keep_blank_values=True)
    if hasattr(self,'pre'):
      getattr(self,'pre')()
    for fn in self._formlist('cmd',['default']):
      if not fn.startswith('_'):
        if hasattr(self,fn):
          getattr(self,fn)()
    if hasattr(self,'post'):
      getattr(self,'post')()

  def _form(self,k,default=None):
    form=self.formob
    if k in form:
      v=form[k]
      if type(v) is not list:
        return v.value
      formk=[]
      for i in v:
        if i.filename:
          formk.append([i.filename,i.file])
        else:
          formk.append(i.value)
      return formk
    return default

  def _formkeys(self):
    return self.formob.keys()

  def _formlist(self,k,default=None):
    v=self._form(k,default)
    if type(v) is not list:
      return [v]
    else:
      return v


  def dump(self):
    cookie=Cookie.SimpleCookie()
    cookie['lastvisit'] = time.asctime()
    print cookie

    header='''Content-type: text/html\n\n'''
    print header

    page,head,body=simplepage()

    body.add('h2','Form')
    table=body.add('table')
    for k in self._formkeys():
      for v in self._formlist(k):
        table.add('tr',tag('td',k),tag('td',v))

    body.add('h2','os.environ')
    table=body.add('table')
    for k,v in self.environ.items():
      table.add('tr',tag('td',k),tag('td',v))

    body.add('h2','Cookies')
    if 'HTTP_COOKIE' in os.environ:
      cookie.load(os.environ['HTTP_COOKIE'])
      table=body.add('table')
      for k,v in cookie.items():
        table.add('tr',tag('td',k),tag('td',v.value))
    else:
      body.add('p','No cookies')

    print page

  def default(self):
    self.dump()




def mktable(cols,data):
  t=tag('table')
  for i in cols:
    r=t.add('tr')
    for j in i:
      r.add('th',j)
  for i in data:
    r=t.add('tr')
    for j in i:
      r.add('td',j)
  return t


def mkhidden(name,value=""):
  return tag('input',type='hidden',name=name,value=value)

def mkinput(name,value=""):
  return tag('input',type='text',name=name,value=value)

def mkselect(name,value=[],select=None):
  s=tag('select',name=name)
  for k in value:
    if k==select:
      s.add('option',k,value=k,selected="selected")
    else:
      s.add('option',k,value=k)
  return s

def mksubmit():
  return tag('input',type='submit',value='Submit')


def printform(act,cols,data):
  f=tag('form',action=act,method='post')
  t=f.add('table')
  r=t.add('tr')
  for j in cols:
    r.add('th',j)
  r=t.add('tr')
  for j in cols:
    d=data.get(j,'')
    if hasattr(d,'pop'):
      s=r.add('td').add('select',name=j)
      for k in d:
        s.add('option',k,value=k)
    else:
      r.add('td',tag('input',type='text',name=j,value=d))
  f.add('input',type='submit',value='Submit')
  return f


