#!/usr/bin/python

import re
import sys
import elementtree.ElementTree as ET
from cStringIO import StringIO
import textopng
import py2tex

class wiki:
  lead={'=': 'heading',
        '*': 'ulist',
        '{': 'meta',
        '#': 'olist',
        ';': 'dlist',
        '$': 'equation',
        '[': 'figure',
        '!': 'script',
        '|': 'table',
        ' ': 'code',
         }

  def __init__(self,f,prefix="",pretex="",basedir=""):
    self.prefix=prefix
    self.pretex=pretex
    self.basedir=basedir
    textopng.base=basedir
    self.root=ET.Element("html")
    self.head=ET.SubElement(self.root,"head")
    self.meta=ET.SubElement(self.root,"meta")
    self.body=ET.SubElement(self.root, "body")
    self.div=ET.SubElement(self.body, "div")
    self.cont=self.div
    self.list=[self.div]
    self.scripteof=''
    self.scriptsource=''
    self.scriptprog=''
    self.parse(f)

  def parsescript(self,line):
    a=re.match('!([A-z \-\.]*)(<<)?([A-z]*)',line)
    if a:
      (self.scriptprog,s,self.scripteof)=a.groups()


  def scriptexe(self):
    import os
    import md5
    self.scriptsource=self.scriptsource[0:-1]
    h=md5.new(self.scriptsource).hexdigest()
#    print self.scriptprog+' '+h
#    print self.scripteof
    if not os.path.exists(h+'.out'):
      open(h,'w').write(self.scriptsource)
      out=os.popen(self.scriptprog+' '+h)
      open(h+'.out','w').write(out.read())
    out=open(h+'.out')
    ET.SubElement(self.cont, "pre").text=out.read()
    return

  def  parsetable(self,line):
    if self.cont is self.div:
      self.cont=ET.SubElement(self.cont, "table")
      self.cont.attrib['border']='1'
      self.cont.attrib['cellspacing']='0'
    row=ET.SubElement(self.cont, "tr")
    cols=line.split('|')[1:-1]
    for i in cols:
      if i:
        if i[0]==',':
          if c.attrib.has_key('colspan'):
            c.attrib['colspan']=str(int(c.attrib['colspan'])+1)
          else:
            c.attrib['colspan']='2'
          continue
        else:
          c=ET.SubElement(row, "td")
        if i[-1]==',':
          c.attrib['rowspan']='2'
          i=i[0:-1]
        if i==' '*len(i):
          c.text=u'\xA0'
        else:
          c.text=i
        if i[0]!=' ' and i[-1]==' ':
          c.attrib['align']='left'
        elif i[0]==' ' and i[-1]!=' ':
          c.attrib['align']='right'
        elif i[0]!=' ' and i[-1]!=' ':
          c.attrib['align']='justify'
        elif i==' '*len(i):
          c.text=u'\xA0'
        else :
          c.attrib['align']='center'

  def  parseequation(self,line):
    text=str(py2tex.tex(line[1:-1]))
    fn=textopng.textopng('$$%s$$' % text)
    ET.SubElement(self.cont, "p").text=line
    ET.SubElement(self.cont, "p").text=text
    l=ET.SubElement(self.cont, "center")
    l=ET.SubElement(l,"img")
    l.attrib['alt']=line
    l.attrib['src']=fn

  def parsecode(self,line):
    if self.cont.tag!='pre':
      self.cont=ET.SubElement(self.cont, "pre")
      self.cont.text=line[1:]
    else:
      self.cont.text+='\n'+line[1:]

  def parsefigure(self,line):
    if line[1:-1].split('.')[-1] in ('png','jpg','gif'):
      l=ET.SubElement(self.cont, "img")
      l.attrib['alt']=line[1:-1]
      l.attrib['src']=line[1:-1]
    else:
      l=ET.SubElement(self.cont, "a")
      l.text=line[1:-1]
      l.attrib['href']=line[1:-1]

  def parseulist(self,line):
    level=0
    for c in line:
      if c=='*':
        level+=1
      else:
        break
    while (level>len(self.list)-1):
      self.list.append(ET.SubElement(self.list[-1], "ul"))
    while (level<len(self.list)-1):
      self.list.pop()
    self.cont=ET.SubElement(self.list[-1],"li")
    self.cont.text=line[level:]

  def parseolist(self,line):
    level=0
    for c in line:
      if c=='#':
        level+=1
      else:
        break
    while (level>len(self.list)-1):
      self.list.append(ET.SubElement(self.list[-1], "ol"))
    while (level<len(self.list)-1):
      self.list.pop()
    self.cont=ET.SubElement(self.list[-1],"li")
    self.cont.text=line[level:]

  def parsedlist(self,line):
    level=0
    for c in line:
      if c==';':
        level+=1
      else:
        break
    while (level>len(self.list)-1):
      self.list.append(ET.SubElement(self.list[-1], "dl"))
    while (level<len(self.list)-1):
      self.list.pop()
    d=line[level:].split(':',1)
    ET.SubElement(self.list[-1],"dt").text=d[0]
    self.cont=ET.SubElement(self.list[-1],"dd")
    self.cont.text=d[1]

  def parseheading(self,line):
    level=None
    for i in range(4,0,-1):
      if (line[0:i]=='='*i) and (line[-i:]=='='*i):
        level=i
        line=line[level:-level]
        break
    if level:
      ET.SubElement(self.div, "h%i" %level).text=line
    else:
      self.parseparagraph(line)

  def parsemeta(self,line):
    l=line[1:-1].split(':',1)
    self.meta.attrib[l[0]]=l[1].strip()

  def parseparagraph(self,line):
    if self.cont is self.div:
      self.cont=ET.SubElement(self.div, "p")
      self.cont.text=line
    else:
      self.cont.text+=' '+line

  def parseblank(self,line):
    self.cont=self.div
    self.list=[self.div]



  def parse(self,f):
    for line in f.split('\n'):
      if self.scriptprog=='':
        if line:
          if line[0] in self.lead.keys():
            t=self.lead[line[0]]
          else:
            t='paragraph'
        else:
          t='blank'
#        print "%-10s >%s" % (t,line)
        if hasattr(self,'parse'+t):
          getattr(self,'parse'+t)(line)
      else:
        if line==self.scripteof:
          self.scriptexe()
          self.scripteof=''
          self.scriptsource=''
          self.scriptprog=''
        else:
          self.scriptsource+=line+'\n'

  def __str__(self):
    s=StringIO()
    ET.ElementTree(self.root).write(s)
    s.reset()
    s=s.read()
    s=re.sub("''([^']*)''",r"<em>\1</em>",s)
    s=re.sub("\+([^+]*)\+",r"<strong>\1</strong>",s)
    s=re.sub("-([^\-]*)-",r"<del>\1</del>",s)
    s=re.sub("_([^_]*)_",r"<ins>\1</ins>",s)
    s=re.sub(",,([^,]*),,",r"<code>\1</code>",s)
    s=re.sub("\$([^$]*)\$",self.maketexformula,s)
    s=re.sub("\[\[\[([^\] ]*)\]\]\]",r'<a name="\1"/>',s)
    s=re.sub("\[\[\[([^\] ]*) ([^\]]*)\]\]\]",r'<a name="\1">\2</a>',s)
    s=re.sub("\[([^\] ]*)\]",r'<a href="\1">\1</a>',s)
    s=re.sub("\[([^\] ]*) ([^\]]*)\]",r'<a href="\1">\2</a>',s)
    s=re.sub("\[\[([^\] ]*)\]\]",r'<a href="'+self.prefix+r'"\1">\1</a>',s)
    s=re.sub("\[\[([^\] ]*) ([^\]]*)\]\]",r'<a href='+self.prefix+r'"\1">\2</a>',s)
    pat= re.compile(r'((ftp|http)://[\w-]+(?:\.[\w-]+)*(?:/[~\&\w\-\.?=]*)*)')
    s= pat.sub('<a href="\\1">\\1</a>', s)
    pat = re.compile(r'([\w-]+(?:\.[\w-]+)*@[\w-]+(?:\.[\w-]+)*)')
    s = pat.sub('<a href="/cgi-bin/mail.py?p=\\1">\\1</a>', s)
    pat = re.compile(r'(00[0-9]{4,18})')
    s = pat.sub('<a href="/cgi-bin/sms.py?p=\\1">\\1</a>', s)
    return s

if __name__=='__main__':
  w=wiki(sys.stdin.read())
  print w
  open("o.html",'w').write(str(w))

