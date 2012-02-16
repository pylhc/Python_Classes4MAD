import os
import md5
import shutil

base=os.path.expanduser('~')
base=os.path.join(base,'.texcache')

pre=r"""
\input amstex
\input amssym
%\loadmsam
%\loadmsbm
\loadbold
\UseAMSsymbols
\nopagenumbers
"""
post="""
\end
"""

def textopng(text):
  if not os.path.exists(base):
    os.mkdir(base)
  h=md5.new(text).hexdigest()
  fn=os.path.join(base,h)
  if not os.path.exists(fn+'.png'):
    open(fn+'.tex','w').write("\n".join((pre,text,post)))
#    os.system('(cd %s; ls)' % (base,))
    os.system('(cd %s;tex >/dev/null %s)' % (base,h))
    os.system('(cd %s;dvipng >/dev/null -bg Transparent -T tight -o %s.png %s )' % (base,h,h))
#    if not base:
#      shutil.move(h+'.png',fn+'.png')
#      os.symlink(fn+'.png',h+'.png')
#  else:
#    if not base:
#      os.symlink(fn+'.png',h+'.png')
  return fn+'.png'
