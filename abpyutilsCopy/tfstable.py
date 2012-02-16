from tfsdata import tfsdata
from table import table
from plotfunc import plotfunc

def _pythonname(string):
  string=string.replace('.','_')
  string=string.replace('$','_')
  return string.lower()

class tfstable(table,plotfunc):
  _is_s_begin=False
  def __init__(self,filename,rowattr=True):
    data=tfsdata().fromfile(filename)
    col_names=data.labels
    if 'name' in data.labels:
      elems=map(_pythonname,data['name'])
    else:
      elems=['l%05d' % i for i in range(len(data)) ]
    table.__init__(self,elems,col_names,data=data,rowattr=rowattr)

if __name__=='__main__':
  self=tfstable('twiss.ir5b1.data.gz')
  self.show('mqml','s bet')
  self.show((self.s > 4) & (self.s < 10),'s bet')
  print 'tfstable.py:  test OK'
