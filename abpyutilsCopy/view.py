class view(object):
  _addr={}
  def __init__(self,data,**addr):
    view._addr[self]=addr
    self._addr=addr
    self._data=data
    for k in addr:
      self.__dict__[k]=None
  def __getattribute__(self,k):
    d=view._addr[self]
    if k in d:
      return self._data[d[k]]
    else:
      return object.__getattribute__(self,k)
  def __setattr__(self,k,v):
    d=view._addr[self]
    if k in d:
      self._data[d[k]]=v
    object.__setattr__(self,k,v)

  def _clear(self,data,**addr):
    self.__dict__.clear()
    view._addr[self]=addr
    self._addr=addr
    self._data=data
    for k in addr:
      self.__dict__[k]=None

def test():
  print 'Test view'
  from numpy import zeros,array
  data=zeros((10,10))
  a=view(data,x=1,y=(0,1),z=(slice(1,5),slice(1,2)))
  a.normattr='test'
  a.x=12
  a.y=32
  a.z.T[0]=[1,2,3,4]
  t=[]
  t.append(data[0,1]==32)
  t.append(data[1,1]==1)
  t.append(data[1,0]==12)
  t.append(a.normattr=='test')
  res=reduce(lambda x,y: x and y,t)
  if not res: print 'Error in view'
  return res
