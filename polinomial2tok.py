from numpy import *
from numpy.linalg import inv 


a=array([[1,4,16],[1,2,4],[1,1,1]])
b=array([2**-2,0.5,1])

a_1=inv(a)

r=dot(a_1,b)

print r
print sum(b)
