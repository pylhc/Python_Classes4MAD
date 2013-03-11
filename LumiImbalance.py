import sys
from math import sqrt
def imba(b):
    """ Needs a vector with 16 components:
    B1HIP1
    B1VIP1
    B2HIP1
    B2VIP1      
    B1HIP5
    B1VIP5
    B2HIP5
    B2VIP5
    and its errors. Outputs lumi ratio, relative error
    and equivalent betas at IP1 and IP5"""

    aveH1=(b[0]+b[4])/2.
    aveV1=(b[2]+b[6])/2.
    IP1=sqrt(1./aveH1/aveV1)
    aveH5=(b[8]+b[12])/2.
    aveV5=(b[10]+b[14])/2.
    IP5=sqrt(1./aveH5/aveV5)
    imba=IP5/IP1
    eH1=0.5*sqrt(b[1]**2+b[5]**2)  #absolute error
    eV1=0.5*sqrt(b[3]**2+b[7]**2)
    e1=0.5*sqrt(eH1**2/aveH1**2+eV1**2/aveV1**2) # relative err
    eH5=0.5*sqrt(b[9]**2+b[13]**2)
    eV5=0.5*sqrt(b[11]**2+b[15]**2)
    e5=0.5*sqrt(eH5**2/aveH5**2+eV5**2/aveV5**2)
    rel=e1+e5
    #print ea1, eb1, ea5,eb5, e1, e5
    return imba, rel, 1./IP1,e1, 1./IP5,e5

x=imba([0.71, 0.09,   
      0.67, 0.06,   
      0.72, 0.10,  
      0.71, 0.10,  
      0.66, 0.07,  
      0.65, 0.07,  
      0.73, 0.09,  
      0.68, 0.07]) 



print "K-mod (ratio, relative error)", x


x=imba([0.60, 0.06,
 0.61, 0.02,
 0.59, 0.03,
 0.58, 0.02,
		
 0.60, 0.05,
 0.58, 0.07,
 0.60, 0.02,
 0.58, 0.04])

print "AC dip", x

sys.exit()
"""
 0.71 0.09   
 0.67 0.06   
 0.72 0.10  
 0.71 0.10  

 0.66 0.07  
 0.65 0.07  
 0.73 0.09  
 0.68 0.07
"""
