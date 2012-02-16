import math
def erf(x):
    # save the sign of x
    sign = 1
    if x < 0: 
        sign = -1
    x = abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
    return sign*y # erf(-x) = -erf(x)


def cdf(x):
    return 0.5*(1+erf(x))



#############

from pylab import *
t = arange(-2, 10.0, 0.1)
s = map(cdf,t)
figure(1)
plot(s,t)
show()

