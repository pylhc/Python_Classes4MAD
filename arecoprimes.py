def gcd(a,b):
    while b != 0:
        a, b = b, a%b
    return a



def coprime(a,b):
    r=0.
    if gcd(a,b) == 1:
        r=1.
    return r



for N in range(4,51):
    t=4*(N-2)*(N-1)+4*(N-4)*(N-2)
    for y in range(3,N/2+1):
        for x in range(1,y/2+1):
            t=t+8*coprime(x,y)*(N-2*y)*(N-y)
    print N, t, t%8.


import sys
sys.exit()


for N in range(4,51):
    t=4*(N-2)*(N-1)
    for y in range(1,N/2+1):
        for x in range(1,N/2+1):
            
            if x!=y:
                t=t+2*coprime(x,y)*(N-2*y)*(N-2*x)
    print N, t, t%8.
    
