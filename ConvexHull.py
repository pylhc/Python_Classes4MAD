#!/afs/cern.ch/user/r/rtomas/bin/pythonafs

import sys
from string import split

def orientation(p, q, r):
    ''' >0 if p-q-r are clockwise, <0 if counterclockwise, 0 if colinear. '''
    return (q[1]-p[1])*(r[0]-p[0]) - (q[0]-p[0])*(r[1]-p[1])
def hulls(Points):
    ' Graham scan to find upper and lower convex hulls of a set of 2D points '
    U = [  ]
    L = [  ]
    # the natural sort in Python is lexicographical, by coordinate
    Points.sort( )
    for p in Points:
        while len(U) > 1 and orientation(U[-2], U[-1], p) <= 0:
            U.pop( )
        while len(L) > 1 and orientation(L[-2], L[-1], p) >= 0:
            L.pop( )
        U.append(p)
        L.append(p)
    return U, L





file=sys.argv[1]
col1=int(sys.argv[2])
col2=int(sys.argv[3])
poi=[]
#print "Computing the Convex Hull of", file
for line in open(file,"r"):
    x=float(split(line)[col1-1])
    y=float(split(line)[col2-1])
    #print x,y
    poi.append((x,y))
        



[hup,hdown]= hulls(poi)
#print hup
#print hdown
for p in hup:
    print p[0],p[1]
for i in range(len(hdown),0,-1):
    print hdown[i-1][0], hdown[i-1][1]
