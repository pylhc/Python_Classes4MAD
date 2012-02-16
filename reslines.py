import matplotlib.pyplot as plt
import matplotlib.cm as cm
from numpy import *
from math import atan2
from optparse import OptionParser



parser = OptionParser()
parser.add_option("-p", "--plot", 
		 help="0 for no  python plot",
		 metavar="PLOT", default="1",dest="plot")
parser.add_option("-N", "--Norder", 
		 help="Maximum resonance order to be shown",
		 metavar="ORDER", default="10",dest="order")
parser.add_option("-b", "--border", 
		 help="x and y range to be displayed xmin,xmax,ymin,ymax",
		 metavar="BORDER", default=".245,.335,.28,.335",dest="border")


(options, args) = parser.parse_args()



def mcd(a,b):
    "Computes the largest common factor of a and b"
    x=max(abs(a),abs(b))
    y=min(abs(a),abs(b))
    if x % y == 0:
      return y
    else:
      return mcd(x,x%y)



def reduce(a,b,c):
  "Divide the 3 integers by the largest common factor..."
  if a==0: # Assume b>0
    f=mcd(b,c)
    return 0,abs(b)/f,c*b/abs(b)/f
  elif b==0: # Assume a>0
    f=mcd(a,c)
    return abs(a)/f, 0, c*a/abs(a)/f
  elif c==0: # Assume a>0
    f=mcd(a,b)
    return abs(a)/f, b*a/abs(a)/f, 0
  else:
    f=mcd(mcd(a,b),c)
    return abs(a)/f, b*a/abs(a)/f, c*a/abs(a)/f

  
  
def Intersect(ray, seg):
  "ray is x,z,angle;  seg is two points in x,y (four numbers)"
  L=2
  segray=[[ray[0],ray[1]],[ray[0]-L*cos(ray[2]), ray[1]-L*sin(ray[2])]]
  #print segray
  x12=seg[0][0] - seg[1][0]
  x34=segray[0][0] - segray[1][0]
  y12=seg[0][1] - seg[1][1]
  y34=segray[0][1] - segray[1][1]
  c=x12 * y34 - y12 * x34
  #print "c",c
  if abs(c)<0.000001:
    #No intersection
    return [0, segray[1][0], segray[1][1]]
  else:
    a = seg[0][0]*seg[1][1] - seg[0][1]*seg[1][0]

    b = segray[0][0]*segray[1][1] - segray[0][1]*segray[1][0]
    x = (a * x34 - b * x12) / c
    y = (a * y34 - b * y12) / c
    xmin = min(seg[0][0], seg[1][0])
    xmax = max(seg[0][0], seg[1][0])
    ymin = min(seg[0][1], seg[1][1])
    ymax = max(seg[0][1], seg[1][1])
    e=0.0001
    if x>= (xmin-e) and x<= (xmax+e) and y>=(ymin-e) and y<=(ymax+e):
      return [1,x,y]
    else:
      return [0,segray[1][0], segray[1][1]]







fig = plt.figure()
fig.subplots_adjust(bottom=0.2)
ax = fig.add_subplot(111)

bb=options.border.split(",")
border=[float(bb[0]),float(bb[1]), float(bb[2]),float(bb[3])] #xmin xmax ymin ymax

order=int(options.order)


seg0=[[border[0],border[2]],[border[1],border[3]]]  #diagonal 1
seg1=[[border[0],border[3]],[border[1],border[2]]]  #diagonal 2

f=open("IntersectedResLines.gplot","w")
print >>f, "set parametric"
print >>f, "set samples 10000"
print >>f, "unset key"
print >>f, "plot [-10:10]["+bb[0]+":"+bb[1]+"]"+"["+bb[2]+":"+bb[3]+"]\\"
AllLines={}

for i in (array(range(2*order+1))-order):
  for j in (array(range(2*order+1))-order):
    cord=abs(i)+abs(j)
    
    if cord<=order:
      for N in (array(range(2*cord+1))-cord):

        
        
          if j==0:  #(vertical line)
             #print i,j,N
             angle=pi/2
             x0=N*1.0/i
             ray=[x0,0,angle]
             xx=[x0,x0]
             yy=[0,1]

          else:
             #print i,j,N
             y0=N*1.0/j              #x0=0
             y1=(N-i)*1.0/j          #x1=1
             angle=atan2(i,-j)
             ray=[0,y0, angle]      
             xx=[0,1]
             yy=[y0,y1]
          
          I0=Intersect(ray,seg0)
          I1=Intersect(ray,seg1)
        
          if I0[0]==1 or I1[0]==1:
            x,y,z=reduce(i,j,N)
            line=str(x)+str(y)+str(z)
            if line not in AllLines.keys():
            
              newcord=abs(x)+abs(y)
              AllLines[line]=1
              if y!=0:
                print >>f, "t,("+str(z)+"-("+str(x)+"*t))/"+str(y)+" lt 3 lw "+str(order-newcord+1)+" , \\"
              else:
                print >>f, "("+str(z)+"-("+str(y)+"*t))/"+str(x)+", t lt 3 lw "+str(order-newcord+1)+", \\"
              plt.plot(xx,yy,c='b',linewidth=order-newcord+1)

print >> f, "0,0 lt 0 lw 0"
f.close()

        
plt.scatter([0.28,0.31],[0.31,0.32],c='r')


plt.axis(border)
plt.ylabel('Y')
plt.xlabel('X')

#save('scatter.eps')
if options.plot=="1":
  plt.show()
