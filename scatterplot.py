import matplotlib.pyplot as plt
import matplotlib.cm as cm

import random

fig = plt.figure()
fig.subplots_adjust(bottom=0.2)
ax = fig.add_subplot(111)

N=40
delta=[]
vf=[]
dS=[]
for x in range(N):
  for y in range(N):
    delta.append(x)
    vf.append(y)
    dS.append(x*x+y*y)



#plt.scatter(delta,vf,c=dS,alpha=4,cmap=cm.Paired)

plt.scatter(delta,vf,s=80,c=dS,marker='s', edgecolors='none',cmap=cm.jet)
cb = plt.colorbar()
cb.set_label('counts')
plt.axis([0, 40, 0, 40])
plt.ylabel('Y')
plt.xlabel('X')
plt.title('Scatter')
#save('scatter.eps')
plt.show()
