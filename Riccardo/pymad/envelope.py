from numpy import zeros


class envelope:
  "find aperture give a survey and a twiss table"

#    kbeta=1.1,      # beta beating
#    nsigma=9.5,      # 9.5
#    emit=3.75E-6,   #3.75E-6
#    delta=1.129E-4, # RMS energy spread
#    tol=4.6E-3,  # CO=3mm + dtol=1.6mm
#    deltamax=8E-4,   # for chromaticity measurment
#    betamaxarc=180,  # maximum beta in the arcs
#    dxmaxarc=2,  # maximum beta in the arc
def __init__(s,t, kbeta=1.1, nsigma=9.5, emit=3.75E-6, delta=1.129E-4, tol=4.6E-3, deltamax=8E-4, betamaxarc=180, dxmaxarc=2):
  self.co=zeros([len(s1.x),3],float)
  self.xp=zeros([len(s1.x),3],float)
  self.xm=zeros([len(s1.x),3],float)
  self.yp=zeros([len(s1.x),3],float)
  self.ym=zeros([len(s1.x),3],float)
  self.yp2D=zeros([len(s1.x),3],float)
  self.ym2D=zeros([len(s1.x),3],float)

  for i in xrange(len(s1.x)):
    vro = array([s.x[i],s.y[i],s.z[i]])
    theta,phi,psi = s.theta[i],s.phi[i],s.psi[i]
    betx,bety,dx,dy,x,y= t.betx[i],t.bety[i],t.dx[i],t.dy[i],t.x[i],t.y[i]
    thetam=array([[cos(theta) ,           0,sin(theta)],
         [          0,           1,         0],
         [-sin(theta),           0,cos(theta)]])
    phim=  array([[          1,          0,          0],
        [          0,cos(phi)   ,   sin(phi)],
        [          0,-sin(phi)  ,   cos(phi)]])
    psim=  array([[   cos(psi),  -sin(psi),          0],
        [   sin(psi),   cos(psi),          0],
        [          0,          0,          1]])
    wm=matrixmultiply(thetam,matrixmultiply(phim,psim))
    ex=matrixmultiply(wm,array([1,0,0]))
    ey=matrixmultiply(wm,array([0,1,0]))
    self.co[i]=vro+x * ex + y * ey
    dx+= dxmaxarc*sqrt(bx/betamaxarc)
    dy+= dxmaxarc*sqrt(by/betamaxarc)
    xsize=kbeta* (nsigma*sqrt(betx*emit + (dx*delta)**2) + deltamax*dx)+ tol
    ysize=kbeta* (nsigma*sqrt(bety*emit + (dy*delta)**2) + deltamax*dx)+ tol
    self.xp[i]=self.co[i] + xsize * ex
    self.xm[i]=self.co[i] - xsize * ex
    self.yp[i]=self.co[i] + ysize * ey
    self.ym[i]=self.co[i] - ysize * ey
    self.yp2D[i]=self.co[i] + ysize * ex
    self.ym2D[i]=self.co[i] - ysize * ex

if __name__ == '__main__':

  s1=madtable('survey.lhcb1.data')
  s2=madtable('survey.lhcb2.data')
  t1=madtable('twiss.lhcb1.data')
  t2=madtable('twiss.lhcb2.data')
  ap1=aperture(s1,t1)
  ap2=aperture(s1,t1)

  hold(False)
  figure(figsize=(6,6))
  hold(True)
  plot(ap1.co[:,2],ap1.co[:,0])
  plot(ap1.xp[:,2],ap1.xp[:,0],'g',linewidth=.1)
  plot(ap1.yp2D[:,2],ap1.yp2D[:,0],'r',linewidth=.1)
  plot(ap1.xm[:,2],ap1.xm[:,0],'g',linewidth=.1)
  plot(ap1.ym2D[:,2],ap1.ym2D[:,0],'r',linewidth=.1)
  axis([-10000,10000,-15000,5000])

  savefig('ring.eps',dpi=600)

  plot(s2.z,s2.x)

