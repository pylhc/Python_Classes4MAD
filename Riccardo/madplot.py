from pylab import *
from utils import *

class madplot:
  """
use after:
select,flag=twiss,column=name,s,l,betx,bety,dx,k0l,k1l;
use,sequence=lhcb1,range=IP5/e.DS.R5.B1;
twiss,betx=0.25,bety=0.25,file=twiss.data;

t=madtable('twiss.data')
plotbeta(t)
  """
  def __init__(t,ax,ay1,ay2=[],strs=[],thin=False):
    self.t=t
    self.ax=ax
    self.ay1=ay1
    self.ay2=ay1
    self.strs=strs
    f=figure()
    elemAxes=subplot(111)
    subplots_adjust(right=0.75)
    if thin:
      mask=t.nrange(1E-4,1,abs(t.angle))
      tk0=t.extract(mask,[t.s,1E-2,t.angle,t.k1l])
      mask=t.nrange(1E-4,1,abs(t.k1l))
      tk1=t.extract(mask,[t.s,1E-2,t.k1l,t.k1l])
    else:
      mask=t.nrange(1E-4,1,abs(t.angle)) & t.nrange(1E-4,1,abs(t.angle))
      tk0=t.extract(mask,[t.s-t.l,t.l,t.angle/t.l,t.k1l/t.l])
      mask=t.nrange(1E-4,1,abs(t.k1l)) & t.nrange(1E-4,1,abs(t.k1l))
      tk1=t.extract(mask,[t.s-t.l,t.l,t.k1l/t.l,t.k1l/t.l])

    k0p=bar(tk0[:,0],tk0[:,2],tk0[:,1],color='#aaffaa')
    setp(k0p,facecolor="#e0ffe0",edgecolor="#666666")
    k1p=bar(tk1[:,0],tk1[:,2],tk1[:,1],color='#aaaaff')
    setp(k1p,facecolor="#e0e0ff",edgecolor="#666666")
    setp( elemAxes.yaxis, visible=False)

    betaAxes=gcf().add_axes(elemAxes.get_position(), sharex=elemAxes, frameon=False, label="beta")
    bxp,=plot(t.s,sqrt(t.betx),'k',label=r'$\beta_x$')
    byp,=plot(t.s,sqrt(t.bety),'r',label=r'$\beta_y$')
    axis(ymax=200)
    betaAxes.yaxis.set_label_position('left')
    betaAxes.yaxis.set_ticks_position('left')
    xlabel(r'$s [m]$')
    ylabel(r'$\sqrt\beta [m]$')

    dispAxes=gcf().add_axes(elemAxes.get_position(), sharex=betaAxes, frameon=False, label="disp")
    dispAxes.yaxis.tick_right()
    dispAxes.yaxis.set_label_position('right')
    dispAxes.yaxis.set_ticks_position('right')
    dxp,=plot(t.s,t.dx,'g',label=r'$D_x$')
    ylabel(r"$D [m]$")


    figlegend([bxp,byp,dxp,k0p[0],k1p[0]],[r'$\beta_x$',r'$\beta_y$',r'$D_x$','Bend','Quad'],'upper right')
    draw()


  def update(plot):
    self.t.update()
    self.bxp.set_ydata(sqrt(t.betx))
    self.byp.set_ydata(sqrt(t.bety))
    self.dxp.set_ydata(t.dx)
    self.bxp.set_xdata(t.s)
    self.byp.set_xdata(t.s)
    self.dxp.set_xdata(t.s)

    if thin:
      mask=t.nrange(1E-4,1,abs(t.angle))
      tk0=t.extract(mask,[t.s,1E-2,t.angle,t.k1l])
      mask=t.nrange(1E-4,1,abs(t.k1l))
      tk1=t.extract(mask,[t.s,1E-2,t.k1l,t.k1l])
    else:
      mask=t.nrange(1E-4,1,abs(t.angle)) & t.nrange(1E-4,1,abs(t.angle))
      tk0=t.extract(mask,[t.s-t.l,t.l,t.angle/t.l,t.k1l/t.l])
      mask=t.nrange(1E-4,1,abs(t.k1l)) & t.nrange(1E-4,1,abs(t.k1l))
      tk1=t.extract(mask,[t.s-t.l,t.l,t.k1l/t.l,t.k1l/t.l])

    for i in arange(len(k0p)):
      k0p[i].set_x(float(tk0[i,0]))
      k0p[i].set_width(float(tk0[i,1]))
      k0p[i].set_y(0)
      k0p[i].set_height(float(tk0[i,2]))
    for i in arange(len(k1p)):
      k1p[i].set_x(float(tk1[i,0]))
      k1p[i].set_width(float(tk1[i,1]))
      k1p[i].set_y(0)
      k1p[i].set_height(float(tk1[i,2]))
    draw()
    return (t,bxp,byp,dxp,k0p,k1p)



