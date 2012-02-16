import pylab as _p
import numpy as _n
from utils import container

lglabel={
    'betx':    r'$\beta_x$',
    'bety':    r'$\beta_y$',
    'dx':    r'$D_x [m]$',
    'dy':    r'$D_y [m]$',
    }

axlabel={
    's':       r'$s [m]$',
    'betx':    r'$\beta [m]$',
    'bety':    r'$\beta [m]$',
    'dx':    r'$D [m]$',
    'dy':    r'$D [m]$',
    'x':    r'$co [m]$',
    'y':    r'$co [m]$',
    }


def _mylbl(d,x): return d.get(x,r'$%s$'%x)


class qdplot(container):
  def __repr__(self):
    return object.__repr__(self)
  def update(self):
    for i in self.plots:
      try:
       i.set_ydata(i.get_ydata())
       i.set_xdata(i.get_xdata())
      except:
        pass
    _p.draw()



class plotfunc(object):
  _is_s_begin=True

  def plot(self,yl='',yr='',x='s',idx=slice(None),clist='k r b g c m',lattice=True,newfig=True):
    self._clist=clist.split()
    if newfig:
      f=_p.figure()
    else:
      f=_p.gcf()
    sp=_p.subplot(111)
    _p.subplots_adjust(right=0.72)
    _p.setp( sp.yaxis, visible=False)
    xd=getattr(self,x)[idx]
    out=qdplot(figure=f,subplot=sp,plots=[],lines=[],legends=[],xaxis=xd)
    if lattice:
      self._lattice(out,idx,['kn0l','angle'],"#a0ffa0",'Bend h')
      self._lattice(out,idx,['ks0l'],"#ffa0a0",'Bend v')
      self._lattice(out,idx,['kn1l','k1l'],"#a0a0ff",'Quad')
      self._lattice(out,idx,['hkick'],"#e0a0e0",'Kick h')
      self._lattice(out,idx,['vkick'],"#a0e0e0",'Kick v')
    for i in yl.split():
      self._column(out,idx,i,'left')
    for i in yr.split():
      self._column(out,idx,i,'right')

    _p.xlabel(_mylbl(axlabel,x))

    _p.xlim(min(xd),max(xd))
    _p.figlegend(out.lines,out.legends,'upper right')
    _p.grid()
    _p.draw()
    del self._clist
    if hasattr(self,'_colAxleft'): delattr(self,'_colAxleft')
    if hasattr(self,'_colAxright'): delattr(self,'_colAxright')
    self.currplot=out
    return out


  def _lattice(self,out,idx,names,color,lbl):
    vd=0
    sp,s=out.subplot,out.xaxis
    for i in names:
      if hasattr(self,i): vd=getattr(self,i)[idx]+vd
    if vd is not 0:
      m=_n.abs(vd).max()
      if m>1E-10:
        c=_n.where(abs(vd) > m*1E-4)[0]
        if len(c)>0:
          if _n.all(self.l[c]>0):
            vd[c]=vd[c]/self.l[c]
            m=abs(vd[c]).max()
          vd[c]/=m
          if self._is_s_begin:
            plt=_p.bar(s[c],vd[c],self.l[c],color='#aaffaa')
          else:
            plt=_p.bar(s[c]-self.l[c],vd[c],self.l[c],color='#aaffaa')
  #        _p.setp(plt,facecolor=color,edgecolor="#666666")
          _p.setp(plt,facecolor=color,edgecolor=color)
          if plt:
            out.plots.append(plt)
            out.lines.append(plt[0])
            out.legends.append(lbl)
        _p.ylim(-1.5,1.5)
        _p.setp( sp.yaxis, visible=False)

  def _column(self,out,idx,name,pos='left'):
    fig,sp,s=out.figure,out.subplot,out.xaxis
    if hasattr(self,'_colAx'+pos):
      Ax=getattr(self,'_colAx'+pos)
    else:
      Ax=fig.add_axes(sp.get_position(), sharex=sp, frameon=False, label=name)
      setattr(self,'_colAx'+pos,Ax)
    color=self._clist.pop(0)
    self._clist.append(color)
    bxp,=Ax.plot(s,getattr(self,name)[idx],color,label=_mylbl(lglabel,name))
    Ax.yaxis.set_label_position(pos)
    Ax.yaxis.set_ticks_position(pos)
    _p.ylabel(_mylbl(axlabel,name))
    out.plots.append(bxp)
    out.lines.append(bxp)
    out.legends.append(_mylbl(lglabel,name))

  def plotbeta(self,newfig=True):
    self.plot('betx bety','dx dy',newfig=newfig)

  def plotcross(self,newfig=True):
    self.plot('x y','dx dy',newfig=newfig)

  def plottune(self,lbl=''):
    _p.title(r"${\rm Tune} \quad {\rm vs} \delta$")
    _p.xlabel("$\delta$")
    _p.ylabel("Fractional tune")
    tt=r'$%s \rm{%s}$'
    _p.plot(self.deltap,self.q1-self.q1.round(),label=tt %('Q_x',lbl))
    _p.plot(self.deltap,self.q2-self.q2.round(),label=tt %('Q_y',lbl))
    qx=(self.q1-self.q1.round())[abs(self.deltap)<1E-15][0]
    qy=(self.q2-self.q2.round())[abs(self.deltap)<1E-15][0]
    _p.text(0.0,qx,r"$Q_x$")
    _p.text(0.0,qy,r"$Q_y$")
    _p.grid(True)
    _p.legend()

  def plotbetabeat(self,t1,dp='0.0003'):
    _p.title(r"$\rm{Beta beat: 1 - \beta(\delta=%s)/\beta(\delta=0)}$" % dp)
    _p.ylabel(r"$\Delta\beta/\beta$")
    _p.xlabel(r"$s$")
    _p.plot(self.s,1-t1.betx/self.betx,label=r'$\Delta\beta_x/\beta_x$')
    _p.plot(self.s,1-t1.bety/self.bety,label=r'$\Delta\beta_y/\beta_y$')
    _p.grid(True)
    _p.legend()

  def plotw(self,lbl=''):
    title(r"Chromatic function: %s"%lbl)
  #  ylabel(r"$w=(\Delta\beta/\beta)/\delta$")
    ylabel(r"$w$")
    xlabel(r"$s$")
    plot(self.s,self.wx,label=r'$w_x$')
    plot(self.s,self.wy,label=r'$w_y$')
    grid(True)
    legend()

