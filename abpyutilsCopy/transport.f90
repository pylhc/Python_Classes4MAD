function colname()
character(len=200) :: colname
colname="l kn0l ks0l kn1l ks1l br gr &
& s x px y py t pt dx dpx dy dpy betx alfx mux bety alfy muy"
end function colname

subroutine quad(v,n,r)
! Compute the linear transfer matrix of a
! combined fucntion magnet
implicit none

double precision,intent(in) :: v(n)
integer,intent(in) :: n
double precision,intent(out):: r(6,6)

double precision ::l,kn0l,ks0l,kn1l,ks1l,br,gr
double precision :: cx,sx,cy,sy,dx,dy,kx,ky,skx,sky,fx

l  =v(1)
kn0l=v(2)
ks0l=v(3)
kn1l=v(4)
ks1l=v(5)
br =v(6)
gr =v(7)

r=0

if (l.eq.0d0) then
  r(1,1)=1
  r(1,2)=0
  r(2,1)=-kn1l-kn0l**2
  r(2,2)=1
  r(3,3)=1
  r(3,4)=0
  r(4,3)=kn1l-ks0l**2
  r(4,4)=1
  r(1,6)=0 !to be checked
  r(2,6)=kn0l ! to be checked
  r(5,1)=-r(2,6)
  r(5,2)=-r(1,6)
  r(3,6)=0 !to be checked
  r(4,6)=ks0l ! to be checked
  r(5,3)=-r(4,6)
  r(5,4)=-r(3,6)
  r(5,5)=1
  r(5,6)=l/br**2/gr**2
  r(6,6)=1
else
  kx=(+kn1l/l+(kn0l/l)**2)
  if (abs(kx).lt.1E-10) then
    cx=1-l**2*kx/2
    sx=l-l**3*kx/6
    dx=l**2/2
    fx=l**3/6
  else
    if (kx.gt.0D0) then
      skx=sqrt(kx)
      cx=cos(skx*l)
      sx=sin(skx*l)/skx
      dx =(1-cx)/kx
      fx =(l-sx)/kx
    else
      skx=sqrt(-kx)
      cx=cosh(skx*l)
      sx=sinh(skx*l)/skx
      dx =(1-cx)/kx
      fx =(l-sx)/kx
    endif
  endif

  ky=(-kn1l/l+(ks0l/l)**2)
  if (abs(ky).lt.1E-10) then
    cy=1-l**2*ky/2
    sy=l-l**3*ky/6
    dy=l**2/2
  else
    if (ky.gt.0D0) then
      sky=sqrt(ky)
      cy=cos(sky*l)
      sy=sin(sky*l)/sky
      dy =(1-cy)/ky
    else
      sky=sqrt(-ky)
      cy=cosh(sky*l)
      sy=sinh(sky*l)/sky
      dy =(1-cy)/ky
    endif
  endif

  r(1,1)=cx
  r(1,2)=sx
  r(2,1)=-kx*sx
  r(2,2)=cx
  r(3,3)=cy
  r(3,4)=sy
  r(4,3)=-ky*sy
  r(4,4)=cy
  r(1,6)=dx*kn0l/l
  r(2,6)=sx*kn0l/l
  r(5,1)=-r(2,6)
  r(5,2)=-r(1,6)
  r(3,6)=dy*ks0l/l
  r(4,6)=sy*ks0l/l
  r(5,3)=-r(4,6)
  r(5,4)=-r(3,6)
  r(5,5)=1
  r(5,6)=l/br**2/gr**2 -(kn0l/l)**2 * fx / br**2
  r(6,6)=1

endif

end subroutine quad


subroutine loop(v,m,n,o)
implicit none

double precision,intent(inout) :: v(n+o,m)
integer,intent(in) :: n,m,o

double precision:: z(6)
double precision:: betx1,alfx1,mux1,bety1,alfy1,muy1
double precision:: betx2,alfx2,mux2,bety2,alfy2,muy2
double precision:: r11,r12,r21,r22,tmp1,tmp2
!double precision:: kq
double precision:: r(6,6),vt(n)
integer :: i


do i=1,m-1
  vt=v(1:n,i)
!  write(*,*) vt
  call quad(vt,n,r)
!  write(*,'(6e12.4)') r
!  write(*,*)

  ! Track s
  v(n+1,i+1)=v(n+1,i)+v(1,i)

  ! Track orbit
  v(n+2:n+7,i+1)=matmul(r,v(n+2:n+7,i))
!  ! Track dispersion
  z(1:4)=v(n+8:n+11,i)
  z(5)=0D0
  z(6)=1D0
  z=matmul(r,z)
  v(n+8:n+11,i+1)=z(1:4)
!  ! Track Twiss
  betx1 =v(n+12,i)
  alfx1 =v(n+13,i)
  mux1  =v(n+14,i)
  bety1 =v(n+15,i)
  alfy1 =v(n+16,i)
  muy1  =v(n+17,i)
!  kq=abs(v(4,i))
!  v(n+18,i)=betx1+alfx1**2/betx1/abs(kq)
!  v(n+19,i)=bety1+alfy1**2/bety1/abs(kq)
  r11=r(1,1)
  r12=r(1,2)
  r21=r(2,1)
  r22=r(2,2)
  tmp1  =( r11*betx1 - r12*alfx1 )
  tmp2  =( r21*betx1 - r22*alfx1 )
  betx2 =( tmp1**2 + r12**2 ) / betx1
  alfx2=-( tmp1*tmp2 + r12*r22 ) / betx1
  mux2=mux1 + atan2( r12,tmp1 ) / 8/atan(1D0)
  r11=r(3,3)
  r12=r(3,4)
  r21=r(4,3)
  r22=r(4,4)
  tmp1  =( r11*bety1 - r12*alfy1 )
  tmp2  =( r21*bety1 - r22*alfy1 )
  bety2 =( tmp1**2 + r12**2 ) / bety1
  alfy2=-( tmp1*tmp2 + r12*r22 ) / bety1
  muy2=muy1 + atan2( r12,tmp1 ) / 8/atan(1D0)
  v(n+12,i+1) =betx2
  v(n+13,i+1) =alfx2
  v(n+14,i+1) =mux2
  v(n+15,i+1) =bety2
  v(n+16,i+1) =alfy2
  v(n+17,i+1) =muy2
end do

end subroutine loop


