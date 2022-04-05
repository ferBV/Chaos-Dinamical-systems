program Lorenz

 real, parameter :: pi=3.141592d0
 real :: k1x,k1y,k1z
 real :: k2x,k2y,k2z
 real :: k3x,k3y,k3z
 real :: k4x,k4y,k4z
 real :: x,y,z
 real :: t,dt,tend
 real :: sigma,b,rc,rho
 real :: xo,yo,zo
 real :: tn1,tn
 real :: k1x2,k1y2,k1z2
 real :: k2x2,k2y2,k2z2
 real :: k3x2,k3y2,k3z2
 real :: k4x2,k4y2,k4z2
 real :: x2,y2,z2
 real :: x02,y02,z02
 integer :: n, i


 n=10000000d0
 xo=1.0d0
 yo=0.0d0
 zo=0.0d0
 x=xo
 y=yo
 z=zo

 sigma=10.0d0
 b=8./3.
 rc=sigma*(sigma+b+3.)/(sigma-b-1.)
 rho=13.0d0

 t=0.0d0
 dt=0.001d0
 tend=50.0d0

 !write(*,*)sigma,b,rc,rho
!!!!!!!!!!!!!!!


!print *,tn, xn, vn
x02=x+0.01
y02=y+0.5
z02=z+0.9
x2=x02
y2=y02
z2=z02

!d0x=x0-x02
!d0y=y0-y02
!d0z=z0-z02


!!!!!!!!!!!!!!!!!!!!

open(unit=64,file='/home/fer/Documents/Caos/lorenz/seriex2.txt')
open(unit=65,file='/home/fer/Documents/Caos/lorenz/seriey2.txt')
open(unit=66,file='/home/fer/Documents/Caos/lorenz/seriez2.txt')

open(unit=67,file='/home/fer/Documents/Caos/lorenz/LiapunovX.dat')
open(unit=68,file='/home/fer/Documents/Caos/lorenz/LiapunovY.dat')
open(unit=69,file='/home/fer/Documents/Caos/lorenz/LiapunovZ.dat')


!!!!!!!!!!!!!!!!!!!!!


 do i=1,n

  if(t.ge.tend)goto 200
    !if(i.le.2)write(*,*)t,x,y,z
  k1x=-sigma*(x-y)*dt
  k1y=(rho*x-y-z*x)*dt
  k1z=(-b*z+x*y)*dt

       !if(i.le.n)write(*,*)t,k1x,k1y,k1z

  k2x=-sigma*(x+k1x/2.-y-k1y/2.)*dt
  k2y=(rho*(x+k1x/2.)-y-k1y/2.-(z+k1z/2.)*(x+k1x/2.))*dt
  k2z=(-b*(z+k1z/2.)+(x+k1x/2.)*(y+k1y/2.))*dt

  k3x=-sigma*(x+k2x/2.-y-k2y/2.)*dt
  k3y=(rho*(x+k2x/2.)-y-k2y/2.-(z+k2z/2.)*(x+k2x/2.))*dt
  k3z=(-b*(z+k2z/2.)+(x+k2x/2.)*(y+k2y/2.))*dt

  k4x=-sigma*(x+k3x-y-k3y)*dt
  k4y=(rho*(x+k3x)-y-k3y-(z+k3z)*(x+k3x))*dt
  k4z=(-b*(z+k3z)+(x+k3x)*(y+k3y))*dt

  x=x+(k1x+2.*k2x+2.*k3x+k4x)/6.
  y=y+(k1y+2.*k2y+2.*k3y+k4y)/6.
  z=z+(k1z+2.*k2z+2.*k3z+k4z)/6.
  t=t+dt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  k1x2=(-sigma*x2+sigma*y2)*dt
  k1y2=(rho*x2-y2-z2*x2)*dt
  k1z2=(-b*z2+x2*y2)*dt

  k2x2=(-sigma*(x2+k1x2/2)+sigma*(y2+k1y2/2))*dt
  k2y2=(rho*(x2+k1x2/2)-(y2+k1y2/2)-(z2+k1z2/2)*(x2+k1x2/2))*dt
  k2z2=(-b*(z2+k1z2/2)+(x2+k1x2/2)*(y2+k1y2/2))*dt

  k3x2=(-sigma*(x2+k2x2/2)+sigma*(y2+k2y2/2))*dt
  k3y2=(rho*(x2+k2x2/2)-(y2+k2y2/2)-(z2+k2z2/2)*(x2+k2x2/2))*dt
  k3z2=(-b*(z2+k2z2/2)+(x2+k2x2/2)*(y2+k2y2/2))*dt

  k4x2=(-sigma*(x2+k3x2)+sigma*(y2+k3y2))*dt
  k4y2=(rho*(x2+k3x2)-(y2+k3y2)-(z2+k3z2)*(x2+k3x2))*dt
  k4z2=(-b*(z2+k3z2)+(x2+k3x2)*(y2+k3y2))*dt

  x2=x2+(k1x2+2*k2x2+2*k3x2+k4x2)/6
  y2=y2+(k1y2+2*k2y2+2*k3y2+k4y2)/6
  z2=z2+(k1z2+2*k2z2+2*k3z2+k4z2)/6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  lambx=log(abs(x2-x)/abs(0.5))/t
  lamby=log(abs(y2-y)/abs(0.3))/t
  lambz=log(abs(z2-z)/abs(0.9))/t

!!!!!!!!!!!!!!!!
      open(unit=60,file='/home/fer/Documents/Caos/lorenz/seriet.txt')
      open(unit=61,file='/home/fer/Documents/Caos/lorenz/seriex.txt')
      open(unit=62,file='/home/fer/Documents/Caos/lorenz/seriey.txt')
      open(unit=63,file='/home/fer/Documents/Caos/lorenz/seriez.txt')
       write(60,98)t
       write(61,98)x
       write(62,98)y
       write(63,98)z
98    format(f21.17)

!!!!!!!!!!!!!!!!!

      write(64,98)x2
      write(65,98)y2
      write(66,98)z2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    open(unit=80,file='/home/fer/Documents/Caos/lorenz/poincare.txt')
    if(z<28.20.and.z>27.80)write(80,*)x,y

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  write(67,*)t,lambx
  write(68,*)t,lamby
  write(69,*)t,lambz




enddo


 200 continue
 stop
 end program

