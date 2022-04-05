program Lorenz
!Declaramos variables
real :: k1x,k1y,k1x2, k1y2
real :: k2x,k2y,k2x2, k2y2
real :: k3x,k3y,k3x2, k3y2
real :: k4x,k4y,k4x2, k4y2
real :: x,y,x2,y2
real :: t,dt,tend
real :: alpha, omega, gama, beta, m, k
real :: x0,y0, x02, y02
integer :: n, i

! Asignamos valores de las condiciones iniciales
x0=10.0d0
y0=20.0d0
x=x0
y=y0

! Valor de los parametros
alpha=0.1
k=1.0d0
beta=0.1d0
gama=1.0d0
omega=2.0d0
m=1.0d0

! t es tiempo inicial, tend tiempo final, dt el tamaño del intervalo y n el numero de pasos
t=0.0d0
dt=0.001d0
tend=0.5d0
n=1000d0

dlambda = 0.30d0

x02=x+dlambda
y02=y+dlambda
x2=x02
y2=y02

! Open abre el archivo y si no existe lo crea,
!Almacenamos los valores de x
open(unit=11,file='/home/fer/Documents/Caos/oscilador/seriet.txt')
open(unit=12,file='/home/fer/Documents/Caos/oscilador/seriex.txt')
open(unit=13,file='/home/fer/Documents/Caos/oscilador/seriey.txt')

!Archivos para almacenar los valores x2, y2, z2
open(unit=22,file='/home/fer/Documents/Caos/oscilador/seriex2.txt')
open(unit=23,file='/home/fer/Documents/Caos/oscilador/seriey2.txt')

!Archivos para almacenar el coeficiente de Liapunov
open(unit=32,file='/home/fer/Documents/Caos/oscilador/LiapunovX.dat')
open(unit=33,file='/home/fer/Documents/Caos/oscilador/LiapunovY.dat')

!Iniciamos un ciclo
 do i=1,n
  ! equivalente a: if(t >= tend) ---> termina el programa
  if(t.ge.tend)goto 200

  ! Evaluamos las correspondientes k
  k1x=(y/m)*dt
  k1y=((-m* (beta**2) *x)+ (3*alpha* (x**2) * cos(omega*t) ))*dt

  k2x=((y+ k1x /2)/m)*dt
  k2y=((-m* (beta**2) *(x+ k1x/2))+ (3*alpha* ((x+ k1x/2)**2) * cos(omega*t) ))*dt

  k3x=((y+k2x/2)/m)*dt
  k3y=((-m* (beta**2) *(x+ k2x/2))+ (3*alpha* ((x+ k2x/2)**2) * cos(omega*t) ))*dt

  k4x=((y+k3x/2)/m)*dt
  k4y=((-m* (beta**2) *(x+ k3x/2))+ (3*alpha* ((x+ k3x/2)**2) * cos(omega*t) ))*dt

  x=x+(k1x+2*k2x+2*k3x+k4x)/6
  y=y+(k1y+2*k2y+2*k3y+k4y)/6

  t=t+dt
  Print *, x, y, t

   write(11,98)t
   write(12,98)x
   write(13,98)y
98    format(f21.17)

  k1x2=(y2/m)*dt
  k1y2=((-m* (beta**2) *x2)+ (3*alpha* (x2**2) * cos(omega*t) ))*dt

  k2x2=((y2+ k1x2 /2)/m)*dt
  k2y2=((-m* (beta**2) *(x2+ k1x2/2))+ (3*alpha* ((x2+ k1x2/2)**2) * cos(omega*t) ))*dt

  k3x2=((y2+k2x/2)/m)*dt
  k3y2=((-m* (beta**2) *(x2+ k2x2/2))+ (3*alpha* ((x2+ k2x2/2)**2) * cos(omega*t) ))*dt

  k4x2=((y2+k3x2/2)/m)*dt
  k4y2=((-m* (beta**2) *(x2+ k3x2/2))+ (3*alpha* ((x2+ k3x2/2)**2) * cos(omega*t) ))*dt

  x2=x2+(k1x2+2*k2x2+2*k3x2+k4x2)/6
  y2=y2+(k1y2+2*k2y2+2*k3y2+k4y2)/6

  lambx=log(abs(x2-x)/abs(dlmmbda))/t
  lamby=log(abs(y2-y)/abs(dlambda))/t

  write(22,98)x2
  write(23,98)y2

  write(32,*)t,lambx
  write(33,*)t,lamby

enddo
 200 continue
 stop
 end program

