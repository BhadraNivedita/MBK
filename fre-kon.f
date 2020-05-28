	! integrates the frenkel-kontorova model by rk4 method

	implicit none

	integer, parameter :: nn=20,mm=nn/2,iter=1000
	double precision   :: xx(nn),yy(nn),dy(nn)
	double precision   :: tau,dtau
	integer            :: i,j,k,l,m,n
	double precision   :: Lam,KK,g0,g1,gam,pi

	common Lam,KK,g0,g1,gam

	open (unit=5,file='fre-kon.dat',status='unknown')
	
	Lam=1.0d0; KK=1.50d0
 	g0=10.0d0; g1=320.0d0; gam=40.0d0; pi=3.1416;

	tau=0.0d0; dtau=0.01d0
	

	yy(1)=pi+0.002d0
	yy(2)=pi+0.003d0
	yy(3)=pi+0.004d0
	yy(4)=pi+0.005d0
	yy(5)=pi+0.006d0
	yy(6)=pi+0.007d0
	yy(7)=pi+0.008d0
	yy(8)=pi+0.009d0
	yy(9)=pi+0.0023d0
	yy(10)=pi+0.0024d0
		
	
	yy(11)=0.001d0
	yy(12)=0.001d0
	yy(13)=0.001d0
	yy(14)=0.001d0
	yy(15)=0.001d0
	yy(16)=0.001d0
	yy(17)=0.001d0
	yy(18)=0.001d0
	yy(19)=0.001d0
	yy(20)=0.001d0
	
	
	do i=1,iter	! iterations

	call derivs(tau,yy,dy)
        call rk4(yy,dy,nn,tau,dtau,xx,derivs)
	yy=xx
	tau = tau + dtau
        write(5,'(1000(xf12.6))') (yy(l),l=1,nn)
	end do		! iterations

	end

	! calculates the time derivative 

        subroutine derivs(tau,yy,dy)

        integer, parameter :: nn=20,mm=nn/2
	double precision   :: yy(nn),dy(nn)
	double precision   :: tau,dtau
	integer            :: i,j,k,l,m,n
        double precision   :: Lam,KK,g0,g1,gam

	common Lam,KK,g0,g1,gam

        do i=1,mm
	dy(i) = Lam*KK*yy(mm+i)
        end do

	dy(mm+1) = (Lam/KK)*(dsin(yy(mm)-yy(1))-dsin(yy(1)-yy(2)))
     .                 -(g0+g1*dcos(gam*tau))*dsin(yy(1))/Lam

	do i= 2,mm-1
	dy(mm+i) = (Lam/KK)*(dsin(yy(i-1)-yy(i))-dsin(yy(i)-yy(i+1)))
     .                 -(g0+g1*dcos(gam*tau))*dsin(yy(i))/Lam
	end do

	dy(nn) = (Lam/KK)*(dsin(yy(mm-1)-yy(mm))-dsin(yy(mm)-yy(1)))
     .                 -(g0+g1*dcos(gam*tau))*dsin(yy(mm))/Lam

        print *, 'The velocity on impact is', dy(nn)
     
        return
        end subroutine

	! performs the rk4 integration

	SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs)

	INTEGER n,NMAX
	REAL*8 h,x,dydx(n),y(n),yout(n)
	EXTERNAL derivs
	PARAMETER (NMAX=2000)
	INTEGER i
	REAL*8 h6,hh,xh,dym(NMAX),dyt(NMAX),yt(NMAX)

	hh=h*0.5
	h6=h/6.
	xh=x+hh

	do i=1,n
        yt(i)=y(i)+hh*dydx(i)
	end do

	call derivs(xh,yt,dyt)
	do i=1,n
        yt(i)=y(i)+hh*dyt(i)
	end do

	call derivs(xh,yt,dym)
	do i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
	end do

	call derivs(x+h,yt,dyt)
	do i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i))
	end do

	return
	END
