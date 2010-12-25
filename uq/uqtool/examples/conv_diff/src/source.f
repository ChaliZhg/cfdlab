c     residual source term
      subroutine source(nc, q, res, xc,dx,a)
      implicit none
      include 'param.h'
      integer :: nc
      real :: q(nc), res(nc),xc(nc),dx,a

      integer :: i
      real :: x,u,ux,uxx,s

      do i=1,nc
        x=xc(i)
        u = 10.0*x*(1.0-x)*sin(a*x)
        ux = 10.0*a*(1.0-x)*x*cos(a*x) + 10.0*(1.0-x)*sin(a*x) -
     1       10.0*x*sin(a*x)
        uxx= -20.0*(a*x*cos(a*x) + sin(a*x)) +
     1        10.0*(1.0-x)*(2.0*a*cos(a*x) - a**2*x*sin(a*x))
        s  = u*ux - uxx

        res(i)=res(i)-s*dx


      enddo

      return
      end
