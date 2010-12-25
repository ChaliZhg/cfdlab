c     residual source term
      subroutine source(nc, q, res, xc,dx,a)
      implicit none
      include 'param.h'
      integer :: nc
      real :: q(nc), res(nc),xc(nc),dx,a

      integer :: i
      real :: s,x

      do i=1,nc
        x=xc(i)
        s = -2.0*a*(1.0-2.0*x)*cos(a*x) + 2.0*sin(a*x) + 
     1      a**2 * (1.0-x)*x*sin(a*x) + 
     2      (1.0-x)*x*sin(a*x)*(a*(1.0-x)*x*cos(a*x) +
     3      (1.0-x)*sin(a*x) - x*sin(a*x))
        s = 10.0*s

        res(i)=res(i)-s*dx


      enddo

      return
      end
