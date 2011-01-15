c     residual source term
      subroutine source(nc, q, res, xc,xv,a)
      implicit none
      include 'param.h'
      integer :: nc
      real :: q(nc),res(nc),xc(nc),xv(nc+1),a

      integer :: i
      real :: x,u,ux,uxx,s,x1,x2

      do i=1,nc
        x=xc(i)
        x1=xv(i)
        x2=xv(i+1)

        s  = 10*sin(a*x1) - 10*sin(a*x2) - 50*x1**2*sin(a*x1)**2 + 
     <       100*x1**3*sin(a*x1)**2 - 50*x1**4*sin(a*x1)**2 
     <       + 50*x2**2*sin(a*x2)**2 - 100*x2**3*sin(a*x2)**2 + 
     <   50*x2**4*sin(a*x2)**2 - 20*x1*sin(a*x1) + 20*x2*sin(a*x2)
     <    - 10*a*x1**2*cos(a*x1) + 10*a*x2**2*cos(a*x2) 
     <    + 10*a*x1*cos(a*x1) - 10*a*x2*cos(a*x2)
        res(i)=res(i)-s

      enddo

      return
      end
