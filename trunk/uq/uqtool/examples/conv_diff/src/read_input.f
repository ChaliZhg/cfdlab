      subroutine read_input(mode, xi_1, CFL, niter, tol_conv)
      implicit none
      include 'param.h'

      integer :: mode, niter, i
      real    :: xi_1, CFL,tol_conv

      open(10, file='param.in', status='old')
      read(10,*)mode
      read(10,*)xi_1 
      read(10,*)CFL
      read(10,*)niter
      read(10,*)tol_conv
      close(10)
      print*,'Read  input from param.in'

      ark(1) = 0.0
      ark(2) = 3.0/4.0
      ark(3) = 1.0/3.0

      do i=1,3
         brk(i) = 1.0 - ark(i)
      enddo

      end
