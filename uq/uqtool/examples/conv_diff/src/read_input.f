      subroutine read_input(mode, xi_1, CFL, niter, tol_conv)
      implicit none
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

      end
