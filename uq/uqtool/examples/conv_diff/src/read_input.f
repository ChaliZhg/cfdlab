      subroutine read_input(mode, nc, xi_1, CFL, niter)
      implicit none
      integer :: mode, nc, niter
      real    :: xi_1, CFL

      print*,'Reading input from param.h ...'
      open(10, file='param.in', status='old')
      read(10,*)mode
      read(10,*)nc 
      read(10,*)xi_1 
      read(10,*)CFL
      read(10,*)niter
      close(10)

      end
