c     Make uniform grid
      program grid
      implicit none
      real, allocatable :: xc(:)
      real    dx
      integer nc
      integer i

      write(*,*) 'Number of cells'
      read(*,*) nc

      allocate( xc(nc) )

c     Set up mesh

      dx = 1.0/nc

      write(*,*)'Writing into grid.dat ...'
      open(10, file='grid.dat')
      write(10,*) nc
      do i=1,nc
         xc(i)=(float(i)-0.5d0)*dx
         write(10,'(e24.14)') xc(i)
      enddo
      close(10)

      deallocate(xc)

      stop
      end
