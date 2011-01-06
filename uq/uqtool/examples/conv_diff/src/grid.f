c     Make uniform grid
      program grid
      implicit none
      real, allocatable :: xv(:)
      real    dx
      integer nc
      integer i

      write(*,*) 'Number of cells'
      read(*,*) nc

      allocate( xv(nc+1) )

c     Set up mesh

      dx = 1.0/nc

      write(*,*)'Writing into grid.dat ...'
      open(10, file='grid.dat')
      write(10,*) nc
      do i=1,nc+1
         xv(i)=(i-1)*dx
         write(10,'(e24.14)') xv(i)
      enddo
      close(10)

      deallocate(xv)

      stop
      end
