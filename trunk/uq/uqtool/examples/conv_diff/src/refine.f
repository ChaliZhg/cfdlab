      program main
      implicit none
      real,allocatable :: error(:)
      integer :: nc, nc_max, nc_ref
      real,allocatable :: xv(:)
      real :: refine_fraction
      integer :: i
      integer,allocatable :: idx(:)
      logical,allocatable :: to_refine(:)

      nc_max = 100
      refine_fraction = 0.05

c     Read current grid
      print*,'Reading current grid'
      open(10,file='grid.dat',status='old')
      read(10,*) nc
      allocate(xv(nc+1))
      do i=1,nc+1
         read(10,*) xv(i)
      enddo
      close(10)

      allocate(error(nc))

c     Read error indicator
      open(10,file='error0.dat',status='old')
      do i=1,nc
         read(10,*) error(i)
      enddo
      close(10)

      error = abs(error)

c     Dont refine if upper limit of nc is reached
      if(nc.ge.nc_max) stop

c     Number of cells to refine
      nc_ref = refine_fraction * nc

c     Find nc_ref cells with largest error
      allocate(idx(nc))
      do i=1,nc
         idx(i) = i
      enddo
      call ssort(error, idx, nc, 0)

c     Flag nc_ref cells for refinement
      allocate(to_refine(nc))
      to_refine = .FALSE.
      do i=1,nc_ref
         to_refine(idx(i)) = .TRUE.
      enddo

c     Number of new cells
      print*,'Writing new grid, nc =',nc+nc_ref
      open(10,file='grid.dat')
      write(10,*) nc + nc_ref
      do i=1,nc
         write(10,'(e24.14)') xv(i)
         if(to_refine(i))then
            print*,'Refining ',i,' error =',error(idx(i))
            write(10,'(e24.14)') 0.5*(xv(i) + xv(i+1))
         endif
      enddo
      write(10,'(e24.14)') xv(nc+1)
      close(10)

      open(10,file='cell_map.dat')
      write(10,*) nc+nc_ref, nc+nc_ref
      do i=1,nc+nc_ref
         write(10,*) i
      enddo
      close(10)

      stop
      end
