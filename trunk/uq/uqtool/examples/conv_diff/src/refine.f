      program main
      implicit none
      integer, parameter :: nc_max = 1000
      real :: xv_h(nc_max)
      real :: xv_hh(nc_max)
      real,allocatable :: error(:)
      integer :: nc, nc_ref, nc_h, nc_hh
      real,allocatable :: xv(:)
      real :: refine_fraction
      integer :: i, c
      logical,allocatable :: to_refine(:)
      real :: dx1, dx2, dx_min

c     Fraction of cells to refine
      refine_fraction = 0.05
c     Lower bound on cell size
      dx_min = 1.0/50.0

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

c     Dont refine if upper limit of nc is reached
      if(nc.ge.nc_max) stop

c     Flag nc_ref cells for refinement
      allocate(to_refine(nc))
      to_refine = .FALSE.
      nc_ref = 0
      do i=1,nc
         if(xv(i+1)-xv(i).gt.dx_min .and. error(i).gt.0)then
            to_refine(i) = .TRUE.
            nc_ref = nc_ref + 1
         endif
      enddo

c     New number of cells
      nc_h = nc + nc_ref

c     Refine grid
      c = 1
      do i=1,nc
         xv_h(c) = xv(i)
         c = c + 1

         if(to_refine(i))then
            print*,'Refining ',i
            xv_h(c) = 0.5*(xv(i) + xv(i+1))
            c = c + 1
         endif
      enddo
      xv_h(c) = xv(nc+1)
      if(c-1.ne.nc_h) stop "c not equal to nc_h"

      deallocate(to_refine)
      allocate(to_refine(nc_h))
      to_refine = .false.
      do i=1,nc_h-1
         dx1 = xv_h(i+1) - xv_h(i)
         dx2 = xv_h(i+2) - xv_h(i+1)
         if(dx1.gt.2.0*dx2) to_refine(i)   = .true.
         if(dx2.gt.2.0*dx1) to_refine(i+1) = .true.
      enddo

c     Refine grid
      c = 1
      do i=1,nc_h
         xv_hh(c) = xv_h(i)
         c = c + 1

         if(to_refine(i))then
            xv_hh(c) = 0.5*(xv_h(i) + xv_h(i+1))
            c = c + 1
         endif
      enddo
      xv_hh(c) = xv_h(nc_h+1)
      nc_hh = c - 1

c     Number of new cells
      print*,'Writing new grid, nc =',nc_hh
      open(10,file='grid.dat')
      write(10,*) nc_hh
      do i=1,nc_hh+1
         write(10,'(e24.14)') xv_hh(i)
      enddo
      close(10)

      open(10,file='cell_map.dat')
      write(10,*) nc_hh, nc_hh
      do i=1,nc_hh
         write(10,*) i
      enddo
      close(10)

      stop
      end
