      program main
      implicit none
      integer :: nc, nc_h, i, c
      real,allocatable :: xv(:), xc(:), dx(:)
      real,allocatable :: xv_h(:)
      integer,allocatable :: cell_map(:)

c     Read current grid
      open(10,file='grid.dat',status='old')
      read(10,*) nc
      allocate(xv(nc+1),xc(nc),dx(nc))
      do i=1,nc+1
         read(10,*) xv(i)
      enddo
      close(10)
      do i=1,nc
         xc(i) = 0.5*(xv(i) + xv(i+1))
         dx(i) = xv(i+1) - xv(i)
      enddo

c     Refine grid
      nc_h = 2 * nc
      allocate(xv_h(nc_h+1),cell_map(nc_h))
      c = 1
      do i=1,nc
         xv_h(c) = xv(i)
         cell_map(c) = i
         c = c + 1

         xv_h(c) = 0.5 * (xv(i) + xv(i+1))
         cell_map(c) = i
         c = c + 1
      enddo
      xv_h(nc_h+1) = xv(nc)

      call reconstruct(nc,xv,xc,dx,'primal.dat')
      call reconstruct(nc,xv,xc,dx,'adjoint.dat')
      call reconstruct(nc,xv,xc,dx,'dprimal.dat')
      call reconstruct(nc,xv,xc,dx,'dadjoint.dat')

c     Overwrite grid with refined grid
      write(*,*) 'Writing refined grid into grid.dat ...'
      open(10,file='grid.dat')
      write(10,*) nc_h
      do i=1,nc_h+1
         write(10,'(e24.14)') xv_h(i)
      enddo
      close(10)

c     overwrite cell_map with new map
      write(*,*) 'Writing new cell_map'
      open(10,file='cell_map.dat')
      write(10,*) nc, nc_h
      do i=1,nc_h
         write(10,*) cell_map(i)
      enddo
      close(10)

      stop
      end
c    
c     Reconstruct solution in filename and write back into same file
c    
      subroutine reconstruct(nc,xv,xc,dx,filename)
      implicit none
      integer :: nc
      real    :: xv(*), xc(*), dx(*)
      character :: filename*(*)

      integer :: i
      real    :: q(nc), grads(nc)

      print*,'Reconstructing ', filename

      open(10,file=filename,status='old')
      do i=1,nc
         read(10,*) q(i)
      enddo
      close(10)

      call calc_grad(nc,xc,xv,q,grads)

c     Overwrite filename with reconstructed values on fine mesh
      open(10, file=filename)
      do i=1,nc
         write(10,'(e24.14)') q(i)-0.25*grads(i)*dx(i)
         write(10,'(e24.14)') q(i)+0.25*grads(i)*dx(i)
        enddo
      close(10)

      end
