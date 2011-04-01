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
      real    :: q1, q2

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
         if(i.eq.1)then
            call do_reconstruct(q(i),   xv(i),   xv(i+1),
     1                          q(i+1), xv(i+1), xv(i+2),
     2                          q(i+2), xv(i+2), xv(i+3), q1, q2)
         else if(i.eq.nc)then
            call do_reconstruct(q(i),   xv(i),   xv(i+1),
     1                          q(i-1), xv(i-1), xv(i),
     2                          q(i-2), xv(i-2), xv(i-1), q1, q2)
         else
            call do_reconstruct(q(i),   xv(i),   xv(i+1),
     1                          q(i-1), xv(i-1), xv(i),
     2                          q(i+1), xv(i+1), xv(i+2), q1, q2)
         endif
c        write(10,'(e24.14)') q(i)-0.25*grads(i)*dx(i)
c        write(10,'(e24.14)') q(i)+0.25*grads(i)*dx(i)
         write(10,'(e24.14)') q1
         write(10,'(e24.14)') q2
      enddo
      close(10)

      end
c
c Fit 
c     q(x) = q0 + a (x - x0) + b[ (x - x0)^2 - h0^2/12]
c so that cell averages are preserved
c
      subroutine do_reconstruct(q0, x01, x02, 
     1                          q1, x11, x12, 
     2                          q2, x21, x22,
     3                          qr1, qr2)
      implicit none

      real :: q0, x01, x02, q1, x11, x12, q2, x21, x22, qr1, qr2
      real :: h0, h1, h2
      real :: mat(2,2), rhs(2)
      real :: x, x0, fun1, fun2
      real :: a, b, det

      fun1(x) = (x - x0)**2 / 2.0
      fun2(x) = (x - x0)**3 / 3.0 - h0**2 * x / 12.0

      x0 = 0.5 * (x01 + x02)
      h0 = x02 - x01
      h1 = x12 - x11
      h2 = x22 - x21

      mat(1,1) = fun1(x12) - fun1(x11)
      mat(1,2) = fun2(x12) - fun2(x11)

      mat(2,1) = fun1(x22) - fun1(x21)
      mat(2,2) = fun2(x22) - fun2(x21)

      rhs(1) = (q1 - q0) * h1
      rhs(2) = (q2 - q0) * h2

      det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
      a   = ( mat(2,2) * rhs(1) - mat(1,2) * rhs(2) ) / det
      b   = (-mat(2,1) * rhs(1) + mat(1,1) * rhs(2) ) / det

c     Average over two cells
      qr1 = q0 * h0/2.0 + a * (fun1(x0) - fun1(x01)) +
     1                    b * (fun2(x0) - fun2(x01))
      qr1 = qr1 * 2.0 / h0

      qr2 = q0 * h0/2.0 + a * (fun1(x02) - fun1(x0)) +
     1                    b * (fun2(x02) - fun2(x0))
      qr2 = qr2 * 2.0 / h0

      end
