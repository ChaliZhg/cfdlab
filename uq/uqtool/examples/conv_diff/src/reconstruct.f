      program main
      implicit none
      integer :: nc, i
      real,allocatable :: xv(:), xc(:), dx(:)

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

      call reconstruct(nc,xv,xc,dx,'primal.dat')
      call reconstruct(nc,xv,xc,dx,'adjoint.dat')
      call reconstruct(nc,xv,xc,dx,'dprimal.dat')
      call reconstruct(nc,xv,xc,dx,'dadjoint.dat')

      stop
      end

c     Reconstruct solution in filename and write back into same file
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

      open(10, file=filename)
      do i=1,nc
         write(10,'(e24.14,1x,e24.14)') xc(i)-0.25d0*dx(i),
     <   q(i)-0.25*grads(i)*dx(i)
         write(10,'(e24.14,1x,e24.14)') xc(i)+0.25d0*dx(i),
     <   q(i)+0.25*grads(i)*dx(i)
        enddo
      close(10)

      end
