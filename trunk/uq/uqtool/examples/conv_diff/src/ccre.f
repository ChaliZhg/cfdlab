      program main
      implicit none
      integer :: nc, i
      real,allocatable :: v(:), R(:), AR(:), dq(:), dv(:)

      open(10,file='grid.dat',status='old')
      read(10,*) nc
      close(10)

      allocate( v(nc) )
      allocate( R(nc) )
      allocate(AR(nc) )
      allocate(dq(nc) )
      allocate(dv(nc) )

c     Read adjoint
c     Read primal residual
c     Read adjoint residual
c     Read dprimal
c     Read dadjoint
      open(10,file='adjoint.dat',   status='old')
      open(11,file='p_residual.dat',status='old')
      open(12,file='a_residual.dat',status='old')
      open(13,file='dprimal.dat',   status='old')
      open(14,file='dadjoint.dat',  status='old')
      do i=1,nc
         read(10,*) v(i)
         read(11,*) R(i)
         read(12,*) AR(i)
         read(13,*) dq(i)
         read(14,*) dv(i)
      enddo

      close(10)
      close(11)
      close(12)
      close(13)
      close(14)

      open(10,file='VdotR.dat')
      open(11,file='RE.dat')
      do i=1,nc
         write(10,'(e24.14)') v(i)*R(i)
         write(11,'(e24.14)') dv(i)*R(i) + dq(i)*AR(i)
      enddo
      close(10)
      close(11)

      stop
      end
