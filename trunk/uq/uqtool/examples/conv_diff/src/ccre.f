      program main
      implicit none
      integer :: nc, nc_h, i, j
      real,allocatable :: v(:), R(:), AR(:), dq(:), dv(:)
      real,allocatable :: cc(:), re(:)
      integer,allocatable :: cell_map(:)

      write(*,*) 'Reading cell_map.dat'
      open(10,file='cell_map.dat',status='old')
      read(10,*) nc, nc_h
      allocate(cell_map(nc_h) )
      do i=1,nc_h
         read(10,*) cell_map(i)
      enddo
      close(10)

      allocate( v(nc_h) )
      allocate( R(nc_h) )
      allocate(AR(nc_h) )
      allocate(dq(nc_h) )
      allocate(dv(nc_h) )

      allocate(cc(nc))
      allocate(re(nc))

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
      do i=1,nc_h
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

      cc = 0.0
      re = 0.0
      do i=1,nc_h
         j = cell_map(i)
         cc(j) = cc(j) + v(i)*R(i)
         re(j) = re(j) + dv(i)*R(i) + dq(i)*AR(i)
      enddo

      open(10,file='VdotR.dat')
      open(11,file='RE.dat')
      do i=1,nc
         write(10,'(e24.14)') cc(i)
         write(11,'(e24.14)') re(i)
      enddo
      close(10)
      close(11)

      stop
      end
