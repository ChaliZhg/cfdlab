c convert msh mesh to vtx format
c msh file must contain only quads, no lines or points
      implicit none
      integer nv, ncell, i, idum, v1, v2, v3, v4
      real*8 x, y, z


      open(10,file='karman.msh',status='old')
      open(11,file='karman.vtx')

      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*)
      read(10,*) nv
      write(*,*) 'Number of points =', nv
      write(11,*) nv
      do i=1,nv
         read(10,*) idum, x, y, z
         write(11,*) x, y
      enddo

      read(10,*)
      read(10,*)
      read(10,*) ncell
      write(*,*) 'Number of cells =', ncell
      write(11,*) ncell
      do i=1,ncell
         read(10,*) idum,idum,idum,idum,idum,v1,v2,v3,v4
         write(11,*) v1,v2,v3,v4
      enddo

      close(10)
      close(11)

      stop
      end
