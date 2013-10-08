subroutine saveprim(t, pri)
   use comvar
   implicit none

   real    :: t
   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, q, a, m
   character(len=512) :: filename

   filename = 'sol'
   call getfilename(filename, fileid_sol)

   open(10,file=trim(filename))
   write(10,*)'TITLE = "mhd"'
   write(10,*)'VARIABLES = "x", "y", "Density", "Vx", "Vy", "Vz", "Pressure", ', &
              '"Bx", "By", "Bz", "Mach"'
   write(10,*)'ZONE STRANDID=1, SOLUTIONTIME=',t,', I=',nx,', J=',ny,&
              ', DATAPACKING=POINT'

   do j=1,ny
      do i=1,nx
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         q = sqrt(pri(2,i,j)**2 + pri(3,i,j)**2 + pri(4,i,j)**2)
         a = sqrt(gamma*pri(5,i,j)/pri(1,i,j))
         m = q/a
         write(10,'(11E24.14)') x, y, pri(1,i,j), pri(2,i,j), pri(3,i,j), &
                               pri(4,i,j), pri(5,i,j), pri(6,i,j), &
                               pri(7,i,j), pri(8,i,j), m

      enddo
   enddo

   close(10)

end subroutine saveprim
