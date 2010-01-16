subroutine saveprim(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y
   character(len=512) :: filename

   filename = 'sol'
   call getfilename(filename, fileid_sol)

   open(10,file=trim(filename))
   write(10,*)'TITLE = "vortex flow"'
   write(10,*)'VARIABLES = "x", "y", "Density", "Velx", "Vely", "Pressure"'
   write(10,*)'ZONE I=',nx,', J=',ny,', DATAPACKING=POINT'

   do i=1,nx
      do j=1,ny
         x = xmin + (i-1)*dx
         y = ymin + (j-1)*dy

         write(10,'(6E24.14)') x, y, rho(i,j), vex(i,j), vey(i,j), pre(i,j)

      enddo
   enddo

   close(10)

end subroutine saveprim
