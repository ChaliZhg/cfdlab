subroutine savevort(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, v_x, u_y
   character(len=512) :: filename

   filename = 'omg'
   call getfilename(filename, fileid)

   open(10,file=trim(filename))
   write(10,*)'TITLE = "vortex flow"'
   write(10,*)'VARIABLES = "x", "y", "Vorticity"'
   write(10,*)'ZONE I=',nx+1,', J=',ny+1,', DATAPACKING=POINT'

   ! vorticity is computed at vertices (not at cell centers)
   do i=1,nx+1
      do j=1,ny+1
         x = xmin + (i-1)*dx - 0.5*dx
         y = ymin + (j-1)*dy - 0.5*dy

         v_x = 0.5*(vey(i,j-1) + vey(i,j)) - 0.5*(vey(i-1,j-1) + vey(i-1,j))
         v_x = v_x/dx

         u_y = 0.5*(vex(i-1,j) + vex(i,j)) - 0.5*(vex(i-1,j-1) + vex(i,j-1))
         u_y = u_y/dy

         write(10,'(3E24.14)') x, y, v_x - u_y

      enddo
   enddo

   close(10)

end subroutine savevort
