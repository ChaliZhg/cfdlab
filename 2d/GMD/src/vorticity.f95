subroutine vorticity(rho, vex, vey, pre, omg)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)
   real    :: omg( 1:nx+1,  1:ny+1)

   integer :: i, j
   real    :: x, y, v_x, u_y

   ! vorticity is computed at vertices (not at cell centers)
   do i=1,nx+1
      do j=1,ny+1
         x = xmin + (i-1)*dx - 0.5*dx
         y = ymin + (j-1)*dy - 0.5*dy

         v_x = 0.5*(vey(i,j-1) + vey(i,j)) - 0.5*(vey(i-1,j-1) + vey(i-1,j))
         v_x = v_x/dx

         u_y = 0.5*(vex(i-1,j) + vex(i,j)) - 0.5*(vex(i-1,j-1) + vex(i,j-1))
         u_y = u_y/dy

         omg(i,j) = v_x - u_y

      enddo
   enddo


end subroutine vorticity
