subroutine global_quantities(rho, vex, vey, pre, ke, entropy)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)
   real    :: omg( 1:nx+1,  1:ny+1)
   real    :: ke, entropy

   integer :: i, j
   real    :: v_x, u_y

   ke = 0.0
   entropy = 0.0

   ! vorticity is computed at vertices (not at cell centers)
   do i=1,nx
      do j=1,ny
         ke = ke + rho(i,j)*(vex(i,j)**2 + vey(i,j)**2)
         entropy = entropy + rho(i,j) * (log(pre(i,j) / rho(i,j)**gamma))
      enddo
   enddo

   ke = 0.5 * ke * dx * dy
   entropy = entropy * dx * dy

end subroutine global_quantities
