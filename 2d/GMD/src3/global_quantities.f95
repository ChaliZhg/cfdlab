subroutine global_quantities(rho, vex, vey, vez, pre, ke, entropy)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: pre(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: omg( 1:nx+1,  1:ny+1, -1:nz+2)
   real    :: ke, entropy

   integer :: i, j, k

   ke = 0.0
   entropy = 0.0

   ! vorticity is computed at vertices (not at cell centers)
   do i=1,nx
      do j=1,ny
         do k=1,nz
            ke = ke + rho(i,j,k)*(vex(i,j,k)**2 + vey(i,j,k)**2 + vez(i,j,k)**2)
            entropy = entropy + rho(i,j,k) * (log(pre(i,j,k) / rho(i,j,k)**gamma))
         enddo
      enddo
   enddo

   ke = 0.5 * ke * dx * dy * dz
   entropy = entropy * dx * dy * dz

end subroutine global_quantities
