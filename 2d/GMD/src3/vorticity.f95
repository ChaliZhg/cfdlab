subroutine vorticity(rho, vex, vey, vez, pre, omg)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: pre(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: omg( 1:nx+1,  1:ny+1,  1:nz+1)

   integer :: i, j, k
   real    :: v_x, u_y

   ! vorticity is computed at vertices (not at cell centers)
   do i=1,nx+1
      do j=1,ny+1
         do k=1,nz+1
            v_x = 0.25*(vey(i,j-1,k) + vey(i,j,k)) - &
                  0.25*(vey(i-1,j-1,k) + vey(i-1,j,k))
            v_x = v_x/dx

            u_y = 0.25*(vex(i-1,j,k) + vex(i,j,k)) - &
                  0.25*(vex(i-1,j-1,k) + vex(i,j-1,k))
            u_y = u_y/dy

            omg(i,j,k) = v_x - u_y
         enddo
      enddo
   enddo


end subroutine vorticity
