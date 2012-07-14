subroutine cons2prim(con, rho, vex, vey, vez, pre)
   use comvar
   implicit none

   real :: con(nvar, -1:nx+2, -1:ny+2, -1:nz+2)
   real :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: pre(-1:nx+2, -1:ny+2, -1:nz+2)

   integer :: i, j, k

   do i=-1,nx+2
      do j=-1,ny+2
         do k=-1,nz+2
         rho(i,j,k) = con(1,i,j,k)
         vex(i,j,k) = con(2,i,j,k)/con(1,i,j,k)
         vey(i,j,k) = con(3,i,j,k)/con(1,i,j,k)
         vez(i,j,k) = con(4,i,j,k)/con(1,i,j,k)
         pre(i,j,k) = (gamma-1.0)*(con(5,i,j,k) - &
                    0.5*(con(2,i,j,k)**2 + con(3,i,j,k)**2 + con(4,i,j,k)**2)/con(1,i,j,k))
         enddo
      enddo
   enddo

end subroutine cons2prim
