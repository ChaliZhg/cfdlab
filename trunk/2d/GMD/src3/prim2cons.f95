subroutine prim2cons(rho, vex, vey, vez, pre, con)
   use comvar
   implicit none

   real :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: pre(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: con(nvar, -1:nx+2, -1:ny+2, -1:nz+2)

   integer :: i, j, k

   do i=-1,nx+2
      do j=-1,ny+2
         do k=-1,nz+2

            con(1,i,j,k) = rho(i,j,k)
            con(2,i,j,k) = rho(i,j,k)*vex(i,j,k)
            con(3,i,j,k) = rho(i,j,k)*vey(i,j,k)
            con(4,i,j,k) = rho(i,j,k)*vez(i,j,k)
            con(5,i,j,k) = pre(i,j,k)/(gamma-1.0) + &
                        0.5*rho(i,j,k)*(vex(i,j,k)**2 + vey(i,j,k)**2 + vez(i,j,k)**2)
         enddo
      enddo
   enddo

end subroutine prim2cons
