subroutine prim2cons(rho, vex, vey, pre, con)
   use comvar
   implicit none

   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: pre(-1:nx+2, -1:ny+2)
   real :: con(4, -1:nx+2, -1:ny+2)

   integer :: i, j

   do i=-1,nx+2
      do j=-1,ny+2

         con(1,i,j) = rho(i,j)
         con(2,i,j) = rho(i,j)*vex(i,j)
         con(3,i,j) = rho(i,j)*vey(i,j)
         con(4,i,j) = pre(i,j)/(gamma-1.0) + &
                      0.5*rho(i,j)*(vex(i,j)**2 + vey(i,j)**2)
      enddo
   enddo

end subroutine prim2cons
