subroutine timestep(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: speed, eig

   speed = 0.0

   do i=1,nx
      do j=1,ny
         eig = sqrt(vex(i,j)**2 + vey(i,j)**2) + &
               sqrt(gamma*pre(i,j)/rho(i,j))
         speed = max(speed, eig)
      enddo
   enddo

   dt = cfl*dx/speed

   write(*,*)'Time step =', dt

end subroutine timestep
