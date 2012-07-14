subroutine timestep(rho, vex, vey, vez, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: pre(-1:nx+2, -1:ny+2, -1:nz+2)

   integer :: i, j, k
   real    :: speed, eig, ds

   speed = 0.0

   do i=1,nx
      do j=1,ny
         do k=1,nz
            eig = sqrt(vex(i,j,k)**2 + vey(i,j,k)**2 + vez(i,j,k)**2) + &
                  sqrt(gamma*pre(i,j,k)/rho(i,j,k))
            speed = max(speed, eig)
         enddo
      enddo
   enddo

   ds = min(dx, min(dy,dz))
   dt = ds/speed
   if(mu > 0.0)then
      dt = min(dt, ds**2/mu)
   endif

   dt = cfl*dt

   write(*,*)'Time step =', dt

end subroutine timestep
