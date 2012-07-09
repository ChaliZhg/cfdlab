subroutine init_cond_isen(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, r2, Temp, circ, circ1, circ2
   real    :: vx0, vy0

   xmin =-5.0
   xmax = 5.0
   ymin =-5.0
   ymax = 5.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   circ = 5.0
   circ1= (gamma-1.0)*circ**2/(8.0*gamma*M_PI**2)
   circ2= circ/(2.0*M_PI)
   ! Translation velocity of vortex
   vx0  = 0.0
   vy0  = 0.0

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         r2   = x**2 + y**2
         Temp = 1.0 - circ1*exp(1.0-r2)

         rho(i,j) =  Temp**(1.0/(gamma-1.0))
         vex(i,j) = -circ2*y*exp(0.5*(1.0-r2)) + vx0
         vey(i,j) =  circ2*x*exp(0.5*(1.0-r2)) + vy0
         pre(i,j) =  rho(i,j)*Temp
      enddo
   enddo

end subroutine init_cond_isen
