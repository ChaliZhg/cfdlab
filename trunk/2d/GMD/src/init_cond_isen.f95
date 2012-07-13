subroutine init_cond_isen(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, r2, Temp, circ, circ1, circ2, mach_inf
   real    :: x0, y0, theta

   final_time = 50.0

   xmin =-5.0
   xmax = 5.0
   ymin =-5.0
   ymax = 5.0

   xmin =-10.0
   xmax = 10.0
   ymin =-10.0
   ymax = 10.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   x0 = 0.0
   y0 = 0.0
   mach_inf = 0.5
   theta = 0.0
   circ = 5.0

   theta= theta*M_PI/180.0
   circ1= (gamma-1.0)*circ**2*mach_inf**2/(8.0*M_PI**2)
   circ2= circ/(2.0*M_PI)

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         r2   = (x-x0)**2 + (y-y0)**2

         rho(i,j) = (1.0 - circ1*exp(1.0-r2))**(1.0/(gamma-1.0))
         vex(i,j) = mach_inf * (cos(theta) - circ2*(y-y0)*exp(0.5*(1.0-r2)))
         vey(i,j) = mach_inf * (sin(theta) + circ2*(x-x0)*exp(0.5*(1.0-r2)))
         pre(i,j) = rho(i,j)**gamma / gamma
      enddo
   enddo

end subroutine init_cond_isen
