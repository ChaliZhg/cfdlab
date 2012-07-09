subroutine init_cond_multvortdiff(time, rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: time
   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, r2, Temp, circ, circ1, circ2
   real    :: vx0, vy0

   xmin = 0.0
   xmax = 2.0*M_PI
   ymin = 0.0
   ymax = 2.0*M_PI

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         rho(i,j) =  1.0
         vex(i,j) = -cos(x) * sin(y) * exp(-2.0 * mu * time)
         vey(i,j) =  sin(x) * cos(y) * exp(-2.0 * mu * time)
         pre(i,j) =  50.0 - 0.25 * (cos(2.0*x) + cos(2.0*y)) * exp(-4.0 * mu * time)
      enddo
   enddo

end subroutine init_cond_multvortdiff
