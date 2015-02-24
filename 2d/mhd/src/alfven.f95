! Same as in Toth
subroutine alfven(pri)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, alpha, Bpar, Bperp, vperp

   gamma = 5.0/3.0
   alpha = 30.0 * PI/180.0
   Bpar  = 1.0
   final_time = 5.0
   xperiod    = yes
   yperiod    = yes

   xmin = 0.0
   xmax = 1.0/cos(alpha)
   ymin = 0.0
   ymax = 1.0/sin(alpha)

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         vperp = 0.1 * sin(2.0*PI*(x*cos(alpha) + y*sin(alpha)))
         Bperp = vperp

         pri(1,i,j) = 1.0

         pri(2,i,j) =-vperp * sin(alpha)
         pri(3,i,j) = vperp * cos(alpha)
         pri(4,i,j) = 0.1 * cos(2.0*PI*(x*cos(alpha)+y*sin(alpha)))

         pri(5,i,j) = 0.1

         pri(6,i,j) = Bpar * cos(alpha) - Bperp * sin(alpha)
         pri(7,i,j) = Bpar * sin(alpha) + Bperp * cos(alpha)
         pri(8,i,j) = pri(4,i,j)
      enddo
   enddo

end subroutine alfven
