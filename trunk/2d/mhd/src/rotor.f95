! First rotor problem from Toth
subroutine rotor(pri)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, r, r0, r1, f, u0

   gamma = 1.4
   final_time = 0.15
   xperiod    = yes
   yperiod    = yes

   xmin = 0.0
   xmax = 1.0
   ymin = 0.0
   ymax = 1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   r0 = 0.1
   r1 = 0.115
   u0 = 2.0

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         r = sqrt((x-0.5)**2 + (y-0.5)**2)

         if(r.le.r0)then
            pri(1,i,j) = 10.0
            pri(2,i,j) =-u0 * (y-0.5)/r0
            pri(3,i,j) = u0 * (x-0.5)/r0
         elseif(r.ge.r1)then
            pri(1,i,j) = 1.0
            pri(2,i,j) = 0.0
            pri(3,i,j) = 0.0
         else
            f = (r1 - r)/(r1 - r0)
            pri(1,i,j) = 1.0 + 9.0 * f
            pri(2,i,j) =-f * u0 * (y-0.5)/r
            pri(3,i,j) = f * u0 * (x-0.5)/r
         endif

         pri(4,i,j) = 0.0
         pri(5,i,j) = 1.0
         pri(6,i,j) = 5.0/sqrt(4.0*PI)
         pri(7,i,j) = 0.0
         pri(8,i,j) = 0.0
      enddo
   enddo

end subroutine rotor
