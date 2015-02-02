!! http://www.astro.princeton.edu/~jstone/Athena/tests/kh/kh.html
subroutine kelvin_helmholtz(pri)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, r

   gamma      = 1.4
   final_time = 5.0
   xperiod    = yes
   yperiod    = yes

   xmin =-0.5
   xmax = 0.5
   ymin =-0.5
   ymax = 0.5

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         call random_number(r)
         if(abs(y) < 0.25)then
            pri(1,i,j) = 2.0
            pri(2,i,j) = 0.5 + 0.01 * (2.0*r - 1.0)
         else
            pri(1,i,j) = 1.0
            pri(2,i,j) =-0.5 + 0.01 * (2.0*r - 1.0)
         endif
         call random_number(r)
         pri(3,i,j) = 0.01 * (2.0*r - 1.0)
         pri(4,i,j) = 0.0
         pri(5,i,j) = 2.5
         pri(6,i,j) = 0.5
         pri(7,i,j) = 0.0
         pri(8,i,j) = 0.0
      enddo
   enddo

end subroutine kelvin_helmholtz
