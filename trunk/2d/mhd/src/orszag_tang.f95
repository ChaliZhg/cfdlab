! http://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
subroutine orszag_tang(pri)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y

   gamma      = 5.0/3.0
   final_time = 0.5
   xperiod    = yes
   yperiod    = yes

   xmin = 0.0
   xmax = 1.0
   ymin = 0.0
   ymax = 1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy
         pri(1,i,j) = 25.0/(36.0*PI)
         pri(2,i,j) =-sin(2.0*PI*y)
         pri(3,i,j) = sin(2.0*PI*x)
         pri(4,i,j) = 0.0
         pri(5,i,j) = 5.0/(12.0*PI)
         pri(6,i,j) =-sin(2.0*PI*y)/sqrt(4.0*PI)
         pri(7,i,j) = sin(4.0*PI*x)/sqrt(4.0*PI)
         pri(8,i,j) = 0.0
      enddo
   enddo

end subroutine orszag_tang
