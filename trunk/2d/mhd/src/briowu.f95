! http://www.astro.princeton.edu/~jstone/Athena/tests/orszag-tang/pagesource.html
subroutine briowu(pri)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, alpha, Bpar, Bperp, vperp

   gamma = 5.0/3.0
   final_time = 0.2
   xperiod    = no
   yperiod    = yes

   xmin =-1.0
   xmax = 1.0
   dx = (xmax - xmin)/nx
   dy = dx

   ymin = 0.0
   ymax = ny*dy

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         if(x.lt.0.0)then
            pri(1,i,j) = 1.0

            pri(2,i,j) = 0.0
            pri(3,i,j) = 0.0
            pri(4,i,j) = 0.0

            pri(5,i,j) = 1.0

            pri(6,i,j) = 0.75
            pri(7,i,j) = 1.0
            pri(8,i,j) = 0.0
         else
            pri(1,i,j) = 0.125

            pri(2,i,j) = 0.0
            pri(3,i,j) = 0.0
            pri(4,i,j) = 0.0

            pri(5,i,j) = 0.1

            pri(6,i,j) = 0.75
            pri(7,i,j) =-1.0
            pri(8,i,j) = 0.0
         endif
      enddo
   enddo

end subroutine briowu
