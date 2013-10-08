subroutine init_cond(pri, co1)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)
   real    :: co1(nvar, -1:nx+2, -1:ny+2)


   integer :: i, j
   real    :: x, y, q2, B2

   final_time = 0.5

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
         ! Orszag-Tang vortex
         pri(1,i,j) = 1.0
         pri(2,i,j) = sin(PI*y)
         pri(3,i,j) = sin(PI*x)
         pri(4,i,j) = 0.0
         pri(5,i,j) = 0.6
         pri(6,i,j) = -sin(PI*y)
         pri(7,i,j) = sin(2.0*PI*x)
         pri(8,i,j) = 0.0
         ! Convert to conserved variables
         co1(1,i,j)   = pri(1,i,j)
         co1(2:4,i,j) = pri(2:4,i,j) * pri(1,i,j)
         q2 = pri(2,i,j)**2 + pri(3,i,j)**2 + pri(4,i,j)**2
         B2 = pri(6,i,j)**2 + pri(7,i,j)**2 + pri(8,i,j)**2
         co1(5,i,j)   = pri(5,i,j)/(gamma-1.0) + 0.5*pri(1,i,j)*q2 + 0.5*B2
         co1(6:8,i,j) = pri(6:8,i,j)
      enddo
   enddo

end subroutine init_cond
