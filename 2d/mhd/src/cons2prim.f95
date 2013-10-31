subroutine cons2prim(con, pri)
   use comvar
   use omp_lib
   implicit none

   real :: con(nvar, -1:nx+2, -1:ny+2)
   real :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j

   !$omp parallel do
   do i=-1,nx+2
      do j=-1,ny+2
         ! density
         pri(1,i,j) = con(1,i,j)
         ! velocity
         pri(2:4,i,j) = con(2:4,i,j)/con(1,i,j)
         ! pressure
         pri(5,i,j) = (gamma-1.0)*(con(5,i,j) &
                      - 0.5*(con(2,i,j)**2 + con(3,i,j)**2 + con(4,i,j)**2)/con(1,i,j) &
                      - 0.5*(con(6,i,j)**2 + con(7,i,j)**2 + con(8,i,j)**2))
         ! magnetic field
         pri(6:8,i,j) = con(6:8,i,j)
      enddo
   enddo
   !$omp end parallel do

end subroutine cons2prim
