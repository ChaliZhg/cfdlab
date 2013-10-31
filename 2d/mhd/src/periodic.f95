subroutine periodic(con)
   use comvar
   use omp_lib
   implicit none

   real :: con(nvar,-1:nx+2,-1:ny+2)

   integer :: i, j

   if(xperiod == yes)then
      !$omp parallel do
      do j=1,ny
         con(:,-1,  j) = con(:,nx-1,j)
         con(:, 0,  j) = con(:,nx,  j)
         con(:,nx+1,j) = con(:,1,   j)
         con(:,nx+2,j) = con(:,2,   j)
      enddo
      !$omp end parallel do
   endif

   if(yperiod == yes)then
      !$omp parallel do
      do i=1,nx
         con(:,i,  -1) = con(:,i,ny-1)
         con(:,i,   0) = con(:,i,ny  )
         con(:,i,ny+1) = con(:,i,1   )
         con(:,i,ny+2) = con(:,i,2   )
      enddo
      !$omp end parallel do
   endif

   if(xperiod == yes .and. yperiod == yes)then

      ! lower left corner
      con(:, 0, 0) = con(:,nx,  ny)
      con(:,-1, 0) = con(:,nx-1,ny)
      con(:, 0,-1) = con(:,nx,  ny-1)
      con(:,-1,-1) = con(:,nx-1,ny-1)

      ! lower right corner
      con(:,nx+1, 0) = con(:,1,ny)
      con(:,nx+2, 0) = con(:,2,ny)
      con(:,nx+1,-1) = con(:,1,ny-1)
      con(:,nx+2,-1) = con(:,2,ny-1)

      ! upper left corner
      con(:, 0,ny+1) = con(:,nx,  1)
      con(:, 0,ny+2) = con(:,nx,  2)
      con(:,-1,ny+1) = con(:,nx-1,1)
      con(:,-1,ny+2) = con(:,nx-1,2)

      ! upper right corner
      con(:,nx+1,ny+1) = con(:,1,1)
      con(:,nx+2,ny+1) = con(:,2,1)
      con(:,nx+1,ny+2) = con(:,1,2)
      con(:,nx+2,ny+2) = con(:,2,2)

   endif

end subroutine periodic
