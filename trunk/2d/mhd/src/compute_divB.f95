subroutine compute_divB(pri, divB, maxdivB)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)
   real    :: divB( 1:nx+1,  1:ny+1)
   real    :: maxdivB

   integer :: i, j
   real    :: Bx_x, By_y

   ! vorticity is computed at vertices (not at cell centers)
   maxdivB = 0.0
   do i=1,nx+1
      do j=1,ny+1
         Bx_x = 0.5*(pri(6,i,j-1) + pri(6,i,j)) - 0.5*(pri(6,i-1,j-1) + pri(6,i-1,j))
         Bx_x = Bx_x/dx

         By_y = 0.5*(pri(7,i-1,j) + pri(7,i,j)) - 0.5*(pri(7,i-1,j-1) + pri(7,i,j-1))
         By_y = By_y/dy

         divB(i,j) = Bx_x + By_y
         maxdivB = max(maxdivB, divB(i,j))
      enddo
   enddo


end subroutine compute_divB
