subroutine global_quantities(pri, ke, entropy)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)
   real    :: ke, entropy

   integer :: i, j

   ke  = 0.0
   entropy = 0.0

   do i=1,nx
      do j=1,ny
         ke = ke + pri(1,i,j)*(pri(2,i,j)**2 + pri(3,i,j)**2 + pri(4,i,j)**2)
         entropy = entropy + pri(1,i,j) * (log(pri(5,i,j) / pri(1,i,j)**gamma))
      enddo
   enddo

   ke = 0.5 * ke * dx * dy
   entropy = -entropy * dx * dy / (gamma - 1.0)

end subroutine global_quantities
