!------------------------------------------------------------------------------
! x derivative
!------------------------------------------------------------------------------
real function vertex_gradient_x(i, j, k, u)
   use comvar
   implicit none

   real    :: u(-1:nx+2, -1:ny+2, -1:nz+2)
   integer :: i, j, k

   vertex_gradient_x = u(i+1, j,   k  ) - u(i, j,   k  ) + &
                       u(i+1, j+1, k  ) - u(i, j+1, k  ) + &
                       u(i+1, j,   k+1) - u(i, j,   k+1) + &
                       u(i+1, j+1, k+1) - u(i, j+1, k+1)
   vertex_gradient_x = vertex_gradient_x / 4.0
end function vertex_gradient_x

!------------------------------------------------------------------------------
! y derivative
!------------------------------------------------------------------------------
real function vertex_gradient_y(i, j, k, u)
   use comvar
   implicit none

   real    :: u(-1:nx+2, -1:ny+2, -1:nz+2)
   integer :: i, j, k

   vertex_gradient_y = u(i,   j+1, k  ) - u(i,   j, k  ) + &
                       u(i+1, j+1, k  ) - u(i+1, j, k  ) + &
                       u(i,   j+1, k+1) - u(i,   j, k+1) + &
                       u(i+1, j+1, k+1) - u(i+1, j, k+1)
   vertex_gradient_y = vertex_gradient_y / 4.0
end function vertex_gradient_y

!------------------------------------------------------------------------------
! z derivative
!------------------------------------------------------------------------------
real function vertex_gradient_z(i, j, k, u)
   use comvar
   implicit none

   real    :: u(-1:nx+2, -1:ny+2, -1:nz+2)
   integer :: i, j, k

   vertex_gradient_z = u(i,   j,   k+1) - u(i,   j,   k) + &
                       u(i+1, j,   k+1) - u(i+1, j,   k) + &
                       u(i,   j+1, k+1) - u(i,   j+1, k) + &
                       u(i+1, j+1, k+1) - u(i+1, j+1, k)
   vertex_gradient_z = vertex_gradient_z / 4.0
end function vertex_gradient_z
