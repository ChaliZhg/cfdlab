subroutine init_cond(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, r2, Temp, circ, circ1, circ2

   circ = 5.0
   circ1= (gamma-1.0)*circ**2/(8.0*gamma*M_PI**2)
   circ2= circ/(2.0*M_PI)

   do i=1,nx
      do j=1,ny
         x = xmin + (i-1)*dx
         y = ymin + (j-1)*dy

         r2   = x**2 + y**2
         Temp = 1.0 - circ1*exp(1.0-r2)

         rho(i,j) =  Temp**(1.0/(gamma-1.0))
         vex(i,j) = -circ2*y*exp(0.5*(1.0-r2))
         vey(i,j) =  circ2*x*exp(0.5*(1.0-r2))
         pre(i,j) =  rho(i,j)*Temp
      enddo
   enddo

end subroutine init_cond
