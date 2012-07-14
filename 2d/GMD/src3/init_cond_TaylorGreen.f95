! Taylor-Green vortex
! See paper of Shu, Sun, Gottlieb
subroutine init_cond_TaylorGreen(time, rho, vex, vey, vez, pre)
   use comvar
   implicit none

   real    :: time
   real    :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: pre(-1:nx+2, -1:ny+2, -1:nz+2)

   integer :: i, j, k
   real    :: x, y, z, rho0, V0, Mach, c0, T0, p0

   final_time = 5.0

   xmin = 0.0
   xmax = 2.0*M_PI
   ymin = 0.0
   ymax = 2.0*M_PI
   zmin = 0.0
   zmax = 2.0*M_PI

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny
   dz = (zmax - zmin)/nz

   rho0 = 1.0
   p0   = 100.0

   do i=-1,nx+2
      do j=-1,ny+2
         do k=-1,nz+2
            x = xmin + (i-1)*dx + 0.5*dx
            y = ymin + (j-1)*dy + 0.5*dy
            z = zmin + (k-1)*dz + 0.5*dz

            rho(i,j,k) =  rho0
            vex(i,j,k) =  sin(x) * cos(y) * cos(z)
            vey(i,j,k) = -cos(x) * sin(y) * cos(z)
            vez(i,j,k) =  0.0
            pre(i,j,k) =  p0 + (rho0/16.0)*((cos(2.0*x) + cos(2.0*y))*(cos(2.0*z) + 2.0) - 2.0)
         enddo
      enddo
   enddo

end subroutine init_cond_TaylorGreen
