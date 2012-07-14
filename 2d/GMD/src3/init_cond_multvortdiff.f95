! From HIOCFD 2012 workshop
subroutine init_cond_multvortdiff(time, rho, vex, vey, vez, pre)
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

   xmin = -M_PI
   xmax =  M_PI
   ymin = -M_PI
   ymax =  M_PI
   zmin = -M_PI
   zmax =  M_PI

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny
   dz = (zmax - zmin)/nz

   rho0 = 1.0
   V0   = 1.0
   Mach = 0.1
   c0   = V0 / Mach
   T0   = c0*c0/(gas_const * GAMMA)
   p0   = 50.0

   do i=-1,nx+2
      do j=-1,ny+2
         do k=-1,nz+2
            x = xmin + (i-1)*dx + 0.5*dx
            y = ymin + (j-1)*dy + 0.5*dy
            z = zmin + (k-1)*dz + 0.5*dz

            vex(i,j,k) =  V0 * sin(x) * cos(y) * cos(z)
            vey(i,j,k) = -V0 * cos(x) * sin(y) * cos(z)
            vez(i,j,k) =  0.0
            pre(i,j,k) =  p0 + (1.0/16.0) * rho0 * V0 * V0 * (cos(2.0*x) + cos(2.0*y))*(cos(2.0*z) + 2.0)
            rho(i,j,k) =  pre(i,j,k)/(gas_const * T0)
         enddo
      enddo
   enddo

end subroutine init_cond_multvortdiff
