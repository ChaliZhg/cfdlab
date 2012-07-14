subroutine init_cond(rho, vex, vey, vez, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: pre(-1:nx+2, -1:ny+2, -1:nz+2)

   integer :: i, j, k


   !call init_cond_isen(rho, vex, vey, pre)
   !call init_cond_multvortdiff(0.0, rho, vex, vey, vez, pre)
   call init_cond_TaylorGreen(0.0, rho, vex, vey, vez, pre)
   !call init_cond2(rho, vex, vey, pre)

   return

   xmin = 0.0
   xmax = 1.0
   ymin = 0.0
   ymax = 1.0
   zmin = 0.0
   zmax = 1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny
   dz = (zmax - zmin)/nz

   do i=-1,nx+2
      do j=-1,ny+2
         do k=-1,nz+2
            rho(i,j,k) = 1.0
            vex(i,j,k) = 1.0
            vey(i,j,k) = 1.0
            vez(i,j,k) = 1.0
            pre(i,j,k) = 1.0
         enddo
      enddo
   enddo

end subroutine init_cond
