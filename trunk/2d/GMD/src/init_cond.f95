subroutine init_cond(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)


   call init_cond_isen(rho, vex, vey, pre)
   !call init_cond2(rho, vex, vey, pre)

end subroutine init_cond
