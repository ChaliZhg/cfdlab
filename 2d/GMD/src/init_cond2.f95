! Compact vortex
! Taken from AIAA-2001-0606
! M Murayama and K Nakahashi and S Obayashi
! Numerical simulation of vortical flows using vorticity confinement coupled
! with unstructured grid
subroutine init_cond2(rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: machinf, uinf, cinf, rhoinf, pinf
   real    :: Uc, Rc, R0, A, B
   real    :: x, y, r, theta, utheta

   xmin = 0.0
   xmax = 1.0
   ymin = 0.0
   ymax = 1.0

   dx = (xmax - xmin)/nx
   dy = (ymax - ymin)/ny

   machinf = 0.5
   uinf    = 1.0
   cinf    = uinf/machinf
   rhoinf  = 1.0
   pinf    = rhoinf*cinf**2/gamma

   Uc      = uinf
   Rc      = 0.05
   R0      = 10.0*Rc

   A       =-Uc*Rc/(R0**2 - Rc**2)
   B       = Uc*Rc*R0**2/(R0**2 - Rc**2)

   do i=-1,nx+2
      do j=-1,ny+2
         x = xmin + (i-1)*dx + 0.5*dx - 0.5
         y = ymin + (j-1)*dy + 0.5*dy - 0.5

         r    = sqrt(x**2 + y**2)
         theta= atan2(y, x)

         if(r < Rc) then
            utheta   = Uc*r/Rc
         else if(r >= Rc .and. r <= R0) then
            utheta   = A*r + B/r
         else
            utheta = 0.0
         endif

         rho(i,j) = rhoinf
         vex(i,j) =-utheta*sin(theta) + uinf
         vey(i,j) = utheta*cos(theta)
         pre(i,j) = pinf
      enddo
   enddo

end subroutine init_cond2
