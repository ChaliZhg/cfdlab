program main

   use comvar

   implicit none

   real, dimension(:), allocatable :: rho, vex, vey, pre, co0, co1, res
   real, dimension(:), allocatable :: phi, psi

   nx = 51
   ny = 51

   xmin =-5.0
   xmax = 5.0
   ymin =-5.0
   ymax = 5.0

   itmax = 10000
   itsave= 100

   ! options: iroe, irusanov
   fluxtype = iroe

   ! limiter: ford, muscl3, mmod
   limtype = muscl3

   ! fvm or gmd
   scheme = fvm

   fileid = 0

   dx = (xmax - xmin)/(nx-1)
   dy = (ymax - ymin)/(ny-1)

   cfl = 0.4

   ark(1) = 0.0
   ark(2) = 3.0/4.0
   ark(3) = 1.0/3.0

   allocate( rho( (nx+4)*(ny+4) ) )
   allocate( vex( (nx+4)*(ny+4) ) )
   allocate( vey( (nx+4)*(ny+4) ) )
   allocate( pre( (nx+4)*(ny+4) ) )
   allocate( co0( 4*(nx+4)*(ny+4) ) )
   allocate( co1( 4*(nx+4)*(ny+4) ) )

   allocate( phi(4*(nx+1)*(ny+1)) )
   allocate( psi(4*(nx+1)*(ny+1)) )
   allocate( res(4*(nx+2)*(ny+2)) )

   if(scheme==fvm)then
      call solveFVM(rho, vex, vey, pre, co0, co1, res)
   else if(scheme==gmd)then
      call solveGMD(rho, vex, vey, pre, co0, co1, phi, psi, res)
   else
      write(*,*)'Unknown scheme =',scheme
      stop
   endif

end program main
