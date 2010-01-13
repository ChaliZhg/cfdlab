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

   itmax = 5000
   itsave= 100

   ! options: iroe, irusanov
   fluxtype = iroe
   !fluxtype = irusanov

   ! limiter: nolim, mmod
   limtype = nolim

   fileid = 0

   dx = (xmax - xmin)/(nx-1)
   dy = (ymax - ymin)/(ny-1)

   cfl = 0.9

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
   allocate( res(4*nx*ny) )

   call solve(rho, vex, vey, pre, co0, co1, phi, psi, res)

end program main
