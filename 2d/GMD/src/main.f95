program main

   use comvar

   implicit none

   real, dimension(:), allocatable :: rho, vex, vey, pre, co0, co1, res
   real, dimension(:), allocatable :: phi, psi
   real, dimension(:), allocatable :: phid, psid
   real, dimension(:), allocatable :: omg

   nx = 51
   ny = 51

   xmin =-5.0
   xmax = 5.0
   ymin =-5.0
   ymax = 5.0

   itmax = 10000
   itsave= 100

   ! periodicity conditions
   xperiod = yes
   yperiod = yes

   ! options: iroe, irusanov
   fluxtype = iroe

   ! limiter: ford, muscl3, mmod
   limtype = muscl3

   ! fvm or gmd
   scheme = fvm

   ! vorticity confinement
   vconf  = no

   ! file id for saving solution
   fileid_sol = 0
   fileid_omg = 0

   dx = (xmax - xmin)/(nx-1)
   dy = (ymax - ymin)/(ny-1)

   cfl = 0.8

   ark(1) = 0.0
   ark(2) = 3.0/4.0
   ark(3) = 1.0/3.0

   allocate( rho( (nx+4)*(ny+4) ) )
   allocate( vex( (nx+4)*(ny+4) ) )
   allocate( vey( (nx+4)*(ny+4) ) )
   allocate( pre( (nx+4)*(ny+4) ) )
   allocate( omg( (nx+1)*(ny+1) ) )
   allocate( co0( 4*(nx+4)*(ny+4) ) )
   allocate( co1( 4*(nx+4)*(ny+4) ) )

   allocate( phi(4*(nx+1)*(ny+1)) )
   allocate( psi(4*(nx+1)*(ny+1)) )
   allocate( phid(4*(nx+1)*(ny+1)) )
   allocate( psid(4*(nx+1)*(ny+1)) )
   allocate( res(4*(nx+2)*(ny+2)) )

   if(scheme==fvm)then
      call solveFVM(rho, vex, vey, pre, omg, co0, co1, res)
   else if(scheme==gmd)then
      call solveGMD(rho, vex, vey, pre, omg, co0, co1, phi, psi,  &
                    phid, psid, res)
   else
      write(*,*)'Unknown scheme =',scheme
      stop
   endif

end program main
