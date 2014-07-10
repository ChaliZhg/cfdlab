program main

   use comvar

   implicit none

   real, dimension(:), allocatable :: pri, co0, co1, res, divB

   integer :: fid

   nx = 600
   ny = 600

   ! Default value, modified in initial condition function
   final_time = 10.0

   itmax = 500000
   itsave= 200

   ! Material properties
   gas_const = 1.0
   Cp        = gamma * gas_const / (gamma - 1.0)

   ! periodicity conditions
   xperiod = yes
   yperiod = yes

   ! options: 
   ! ient  = muscl type
   ! ifent = entropy stable
   fluxtype = ifent

   ! limiter: ford, muscl3, mmod
   limtype = mmod

   ! testcase: iot=orszag-tang, ikh=kelvin-helmholtz
   test_case = iot

   ! file id for saving solution
   fileid_sol = 0
   fileid_omg = 0

   cfl = 0.25

   nrk    = 3

   if(nrk.eq.2)then
      ark(1) = 0.0
      ark(2) = 0.5
   else if(nrk.eq.3)then
      ark(1) = 0.0
      ark(2) = 3.0/4.0
      ark(3) = 1.0/3.0
   endif
   brk = 1.0 - ark

   allocate( pri( nvar*(nx+4)*(ny+4) ) )
   allocate( co0( nvar*(nx+4)*(ny+4) ) )
   allocate( co1( nvar*(nx+4)*(ny+4) ) )
   allocate( divB( (nx+4)*(ny+4) ) )

   allocate( res(nvar*(nx+2)*(ny+2)) )

   call solveFVM(pri, co0, co1, res, divB)

end program main
