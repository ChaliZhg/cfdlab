module comvar
   implicit none

   integer :: nvar = 8
   integer :: nx, ny
   real    :: xmin, xmax, ymin, ymax, dx, dy, dt
   integer :: itmax
   integer :: itsave
   real    :: cfl
   real    :: final_time

   real    :: gas_const, Cp

   integer :: nrk
   real    :: ark(3), brk(3)

   real,parameter :: gamma = 5.0/3.0
   real :: g1 = sqrt( (gamma-1.0)/gamma )
   real :: g2 = sqrt( 1.0/gamma )
   real :: g3 = sqrt( 0.5/gamma )
   real :: g4 = sqrt( 1.0/(gamma-1.0) )
   real :: PI = 4.0*atan(1.0)

   integer :: fileid_sol, fileid_omg

   integer :: fluxtype
   integer :: ient=1, ifent=2

   integer :: limtype
   integer :: ford=0, muscl3=1, mmod=2

   integer :: scheme
   integer :: fvm=1, gmd=2, kep=3, mvf=4

   integer :: no=0, yes=1

   integer :: xperiod, yperiod

   integer :: iot=1, ikh=2
   integer :: test_case

end module comvar
