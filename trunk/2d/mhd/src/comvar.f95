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

   real,parameter :: PI = 4.0*atan(1.0)

   real    :: gamma, g1, g2, g3, g4

   integer :: fileid_sol, fileid_omg

   integer :: fluxtype
   integer :: ient=1, ifent=2

   integer :: limtype
   integer :: ford=0, muscl3=1, mmod=2

   integer :: scheme
   integer :: fvm=1, gmd=2, kep=3, mvf=4

   integer :: no=0, yes=1

   integer :: xperiod, yperiod

   integer :: iot=1, ikh=2, irotor=3, ialfven=4
   integer :: test_case

end module comvar
