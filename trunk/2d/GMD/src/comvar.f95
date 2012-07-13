module comvar
   implicit none

   integer :: nx, ny
   real    :: xmin, xmax, ymin, ymax, dx, dy, dt
   integer :: itmax
   integer :: itsave
   real    :: cfl
   real    :: final_time

   real    :: mu, Prandtl, gas_const, kth, Cp

   integer :: nrk
   real    :: ark(3)

   real :: gamma = 1.4
   real :: M_PI = 4.0*atan(1.0)

   integer :: fileid_sol, fileid_omg

   integer :: fluxtype
   integer :: iroe=1, irusanov=2, iadv=3, icusp=4, ikep=5, ikepes=6

   integer :: ikepes_diss

   integer :: limtype
   integer :: ford=0, muscl3=1, mmod=2

   integer :: scheme
   integer :: fvm=1, gmd=2, kep=3, mvf=4

   integer :: vconf

   integer :: no=0, yes=1

   integer :: xperiod, yperiod

end module comvar
