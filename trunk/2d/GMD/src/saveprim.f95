subroutine saveprim(t, rho, vex, vey, pre)
   use comvar
   implicit none

   real    :: t
   real    :: rho(-1:nx+2, -1:ny+2)
   real    :: vex(-1:nx+2, -1:ny+2)
   real    :: vey(-1:nx+2, -1:ny+2)
   real    :: pre(-1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: x, y, q, a, m
   character(len=512) :: filename

   filename = 'sol'
   call getfilename(filename, fileid_sol)

   open(10,file=trim(filename))
   write(10,*)'TITLE = "vortex flow"'
   write(10,*)'VARIABLES = "x", "y", "Density", "Velx", "Vely", "Pressure", "Mach"'
   write(10,*)'ZONE STRANDID=1, SOLUTIONTIME=',t,', I=',nx,', J=',ny,&
              ', DATAPACKING=POINT'

   do j=1,ny
      do i=1,nx
         x = xmin + (i-1)*dx + 0.5*dx
         y = ymin + (j-1)*dy + 0.5*dy

         q = sqrt(vex(i,j)**2 + vey(i,j)**2)
         a = sqrt(gamma*pre(i,j)/rho(i,j))
         m = q/a
         write(10,'(7E24.14)') x, y, rho(i,j), vex(i,j), vey(i,j), pre(i,j), m

      enddo
   enddo

   close(10)

end subroutine saveprim
