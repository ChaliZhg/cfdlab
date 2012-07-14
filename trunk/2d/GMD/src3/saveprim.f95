subroutine saveprim(t, rho, vex, vey, vez, pre)
   use comvar
   implicit none

   real    :: t
   real    :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real    :: pre(-1:nx+2, -1:ny+2, -1:nz+2)

   integer :: i, j, k
   real    :: x, y, z, q, a, m
   character(len=512) :: filename

   filename = 'sol'
   call getfilename(filename, fileid_sol)

   open(10,file=trim(filename))
   write(10,*)'TITLE = "vortex flow"'
   write(10,*)'VARIABLES = "x", "y", "z", "Density", "U", "V", "W", "Pressure", "Mach"'
   write(10,*)'ZONE STRANDID=1, SOLUTIONTIME=',t,', I=',nx,', J=',ny,&
              ', K=',nz,', DATAPACKING=POINT'

   do k=1,nz
      do j=1,ny
         do i=1,nx

            x = xmin + (i-1)*dx + 0.5*dx
            y = ymin + (j-1)*dy + 0.5*dy
            z = zmin + (k-1)*dz + 0.5*dz

            q = sqrt(vex(i,j,k)**2 + vey(i,j,k)**2 + vez(i,j,k)**2)
            a = sqrt(gamma*pre(i,j,k)/rho(i,j,k))
            m = q/a
            write(10,'(7E24.14)') x, y, z, rho(i,j,k), vex(i,j,k), vey(i,j,k), vez(i,j,k), pre(i,j,k), m

         enddo
      enddo
   enddo

   close(10)

end subroutine saveprim
