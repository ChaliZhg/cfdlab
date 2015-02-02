subroutine timestep(pri)
   use comvar
   use omp_lib
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)

   integer :: i, j
   real    :: c2, b2, f, bx, by, cfx, cfy, eigx, eigy

   dt = 1.0e20

   !$omp parallel do private(c2,b2,f,bx,cfx,eigx,by,cfy,eigy) shared(dt)
   do i=1,nx
      do j=1,ny
         c2= gamma*pri(5,i,j)/pri(1,i,j)
         b2= (pri(6,i,j)**2 + pri(7,i,j)**2 + pri(8,i,j)**2) / pri(1,i,j)
         f = c2 + b2

         bx= pri(6,i,j)/sqrt(pri(1,i,j))
         cfx = sqrt(0.5*(f + sqrt(f*f - 4.0*c2*bx*bx)))
         eigx= abs(pri(2,i,j)) + cfx

         by= pri(7,i,j)/sqrt(pri(1,i,j))
         cfy = sqrt(0.5*(f + sqrt(f*f - 4.0*c2*by*by)))
         eigy= abs(pri(3,i,j)) + cfy

         dt = min(dt, 1.0/(eigx/dx + eigy/dy))
      enddo
   enddo
   !$omp end parallel do

   dt = cfl*dt

   write(*,*)'Time step =', dt

end subroutine timestep
