subroutine init_cond(pri, co1)
   use comvar
   implicit none

   real    :: pri(nvar, -1:nx+2, -1:ny+2)
   real    :: co1(nvar, -1:nx+2, -1:ny+2)


   integer :: i, j
   real    :: q2, B2
   
   if(test_case == iot)then
      call orszag_tang(pri)
   elseif(test_case == ikh)then
      call kelvin_helmholtz(pri)
   elseif(test_case == irotor)then
      call rotor(pri)
   elseif(test_case == ialfven)then
      call alfven(pri)
   else
      print*,'Unknown test case !!!'
      stop
   endif

   ! Gas constants stored in comvar
   g1 = sqrt( (gamma-1.0)/gamma )
   g2 = sqrt( 1.0/gamma )
   g3 = sqrt( 0.5/gamma )
   g4 = sqrt( 1.0/(gamma-1.0) )

   do i=-1,nx+2
      do j=-1,ny+2
         ! Convert to conserved variables
         co1(1,i,j)   = pri(1,i,j)
         co1(2:4,i,j) = pri(2:4,i,j) * pri(1,i,j)
         q2 = pri(2,i,j)**2 + pri(3,i,j)**2 + pri(4,i,j)**2
         B2 = pri(6,i,j)**2 + pri(7,i,j)**2 + pri(8,i,j)**2
         co1(5,i,j)   = pri(5,i,j)/(gamma-1.0) + 0.5*pri(1,i,j)*q2 + 0.5*B2
         co1(6:8,i,j) = pri(6:8,i,j)
      enddo
   enddo

end subroutine init_cond
