      subroutine costfun(nc,q,cost,dx)
      implicit none
      integer :: nc
      real    :: q(nc),cost,dx(nc)

      integer :: i

      cost = 0.0
      do i=1,nc
         cost = cost + q(i)**2 * dx(i)
      enddo

      print*,'Cost =',cost

      return
      end
