c======================================================================
      subroutine meshless(dt, n, x0, y0, Q0, x1, y1, Q1, Qx1, Qy1)
c======================================================================
c n = number of points req meshless update
c dt(n) = local time-step
c x0(n),y0(n) = coordinates of pts req meshless update
c x1(nnbr,n),y1(nnbr,n) = coordinates of neighbouring points
c Q0(nvar,n)
c Q1(nvar,nnbr,n)
c Qx1(nvar,nnbr,n), Qy1(nvar,nnbr,n) = derivatives at neighbouring pts
c======================================================================
      implicit none
      include 'misc.h'
      integer n
      real    dt(*), x0(*), y0(*), Q0(nvar,*), x1(nnbr,*), y1(nnbr,*),
     +        Q1(nvar,nnbr,*), Qx1(nvar,nnbr,*), Qy1(nvar,nnbr,*)

c     local variables
      integer i

      do i=1,n
c        find coefficients
c        find q-derivatives
c        find xyflux
c        for each neighbour
c           find left-right states
c           find flux
c           update residue
c        update solution
      enddo

      return
      end
