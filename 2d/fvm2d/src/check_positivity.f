C.....Check whether density, pressure, and turbulent viscosity are
C.....non-negative
      subroutine check_positivity(prim, coord, ptype)
      implicit none
      include 'param.h'
      integer          ptype(npmax)
      double precision prim(nvar,npmax), coord(2,npmax)

      integer          pt

      do pt=1,np

      if(prim(1,pt) .le. 0.0d0 .or. prim(4,pt) .le. 0.0d0 .or.
     &   prim(5,pt) .lt. 0.0d0)then
            print*,'Density/pressure is negative at point ', pt
            print*,'Type    = ', ptype(pt)
            print*,'Coord   = ', coord(1,pt), coord(2,pt)
            print*,'Density = ', prim(1,pt)
            print*,'u vel   = ', prim(2,pt)
            print*,'v vel   = ', prim(3,pt)
            print*,'Pressure= ', prim(4,pt)
            print*,'Visc    = ', prim(5,pt)
            stop
      endif

      enddo

      return
      end
