C.....Convert primitive to conserved variables
      subroutine prim_to_con(pt, prim, con)
      implicit none
      include 'param.h'
      integer          pt
      double precision prim(nvar,npmax), con(nvar)

      con(1) = prim(1,pt)
      con(2) = prim(1,pt)*prim(2,pt)
      con(3) = prim(1,pt)*prim(3,pt)
      con(4) = prim(4,pt)/GAMMA1 
     &        + 0.5d0*prim(1,pt)*(prim(2,pt)**2 + prim(3,pt)**2)

      return
      end

C.....Convert conserved to primitive variables
      subroutine con_to_prim(con, prim)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), con(nvar,npmax)

      integer          pt

      do pt=1,np
            prim(1,pt) = con(1,pt)
            prim(2,pt) = con(2,pt)/con(1,pt)
            prim(3,pt) = con(3,pt)/con(1,pt)
            prim(4,pt) = GAMMA1*(con(4,pt) -
     &                   0.5d0*(con(2,pt)**2+con(3,pt)**2)/con(1,pt))
      enddo

      return
      end
