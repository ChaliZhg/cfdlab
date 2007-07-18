C.....Calculate L2 residue based on density
C.....Normalize with residue of first iteration, res1
      subroutine residue(prim, prim_old, dt)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), prim_old(nvar,npmax), dt(npmax)

      integer          i
      double precision gra, dif

      res   = 0.0d0
      resi  = 0.0d0
      iresi = 0
      do i=1,np
            gra = prim(1,i) - prim_old(1,i)
            res = res + gra**2
            dif = dabs(gra)
            if(dif .gt. resi)then
                  resi = dif
                  iresi= i
            endif
      enddo
      res = dsqrt(res/np)

      if(iter .eq. iterlast)then
            res1 = res
            print*,'Residue in first iteration = ',res1
            if(res1 .eq. 0.0d0) stop
      endif

      res = res/res1

      return
      end
