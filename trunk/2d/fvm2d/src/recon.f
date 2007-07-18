c------------------------------------------------------------------------------
c Limited reconstruction using van-albada limiter
c------------------------------------------------------------------------------
      subroutine recon(e1, e2, coord, prim, qx, qy, priml, primr)
      implicit none
      include 'param.h'
      integer          e1, e2
      double precision coord(2,npmax), prim(nvar,npmax), qx(nvar,npmax),
     &                 qy(nvar,npmax), priml(nvar), primr(nvar)

      integer          i
      double precision dx, dy, dqp, dql, si, dqr, sj, limit_albada

      dx = coord(1,e2) - coord(1,e1)
      dy = coord(2,e2) - coord(2,e1)

      if(ilimit .eq. yes)then

         do i=1,nvar-1
            dqp = prim(i,e2) - prim(i,e1)

            dql = ALBADA21*(dx*qx(i,e1) + dy*qy(i,e1)) + ALBADA22*dqp
            si  = limit_albada(dql, dqp)
            priml(i) = prim(i,e1) + 0.5d0*si

            dqr = ALBADA21*(dx*qx(i,e2) + dy*qy(i,e2)) + ALBADA22*dqp
            sj  = limit_albada( dqr, dqp )
            primr(i) = prim(i,e2) - 0.5d0*sj
         enddo

      else

         do i=1,nvar-1
            dqp = prim(i,e2) - prim(i,e1)

            dql = dx*qx(i,e1) + dy*qy(i,e1)
            si  = ALBADA11*dql + ALBADA12*dqp
            priml(i) = prim(i,e1) + 0.5d0*si

            dqr = dx*qx(i,e2) + dy*qy(i,e2)
            sj  = ALBADA11*dqr + ALBADA12*dqp
            primr(i) = prim(i,e2) - 0.5d0*sj
         enddo

      endif


      return
      end

c------------------------------------------------------------------------------
c van-albada limiter
c------------------------------------------------------------------------------
      double precision function limit_albada(ul, ur)
      implicit none
      include 'param.h'
      double precision ul, ur

      double precision top, bot

      if( ul*ur .le. 0.0d0)then
         limit_albada = 0.0d0
      else
         top = (ul**2 + EPSILON)*ur + (ur**2 + EPSILON)*ul
         bot = ul**2 + ur**2 + 2.0d0*EPSILON
         limit_albada = top/bot
      endif

      return
      end
