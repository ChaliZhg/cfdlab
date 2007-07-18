C Calculates the flux divergence, which is also the steady state residual
      subroutine divergence_kfvs(edge, ds, prim, coord, qx, qy, divf)
      implicit none
      include 'param.h'
      integer          edge(2,nemax)
      double precision ds(2,nemax), prim(nvar,npmax), 
     &                 coord(2,npmax), qx(nvar,npmax), qy(nvar,npmax), 
     &                 divf(nvar,npmax)

      integer          iv, ie, e1, e2
      double precision flux(nvar), priml(nvar), primr(nvar)

      do ie=1,ne
            e1 = edge(1,ie)
            e2 = edge(2,ie)
            call recon(e1, e2, coord, prim, qx, qy, priml, primr)
            call kfvs(ds(1,ie), ds(2,ie), priml, primr, flux)
            do iv=1,nvar-1
                  divf(iv,e1) = divf(iv,e1) + flux(iv)
                  divf(iv,e2) = divf(iv,e2) - flux(iv)
            enddo
      enddo

      return
      end

C.....Kinetic split fluxes
      subroutine kfvs(nx, ny, priml, primr, flux)
      implicit none
      include 'param.h'
      double precision nx, ny, priml(nvar), primr(nvar), flux(nvar)

      integer          i
      double precision rl, ul, vl, pl, el, rr, ur, vr, pr, er,
     &                 ql2, qr2, unl, unr, vnl, vnr,
     &                 length, ct, st, Al, Bl, Ar, Br,
     &                 sl, betal, sr, betar, ERRF,
     &                 Fp(4), Fm(4), Ff(4)

      length= dsqrt(nx**2 + ny**2)
      ct    = nx/length
      st    = ny/length

C     Left state
      rl = priml(1)
      ul = priml(2)
      vl = priml(3)
      pl = priml(4)
      ql2= ul**2 + vl**2
      el = pl/GAMMA1 + 0.5d0*rl*ql2

C     Right state
      rr = primr(1)
      ur = primr(2)
      vr = primr(3)
      pr = primr(4)
      qr2= ur**2 + vr**2
      er = pr/GAMMA1 + 0.5d0*rr*qr2

C     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

      vnl =-ul*st + vl*ct
      vnr =-ur*st + vr*ct

c     Positive flux
      betal = 0.5d0*rl/pl
      sl    = unl*dsqrt(betal)
      Al    = 0.5d0*(1.0d0 + ERRF(sl))
      Bl    = 0.5d0*dexp(-sl**2)/dsqrt(M_PI*betal)

      Fp(1) = rl*(unl*Al + Bl)
      Fp(2) = (pl + rl*unl**2)*Al + rl*unl*Bl
      Fp(3) = rl*(unl*Al + Bl)*vnl
      Fp(4) = (el + pl)*unl*Al + (el + 0.5d0*pl)*Bl

c     Negative flux
      betar = 0.5d0*rr/pr
      sr    = unr*dsqrt(betar)
      Ar    = 0.5d0*(1.0d0 - ERRF(sr))
      Br    = 0.5d0*dexp(-sr**2)/dsqrt(M_PI*betar)

      Fm(1) = rr*(unr*Ar - Br)
      Fm(2) = (pr + rr*unr**2)*Ar - rr*unr*Br
      Fm(3) = rr*(unr*Ar - Br)*vnr
      Fm(4) = (er + pr)*unr*Ar - (er + 0.5d0*pr)*Br

c     Total flux
      do i=1,4
            Ff(i) = Fp(i) + Fm(i)
      enddo

      flux(1) = Ff(1)*length
      flux(2) = (ct*Ff(2) - st*Ff(3))*length
      flux(3) = (st*Ff(2) + ct*Ff(3))*length
      flux(4) = Ff(4)*length

      return
      end
