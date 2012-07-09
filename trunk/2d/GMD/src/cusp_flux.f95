subroutine cusp_flux(ct, st, conl, conr, flux, dflux)
      use comvar
      implicit none

      real :: ct, st, conl(4), conr(4), flux(4), dflux(4)

      integer :: i
      real    :: rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr, &
                 ua, va, qa2, aa2, aa, ha, &
                 ql2, qr2, rl12, rr12, rd, &
                 unl, unr, una, vna, Fc(4), Fd(4), ma, f1, f2

      rl = conl(1)
      ul = conl(2)/conl(1)
      vl = conl(3)/conl(1)
      pl = (GAMMA-1.0)*(conl(4) - 0.5*(conl(2)**2 + conl(3)**2)/conl(1))
      ql2= ul**2 + vl**2
      al2= GAMMA*pl/rl
      hl = al2/(GAMMA-1.0) + 0.5*ql2

      rr = conr(1)
      ur = conr(2)/conr(1)
      vr = conr(3)/conr(1)
      pr = (GAMMA-1.0)*(conr(4) - 0.5*(conr(2)**2 + conr(3)**2)/conr(1))
      qr2= ur**2 + vr**2
      ar2= GAMMA*pr/rr
      hr = ar2/(GAMMA-1.0) + 0.5*qr2

!     Rotated velocity
      unl = ul*ct + vl*st
      unr = ur*ct + vr*st

!     Centered flux
      Fc(1) = 0.5 * (rl*unl    + rr*unr)
      Fc(2) = 0.5 * (pl*ct + rl*ul*unl + pr*ct + rr*ur*unr)
      Fc(3) = 0.5 * (pl*st + rl*vl*unl + pr*st + rr*vr*unr)
      Fc(4) = 0.5 * (rl*hl*unl + rr*hr*unr)

!     Roe average
      rl12 = sqrt(rl)
      rr12 = sqrt(rr)
      rd   = 1.0/(rl12 + rr12)

      ua   = (ul*rl12 + ur*rr12)*rd
      va   = (vl*rl12 + vr*rr12)*rd
      ha   = (hl*rl12 + hr*rr12)*rd
      qa2  = ua**2 + va**2
      aa2  = (GAMMA-1.0)*(ha - 0.5*qa2)

      if(aa2 .le. 0.0)then
         print*,'Sonic speed is negative'
         print*,'Left/right conserved values'
         print*,conl(:)
         print*,conr(:)
         print*
         print*,'Left/right primitive values'
         print*,rl,ul,vl,pl
         print*,rr,ur,vr,pr
         stop
      endif

      aa  = sqrt(aa2)
      una = ua*ct + va*st
      ma  = una/aa
      if(abs(ma) > 1.0)then
         f1 = abs(ma)
         f2 = sign(1.0, ma)
      else
         f1 = 0.25 + ma**2 - 0.25*ma**4
         f2 = 0.5 * ma * (3.0 - ma**2)
      endif

!     Difference of conserved variables
      Fd(1) = 0.5 * f1 * aa * (rr    - rl)
      Fd(2) = 0.5 * f1 * aa * (rr*ur - rl*ul)
      Fd(3) = 0.5 * f1 * aa * (rr*vr - rl*vl)
      Fd(4) = 0.5 * f1 * aa * (rr*hr - rl*hl)

      dflux   = -Fd

      Fd(2) = Fd(2) +  0.5 * f2 * (pr - pl) * ct
      Fd(3) = Fd(3) +  0.5 * f2 * (pr - pl) * st

!     Total flux
      flux    =  Fc - Fd

end subroutine cusp_flux
