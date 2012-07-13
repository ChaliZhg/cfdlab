subroutine kep_flux(ct, st, conl, conr, flux, dflux)
      use comvar
      implicit none

      real :: ct, st, conl(4), conr(4), flux(4), dflux(4)

      integer :: i
      real    :: rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr, &
                 ql2, qr2, un, Fc(4), Fd(4), &
                 r, u, v, p, h

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

      r  = 0.5*(rl + rr)
      u  = 0.5*(ul + ur)
      v  = 0.5*(vl + vr)
      p  = 0.5*(pl + pr)
      h  = 0.5*(hl + hr)

!     Rotated velocity
      un = u*ct + v*st

!     Centered flux
      Fc(1) = r*un
      Fc(2) = p*ct + u*Fc(1)
      Fc(3) = p*st + v*Fc(1)
      Fc(4) = h * Fc(1)

!     Total flux
      flux    =  0.5*Fc
      dflux   =  0.0

end subroutine kep_flux
