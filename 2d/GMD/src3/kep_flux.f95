subroutine kep_flux(cx, cy, cz, conl, conr, flux)
      use comvar
      implicit none

      real :: cx, cy, cz, conl(nvar), conr(nvar), flux(nvar)

      integer :: i
      real    :: rl, ul, vl, wl, pl, al2, hl, rr, ur, vr, wr, pr, ar2, hr, &
                 ql2, qr2, un, r, u, v, w, p, h

      rl = conl(1)
      ul = conl(2)/conl(1)
      vl = conl(3)/conl(1)
      wl = conl(4)/conl(1)
      pl = (GAMMA-1.0)*(conl(5) - 0.5*(conl(2)**2 + conl(3)**2 + conl(4)**2)/conl(1))
      ql2= ul**2 + vl**2 + wl**2
      al2= GAMMA*pl/rl
      hl = al2/(GAMMA-1.0) + 0.5*ql2

      rr = conr(1)
      ur = conr(2)/conr(1)
      vr = conr(3)/conr(1)
      wr = conr(4)/conr(1)
      pr = (GAMMA-1.0)*(conr(5) - 0.5*(conr(2)**2 + conr(3)**2 + conr(4)**2)/conr(1))
      qr2= ur**2 + vr**2 + wr**2
      ar2= GAMMA*pr/rr
      hr = ar2/(GAMMA-1.0) + 0.5*qr2

      r  = 0.5*(rl + rr)
      u  = 0.5*(ul + ur)
      v  = 0.5*(vl + vr)
      w  = 0.5*(wl + wr)
      p  = 0.5*(pl + pr)
      h  = 0.5*(hl + hr)

!     Rotated velocity
      un = u*cx + v*cy + w*cz

!     Centered flux
      flux(1) = r*un
      flux(2) = p*cx + u*flux(1)
      flux(3) = p*cy + v*flux(1)
      flux(4) = p*cz + w*flux(1)
      flux(5) = h * flux(1)

end subroutine kep_flux
