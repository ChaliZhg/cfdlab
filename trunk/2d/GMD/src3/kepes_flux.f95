subroutine kepes_flux(cx, cy, cz, conl, conr, flux)
      use comvar
      implicit none

      real :: cx, cy, cz, conl(nvar), conr(nvar), flux(nvar)

      integer :: i
      real    :: rl, ul, vl, wl, pl, al2, hl, rr, ur, vr, wr, pr, ar2, hr, &
                 ql2, qr2, un, Fc(nvar), &
                 r, u, v, w, p, ra, ba, bl, br, beta, q2
      real    :: logavg

      rl = conl(1)
      ul = conl(2)/conl(1)
      vl = conl(3)/conl(1)
      wl = conl(4)/conl(1)
      pl = (GAMMA-1.0)*(conl(5) - 0.5*(conl(2)**2 + conl(3)**2 + conl(4)**2)/conl(1))
      ql2= ul**2 + vl**2 + wl**2
      al2= GAMMA*pl/rl
      hl = al2/(GAMMA-1.0) + 0.5*ql2
      bl = 0.5 * rl / pl

      rr = conr(1)
      ur = conr(2)/conr(1)
      vr = conr(3)/conr(1)
      wr = conr(4)/conr(1)
      pr = (GAMMA-1.0)*(conr(5) - 0.5*(conr(2)**2 + conr(3)**2 + conr(4)**2)/conr(1))
      qr2= ur**2 + vr**2 + wr**2
      ar2= GAMMA*pr/rr
      hr = ar2/(GAMMA-1.0) + 0.5*qr2
      br = 0.5 * rr / pr

      r  = logavg(rl, rr)
      beta = logavg(bl, br)
      u  = 0.5*(ul + ur)
      v  = 0.5*(vl + vr)
      w  = 0.5*(wl + wr)
      q2 = 0.5*(ql2 + qr2)

      ra = 0.5*(rl + rr)
      ba = 0.5*(bl + br)
      p  = 0.5 * ra / ba

!     Rotated velocity
      un = u*cx + v*cy + w*cz

!     Centered flux
      Fc(1) = r*un
      Fc(2) = p*cx + u*Fc(1)
      Fc(3) = p*cy + v*Fc(1)
      Fc(4) = p*cz + w*Fc(1)
      Fc(5) = 0.5 * (1.0/((gamma-1.0)*beta) - q2) * Fc(1) + u*Fc(2) + v*Fc(3) + w*Fc(4)

      flux  = Fc

end subroutine kepes_flux
