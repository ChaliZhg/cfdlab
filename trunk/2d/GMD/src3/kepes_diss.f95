subroutine kepes_diss(ax, ay, az, conjm1, conj, conjp1, conjp2, dflux)
   use comvar
   implicit none

   real :: ax, ay, az
   real :: conjm1(nvar), conj(nvar), conjp1(nvar), conjp2(nvar), dflux(nvar)
   real :: u, v, w, pjm1, pj, pjp1, pjp2, lam, nuj, nujp1, nu, eps2, eps4, drho
   real :: uj, ujp1, ujp2, ujm1
   real :: vj, vjp1, vjp2, vjm1
   real :: wj, wjp1, wjp2, wjm1
   real :: Tj, Tjp1, Tjp2, Tjm1
   real :: du, dv, dw, dTemp
   real :: r, a, bl, br, beta, udotu
   real :: K2, K4
   real :: logavg

   K2 = 0.5
   K4 = 0.01

   ujm1 = conjm1(2)/conjm1(1)
   uj   = conj(2)  /conj(1)
   ujp1 = conjp1(2)/conjp1(1)
   ujp2 = conjp2(2)/conjp2(1)

   vjm1 = conjm1(3)/conjm1(1)
   vj   = conj(3)  /conj(1)
   vjp1 = conjp1(3)/conjp1(1)
   vjp2 = conjp2(3)/conjp2(1)

   wjm1 = conjm1(4)/conjm1(1)
   wj   = conj(4)  /conj(1)
   wjp1 = conjp1(4)/conjp1(1)
   wjp2 = conjp2(4)/conjp2(1)

   pj  = (gamma-1.0)*(conj(5)   - 0.5*(conj(2)**2   + conj(3)**2   + conj(4)**2  )/conj(1))
   pjm1= (gamma-1.0)*(conjm1(5) - 0.5*(conjm1(2)**2 + conjm1(3)**2 + conjm1(4)**2)/conjm1(1))
   pjp1= (gamma-1.0)*(conjp1(5) - 0.5*(conjp1(2)**2 + conjp1(3)**2 + conjp1(4)**2)/conjp1(1))
   pjp2= (gamma-1.0)*(conjp2(5) - 0.5*(conjp2(2)**2 + conjp2(3)**2 + conjp2(4)**2)/conjp2(1))

   Tjm1 = pjm1/(gas_const * conjm1(1))
   Tj   = pj  /(gas_const * conj(1))
   Tjp1 = pjp1/(gas_const * conjp1(1))
   Tjp2 = pjp2/(gas_const * conjp2(1))
   
   u = 0.5 * (uj + ujp1)
   v = 0.5 * (vj + vjp1)
   w = 0.5 * (wj + wjp1)
   r = 0.5 * (conj(1) + conjp1(1))

   bl    = 0.5 * conj(1) / pj
   br    = 0.5 * conjp1(1) / pjp1
   beta  = logavg(bl, br)
   udotu = uj * ujp1 + vj * vjp1 + wj * wjp1

   a     = sqrt(0.5*gamma/beta)
   lam   = abs(u)*ax + abs(v)*ay + abs(w)*az + a
   nuj   = abs(pjm1 - 2.0*pj + pjp1)/(pjm1 + 2.0*pj + pjp1)
   nujp1 = abs(pj - 2.0*pjp1 + pjp2)/(pj + 2.0*pjp1 + pjp2)
   nu    = max(nuj, nujp1)
   eps2  = min(1.0, K2 * nu)
   eps4  = max(0.0, K4 - eps2)

   drho  =   eps2 * (conjp1(1) - conj(1)) &
           - eps4 *(conjp2(1) - 3.0*conjp1(1) + 3.0*conj(1) - conjm1(1))
   du    =   eps2 * (ujp1 - uj) - eps4 *(ujp2 - 3.0*ujp1 + 3.0*uj - ujm1)
   dv    =   eps2 * (vjp1 - vj) - eps4 *(vjp2 - 3.0*vjp1 + 3.0*vj - vjm1)
   dw    =   eps2 * (wjp1 - wj) - eps4 *(wjp2 - 3.0*wjp1 + 3.0*wj - wjm1)
   dTemp =   eps2 * (Tjp1 - Tj) - eps4 *(Tjp2 - 3.0*Tjp1 + 3.0*Tj - Tjm1)

   !dflux(1) = lam * drho
   !dflux(2) = u * dflux(1)
   !dflux(3) = v * dflux(1)
   !dflux(4) = w * dflux(1)
   !dflux(5) = 0.5*( 1.0/((gamma-1.0)*beta) + udotu ) * dflux(1)

   dflux(1) = drho
   dflux(2) = u * drho + r * du
   dflux(3) = v * drho + r * dv
   dflux(4) = w * drho + r * dw
   dflux(5) = 0.5*( 1.0/((gamma-1.0)*beta) + udotu ) * drho + &
              r *(u*du + v*dv + w*dw) + 0.5 * r /(gamma-1.0) * dTemp

   dflux = 0.5 * lam * dflux

end subroutine kepes_diss
