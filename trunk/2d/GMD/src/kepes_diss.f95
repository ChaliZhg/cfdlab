subroutine kepes_diss(ax, ay, conjm1, conj, conjp1, conjp2, dflux)
   use comvar
   implicit none

   real :: ax, ay
   real :: conjm1(4), conj(4), conjp1(4), conjp2(4), flux(4), dflux(4)
   real :: u, v, pjm1, pj, pjp1, pjp2, lam, nuj, nujp1, nu, eps2, eps4, drho
   real :: r, p, a, bl, br, beta, udotu
   real :: K2, K4
   real :: logavg

   K2 = 0.0
   K4 = 0.01

   u = 0.5 * (conj(2)/conj(1) + conjp1(2)/conjp1(1))
   v = 0.5 * (conj(3)/conj(1) + conjp1(3)/conjp1(1))
   r = 0.5 * (conj(1) + conjp1(1))

   pj= (gamma-1.0)*(conj(4) - 0.5*(conj(2)**2 + conj(3)**2)/conj(1))
   pjm1= (gamma-1.0)*(conjm1(4) - 0.5*(conjm1(2)**2 + conjm1(3)**2)/conjm1(1))
   pjp1= (gamma-1.0)*(conjp1(4) - 0.5*(conjp1(2)**2 + conjp1(3)**2)/conjp1(1))
   pjp2= (gamma-1.0)*(conjp2(4) - 0.5*(conjp2(2)**2 + conjp2(3)**2)/conjp2(1))
   
   p = 0.5*(pj + pjp1)


   bl    = 0.5 * conj(1) / pj
   br    = 0.5 * conjp1(1) / pjp1
   beta  = logavg(bl, br)
   udotu = conj(2)/conj(1) * conjp1(2)/conjp1(1) + &
           conj(3)/conj(1) * conjp1(3)/conjp1(1)

   a     = sqrt(gamma*p/r)
   lam   = abs(u)*ax + abs(v)*ay + a
   nuj   = abs(pjm1 - 2.0*pj + pjp1)/(pjm1 + 2.0*pj + pjp1)
   nujp1 = abs(pj - 2.0*pjp1 + pjp2)/(pj + 2.0*pjp1 + pjp2)
   nu    = max(nuj, nujp1)
   eps2  = min(1.0, K2 * nu)
   eps4  = max(0.0, K4 - eps2)

   drho  =   eps2 * (conjp1(1) - conj(1)) &
           - eps4 *(conjp2(1) - 3.0*conjp1(1) + 3.0*conj(1) - conjm1(1))


   dflux(1) = lam * drho
   dflux(2) = u * dflux(1)
   dflux(3) = v * dflux(1)
   dflux(4) = 0.5*( 1.0/((gamma-1.0)*beta) + udotu ) * dflux(1)


end subroutine kepes_diss
