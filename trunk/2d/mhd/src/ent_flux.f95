!------------------------------------------------------------------------------
! Numerical flux function
! exactly entropy consistent
! left, right, left1, right1 are primitive variables (rho,u1,u2,u3,p,B1,B2,B3)
!------------------------------------------------------------------------------
subroutine ent_flux(lx, ly, left, right, left1, right1, flux)
   use comvar
   implicit none

   real :: lx, ly
   real :: left(nvar), right(nvar), left1(nvar), right1(nvar), flux(nvar)

   real :: rho, u1, u2, u3, ul2, ur2, q2, B1, B2, B3, Bl2, Br2, mB2
   real :: betal, betar, beta, rho_a, beta_a, p, bu1, bu2, bu3
   real :: a, srho, bb1, bb2, bb3, mbb2, cf2, cs2, cf, cs, alpf, alps
   real :: n1, n2, n3, np1, np2, np3, bet, alp, npn1, npn2, npn3, ff
   real :: bn, bnp, s1, s2, s3, sbn, t1, t2, t3
   real :: Kin, Jac(nvar,nvar), Rp(nvar,nvar), R(nvar,nvar), Lambda(nvar)
   real :: sl, sr, ds, ql2, qr2, dv(nvar), LRT(nvar,nvar), Diff
   integer :: i, j, k
   real :: unorm, Bnorm, bunorm
   real :: logavg

   rho = logavg (left(1), right(1))

   u1   = 0.5 * (left(2) + right(2))
   u2   = 0.5 * (left(3) + right(3))
   u3   = 0.5 * (left(4) + right(4))
   unorm= u1*lx + u2*ly

   ul2 = left(2)*left(2) + left(3)*left(3) + left(4)*left(4)
   ur2 = right(2)*right(2) + right(3)*right(3) + right(4)*right(4)
   q2  = 0.5 * (ul2 + ur2)

   B1   = 0.5 * (left(6) + right(6))
   B2   = 0.5 * (left(7) + right(7))
   B3   = 0.5 * (left(8) + right(8))
   Bnorm= B1*lx + B2*ly

   Bl2 = left(6)*left(6) + left(7)*left(7) + left(8)*left(8)
   Br2 = right(6)*right(6) + right(7)*right(7) + right(8)*right(8)
   mB2  = 0.5*(Bl2 + Br2)
   
   betal = left(1)/(2.0*left(5))
   betar = right(1)/(2.0*right(5))
   beta  = logavg(betal, betar)
   
   rho_a = 0.5*(left(1)+right(1))
   beta_a = 0.5*(betal+betar)
   p   = 0.5 * rho_a / beta_a

   bu1   = (betal*left(2) + betar*right(2))/(betal+betar)
   bu2   = (betal*left(3) + betar*right(3))/(betal+betar)
   bu3   = (betal*left(4) + betar*right(4))/(betal+betar)
   bunorm= bu1*lx + bu2*ly

   flux(1) = rho * unorm
   flux(2) = (p + 0.5*mB2)*lx + u1 * flux(1) - Bnorm * B1
   flux(3) = (p + 0.5*mB2)*ly + u2 * flux(1) - Bnorm * B2
   flux(4) =                    u3 * flux(1) - Bnorm * B3
   flux(6) = bunorm * B1 - bu1 * Bnorm
   flux(7) = bunorm * B2 - bu2 * Bnorm
   flux(8) = bunorm * B3 - bu3 * Bnorm
   flux(5) = (1.0/(2.0*(GAMMA-1.0)*beta) - 0.5*q2) * flux(1) &
      + u1 * flux(2) + u2 * flux(3) + u3 * flux(4) &
      + B1 * flux(6) + B2 * flux(7) + B3 * flux(8) &
      - 0.5*unorm*mB2 + (u1*B1+u2*B2+u3*B3)*Bnorm
      
   ! Add entropy dissipation
   ! normal vector
   n1 = lx
   n2 = ly
   n3 = 0.0

   ! Add entropy dissipation
   a = sqrt(0.5 * GAMMA / beta)
   
   srho = sqrt(rho)
   bb1 = B1/srho
   bb2 = B2/srho
   bb3 = B3/srho
   mbb2 = bb1*bb1 + bb2*bb2 + bb3*bb3
   bn  = bb1*n1  + bb2*n2  + bb3*n3
   cf2 = 0.5*(a*a + mbb2) + 0.5*sqrt( (a*a+mbb2)*(a*a+mbb2) - 4.0*a*a*bn*bn)
   cs2 = 0.5*(a*a + mbb2) - 0.5*sqrt( (a*a+mbb2)*(a*a+mbb2) - 4.0*a*a*bn*bn)
   cf = sqrt(cf2)
   cs = sqrt(abs(cs2)) ! cs may be zero or close to zero
   
   alpf = sqrt( abs((a*a - cs2)/(cf2 - cs2)) )
   alps = sqrt( abs((cf2 - a*a)/(cf2 - cs2)) )
   
   ! vector nperp
   ff = abs(mbb2 - bn*bn)
   if(ff < 1.0e-10) then
      if(abs(ly) < 1.0e-14)then
         np1 = 0.0
         np2 = 1.0/sqrt(2.0)
         np3 = 1.0/sqrt(2.0)
      else
         np1 = 1.0/sqrt(2.0)
         np2 = 0.0
         np3 = 1.0/sqrt(2.0)
      endif
   else
      bet = 1.0/sqrt(ff)
      alp = - bet * bn
      np1 = alp*n1 + bet*bb1
      np2 = alp*n2 + bet*bb2
      np3 = alp*n3 + bet*bb3
   endif
   
   ! Primitive eigenvectors
   
   ! entropy wave
   Rp(1,1) = g1*srho
   Rp(2,1) = 0
   Rp(3,1) = 0
   Rp(4,1) = 0
   Rp(5,1) = 0
   Rp(6,1) = 0
   Rp(7,1) = 0
   Rp(8,1) = 0
   
   ! divergence wave
   Rp(1,2) = 0
   Rp(2,2) = 0
   Rp(3,2) = 0
   Rp(4,2) = 0
   Rp(5,2) = 0
   Rp(6,2) = g2*a
   Rp(7,2) = 0
   Rp(8,2) = 0

   ! alfven waves
   ! np x n
   npn1 = np2*n3 - np3*n2
   npn2 = np3*n1 - np1*n3
   npn3 = np1*n2 - np2*n1

   Rp(1,3) = 0
   Rp(2,3) = g3*a*npn1/srho
   Rp(3,3) = g3*a*npn2/srho
   Rp(4,3) = g3*a*npn3/srho
   Rp(5,3) = 0
   Rp(6,3) = g3*a*npn1
   Rp(7,3) = g3*a*npn2
   Rp(8,3) = g3*a*npn3

   Rp(1,4) =  Rp(1,3)
   Rp(2,4) = -Rp(2,3)
   Rp(3,4) = -Rp(3,3)
   Rp(4,4) = -Rp(4,3)
   Rp(5,4) =  Rp(5,3)
   Rp(6,4) =  Rp(6,3)
   Rp(7,4) =  Rp(7,3)
   Rp(8,4) =  Rp(8,3)

   ! fast magneto acoustic wave
   bnp = bb1*np1 + bb2*np2 + bb3*np3
   s1  = (alpf*a*a*n1 + alps*a*(bnp*n1 - bn*np1))/(srho*cf)
   s2  = (alpf*a*a*n2 + alps*a*(bnp*n2 - bn*np2))/(srho*cf)
   s3  = (alpf*a*a*n3 + alps*a*(bnp*n3 - bn*np3))/(srho*cf)
   
   Rp(1,5) =  g3*alpf*srho
   Rp(2,5) = -g3*s1
   Rp(3,5) = -g3*s2
   Rp(4,5) = -g3*s3
   Rp(5,5) =  g3*alpf*srho*a*a
   Rp(6,5) =  g3*alps*a*np1
   Rp(7,5) =  g3*alps*a*np2
   Rp(8,5) =  g3*alps*a*np3
   
   Rp(1,6) =  Rp(1,5)
   Rp(2,6) = -Rp(2,5)
   Rp(3,6) = -Rp(3,5)
   Rp(4,6) = -Rp(4,5)
   Rp(5,6) =  Rp(5,5)
   Rp(6,6) =  Rp(6,5)
   Rp(7,6) =  Rp(7,5)
   Rp(8,6) =  Rp(8,5)
   
   ! slow magneto acoustic waves
   if(bn.gt.0)then
      sbn = +1.0
   else
      sbn = -1.0
   endif
   t1 = sbn*(alps*a*bn*n1 + alpf*cf*cf*np1)/(srho*cf)
   t2 = sbn*(alps*a*bn*n2 + alpf*cf*cf*np2)/(srho*cf)
   t3 = sbn*(alps*a*bn*n3 + alpf*cf*cf*np3)/(srho*cf)
   
   Rp(1,7) =  g3*alps*srho
   Rp(2,7) = -g3*t1
   Rp(3,7) = -g3*t2
   Rp(4,7) = -g3*t3
   Rp(5,7) =  g3*alps*srho*a*a
   Rp(6,7) = -g3*alpf*a*np1
   Rp(7,7) = -g3*alpf*a*np2
   Rp(8,7) = -g3*alpf*a*np3

   Rp(1,8) =  Rp(1,7)
   Rp(2,8) = -Rp(2,7)
   Rp(3,8) = -Rp(3,7)
   Rp(4,8) = -Rp(4,7)
   Rp(5,8) =  Rp(5,7)
   Rp(6,8) =  Rp(6,7)
   Rp(7,8) =  Rp(7,7)
   Rp(8,8) =  Rp(8,7)
   
   ! Jacobian = d(con)/d(prim)
   Kin = 0.5*(u1*u1 + u2*u2 + u3*u3)

   Jac    = 0.0

   Jac(1,1) = 1.0
   Jac(2,1) = u1
   Jac(3,1) = u2
   Jac(4,1) = u3
   Jac(5,1) = Kin

   Jac(2,2) = rho
   Jac(5,2) = rho*u1

   Jac(3,3) = rho
   Jac(5,3) = rho*u2

   Jac(4,4) = rho
   Jac(5,4) = rho*u3

   Jac(5,5) = g4

   Jac(5,6) = B1
   Jac(5,7) = B2
   Jac(5,8) = B3

   Jac(6,6) = 1.0
   Jac(7,7) = 1.0
   Jac(8,8) = 1.0

   ! Conserved eigenvectors: R = Jac * Rp
   do i=1,nvar
      do j=1,nvar
         R(i,j) = 0.0
         do k=1,nvar
            R(i,j) = R(i,j) + Jac(i,k)*Rp(k,j)
         enddo
      enddo
   enddo
   
   Lambda(1) = abs(unorm)
   Lambda(2) = abs(unorm)
   Lambda(3) = abs(unorm-bn)
   Lambda(4) = abs(unorm+bn)
   Lambda(5) = abs(unorm-cf)
   Lambda(6) = abs(unorm+cf)
   Lambda(7) = abs(unorm-cs)
   Lambda(8) = abs(unorm+cs)
   
   sl = log(left1(5)) - GAMMA * log(left1(1))
   sr = log(right1(5)) - GAMMA * log(right1(1))
   ds = sr - sl
   betal = left1(1)/(2.0*left1(5))
   betar = right1(1)/(2.0*right1(5))
   
   ! Jump in entropy variables
   ql2 = left1(2)*left1(2) + left1(3)*left1(3) + left1(4)*left1(4)
   qr2 = right1(2)*right1(2) + right1(3)*right1(3) + right1(4)*right1(4)
   dv(1) =-ds/(GAMMA-1.0) - (betar*qr2 - betal*ql2)
   dv(2) = 2.0*(betar*right1(2) - betal*left1(2))
   dv(3) = 2.0*(betar*right1(3) - betal*left1(3))
   dv(4) = 2.0*(betar*right1(4) - betal*left1(4))
   dv(5) =-2.0*(betar-betal)
   dv(6) = 2.0*(betar*right1(6) - betal*left1(6))
   dv(7) = 2.0*(betar*right1(7) - betal*left1(7))
   dv(8) = 2.0*(betar*right1(8) - betal*left1(8))
   
   ! LRT = L * R^T
   do i=1,nvar
      do j=1,nvar
         LRT(i,j) = Lambda(i)*R(j,i)
      enddo
   enddo
   
   do i=1,nvar
      Diff = 0.0
      do j=1,nvar
         do k=1,nvar
            Diff = Diff + R(i,j) * LRT(j,k) * dv(k)
         enddo
      enddo
      
      flux(i) = flux(i) - 0.5*Diff
   enddo

end subroutine ent_flux
