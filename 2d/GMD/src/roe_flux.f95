subroutine roe_flux(ct, st, conl, conr, flux, dflux)
      use comvar
      implicit none

      real :: ct, st, conl(4), conr(4), flux(4), dflux(4)

      integer :: i
      real    :: rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr, &
                 ua, va, qa2, aa2, aa, ha, &
                 ql2, qr2, rl12, rr12, rd, &
                 unl, unr, una, vna, Fc(4), Fd(4), &
                 m1, m2, a1, a2, a3, a4, l1, l2, l3, l4, &
                 a1l1, a2l2, a3l3, a4l4, aact, aast, &
                 du1, du2, du3, du4, dl, dr, li(4), limit, &
                 e1, e4, del, ETOL

      ETOL = 0.01 ! parameter for entropy fix
      ETOL = 0.00 ! parameter for entropy fix

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
      Fc(1) = rl*unl            + rr*unr
      Fc(2) = pl*ct + rl*ul*unl + pr*ct + rr*ur*unr
      Fc(3) = pl*st + rl*vl*unl + pr*st + rr*vr*unr
      Fc(4) = rl*hl*unl         + rr*hr*unr

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
      vna =-ua*st + va*ct

!     Eigenvalues
      l1 = abs(una - aa)
      l2 = abs(una)
      l3 = l2
      l4 = abs(una + aa)

!     Eigenvalues with entropy fix
!     e1 = abs(una - aa)
!     l2 = abs(una)
!     l3 = l2
!     e4 = abs(una + aa)

!     del= ETOL*aa
!     if(e1 .lt. del)then
!        l1 = 0.5*(del + e1**2/del)
!     else
!        l1 = e1
!     endif

!     if(e4 .lt. del)then
!        l4 = 0.5*(del + e4**2/del)
!     else
!        l4 = e4
!     endif

!     Difference of conserved variables
      du1 = rr           - rl
      du2 = rr*ur        - rl*ul
      du3 = rr*vr        - rl*vl
      du4 = (rr*hr - pr) - (rl*hl - pl)

!     Amplitudes
      m1 = (ct*du2 + st*du3 - una*du1)/aa
      m2 = (GAMMA-1.0)*(du4 - ua*du2 - va*du3 + qa2*du1)/aa**2

      a4 = 0.5*(m1 + m2)
      a1 = 0.5*(m2 - m1)
      a3 = du1 - a1 - a4
      a2 = ( st*du2 - ct*du3 + vna*du1 )/aa

!     Diffusive flux
      a1l1  = a1*l1
      a2l2  = a2*l2
      a3l3  = a3*l3
      a4l4  = a4*l4
      aact  = aa*ct
      aast  = aa*st

      Fd(1) = a1l1               +               a3l3           + a4l4
      Fd(2) = a1l1*(ua - aact)   + a2l2*aa*st  + a3l3*ua        + &
              a4l4*(ua + aact)
      Fd(3) = a1l1*(va - aast)   - a2l2*aa*ct  + a3l3*va        + &
              a4l4*(va + aast)
      Fd(4) = a1l1*(ha - una*aa) + a2l2*aa*vna + a3l3*0.5*qa2 +   &
              a4l4*(ha + una*aa)

!     Total flux
      flux    =  0.5*( Fc - Fd )
      dflux   = -0.5*Fd

end subroutine roe_flux
