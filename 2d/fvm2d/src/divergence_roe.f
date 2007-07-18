C Calculates the flux divergence, which is also the steady state residual
      subroutine divergence_roe(edge, ds, prim, coord, qx, qy, divf)
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
            call roe(ds(1,ie), ds(2,ie), priml, primr, flux)
            do iv=1,nvar-1
                  divf(iv,e1) = divf(iv,e1) + flux(iv)
                  divf(iv,e2) = divf(iv,e2) - flux(iv)
            enddo
      enddo

      return
      end

c Roe flux function
      subroutine roe(nx, ny, priml, primr, flux)
      implicit none
      include 'param.h'
      double precision nx, ny, priml(nvar), primr(nvar), flux(nvar)

      integer  i
      double precision rl, ul, vl, pl, al2, hl, rr, ur, vr, pr, ar2, hr,
     &                 ua, va, qa2, aa2, aa, ha,
     &                 ql2, qr2, rl12, rr12, rd,
     &                 unl, una, vna,
     &                 length, ct, st, Fc(4), Fd(4),
     &                 m1, m2, a1, a2, a3, a4, l1, l2, l3, l4,
     &                 a1l1, a2l2, a3l3, a4l4, aact, aast,
     &                 du1, du2, du3, du4

      length= dsqrt(nx**2 + ny**2)
      ct    = nx/length
      st    = ny/length

C     Left state
      rl = priml(1)
      ul = priml(2)
      vl = priml(3)
      pl = priml(4)
      al2= GAMMA*pl/rl
      ql2= ul**2 + vl**2
      hl = al2/GAMMA1 + 0.5d0*ql2

C     Right state
      rr = primr(1)
      ur = primr(2)
      vr = primr(3)
      pr = primr(4)
      ar2= GAMMA*pr/rr
      qr2= ur**2 + vr**2
      hr = ar2/GAMMA1 + 0.5d0*qr2

C     Rotated velocity
      unl = ul*ct + vl*st

C     Centered flux
      Fc(1) = rl*unl
      Fc(2) = pl*ct + rl*ul*unl
      Fc(3) = pl*st + rl*vl*unl
      Fc(4) = rl*hl*unl

C     Roe average
      rl12 = dsqrt(rl)
      rr12 = dsqrt(rr)
      rd   = 1.0d0/(rl12 + rr12)

      ua   = (ul*rl12 + ur*rr12)*rd
      va   = (vl*rl12 + vr*rr12)*rd
      ha   = (hl*rl12 + hr*rr12)*rd
      qa2  = ua**2 + va**2
      aa2  = GAMMA1*(ha - 0.5d0*qa2)

      if(aa2 .gt. 0.0d0)then
            aa = dsqrt(aa2)
      else
            print*,'Sonic speed is negative'
            stop
      endif

      una = ua*ct + va*st
      vna =-ua*st + va*ct

C     Eigenvalues with entropy fix
      l1 = dmin1( una - aa, 0.0d0)
      l2 = dmin1( una,      0.0d0)
      l3 = l2
      l4 = dmin1( una + aa, 0.0d0)

c     Difference of conserved variables
      du1 = rr    - rl
      du2 = rr*ur - rl*ul
      du3 = rr*vr - rl*vl
      du4 = (pr - pl)/GAMMA1 + 0.5d0*(rr*qr2  - rl*ql2)

c     Amplitudes
      m1 = (ct*du2 + st*du3 - una*du1)/aa
      m2 = GAMMA1*(du4 - ua*du2 - va*du3 + qa2*du1)/aa**2

      a4 = 0.5d0*(m1 + m2)
      a1 = 0.5d0*(m2 - m1)
      a3 = du1 - a1 - a4
      a2 = ( st*du2 - ct*du3 + vna*du1 )/aa

c     Diffusive flux
      a1l1  = a1*l1
      a2l2  = a2*l2
      a3l3  = a3*l3
      a4l4  = a4*l4
      aact  = aa*ct
      aast  = aa*st

      Fd(1) = a1l1               +               a3l3           + a4l4
      Fd(2) = a1l1*(ua - aact)   + a2l2*aa*st  + a3l3*ua        +
     &        a4l4*(ua + aact)
      Fd(3) = a1l1*(va - aast)   - a2l2*aa*ct  + a3l3*va        +
     &        a4l4*(va + aast)
      Fd(4) = a1l1*(ha - una*aa) + a2l2*aa*vna + a3l3*0.5d0*qa2 +
     &        a4l4*(ha + una*aa)

c     Total flux
      do i=1,4
            flux(i) = ( Fc(i) + Fd(i) )*length
      enddo

      return
      end
