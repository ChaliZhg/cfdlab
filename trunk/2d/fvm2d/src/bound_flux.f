C.....Flux for solid boundary edges
      subroutine solid_flux(pt, prim, dsb, flux)
      implicit none
      include 'param.h'
      integer          pt
      double precision prim(nvar,npmax), dsb(2,npmax), flux(nvar)

      flux(2) = prim(4,pt)*dsb(1,pt)
      flux(3) = prim(4,pt)*dsb(2,pt)

      return
      end

C.....Flux for outflow points, like supersonic outflow
      subroutine outflow_flux(pt, prim, dsb, flux)
      implicit none
      include 'param.h'
      integer          pt
      double precision prim(nvar,npmax), dsb(2,npmax), flux(nvar)

      double precision un, E

      un      = prim(2,pt)*dsb(1,pt) + prim(3,pt)*dsb(2,pt)

      flux(1) = prim(1,pt)*un
      flux(2) = prim(4,pt)*dsb(1,pt) + prim(1,pt)*prim(2,pt)*un
      flux(3) = prim(4,pt)*dsb(2,pt) + prim(1,pt)*prim(3,pt)*un
      E       = prim(4,pt)/GAMMA1 + 0.5d0*prim(1,pt)*( prim(2,pt)**2 +
     &          prim(3,pt)**2 )
      flux(4) = (E + prim(4,pt))*un

      return
      end

C.....Find flux for farfield boundary edges using Riemann invariants
C.....This is not appropriate for unsteady flows
C.....Far-field vortex correction, see Blazek
C.....Chord length is assumed to be 1, otherwise change circ
C.....(xref,yref) is center of vortex, taken to be quarter chord point
      subroutine farfield_flux(pt, coord, prim, dsb, flux)
      implicit none
      include 'param.h'
      integer          pt
      double precision coord(2,npmax), prim(nvar,npmax), dsb(2,npmax), 
     &                 flux(nvar)

      double precision dr, nx, ny, un, ut, a, uninf, utinf, ainf, 
     &                 l1, l2, S, S_int, S_inf, R1_int, R2_int,
     &                 R1_inf, R2_inf, R1, R2, r, u, v, p, E,
     &                 dx, dy, dref, theta, circ, fact1,
     &                 fact2, fact3, fact, uinf, vinf, qinf, 
     &                 pinf1, pinf2, pinf, rinf

      if(mach_inf .lt. 1.0d0 .and. vortex .eq. yes)then
            dx    = coord(1,pt) - xref
            dy    = coord(2,pt) - yref
            dref  = dsqrt(dx**2 + dy**2)
            theta = datan2(dy, dx)
            circ  = 0.5d0*q_inf*Cl
            fact1 = circ*dsqrt(1.0d0 - mach_inf**2)
            fact2 = 2.0d0*M_PI*dref
            fact3 = 1.0d0 - (mach_inf*dsin(theta - aoa))**2
            fact  = fact1/(fact2*fact3)
            uinf  = u_inf + fact*dsin(theta)
            vinf  = v_inf - fact*dcos(theta)
            qinf  = dsqrt(uinf**2 + vinf**2)
            pinf1 = p_inf**(GAMMA1/GAMMA)
            pinf2 = 0.5d0*(GAMMA1/GAMMA)*r_inf*(q_inf**2 -
     &              qinf**2)/p_inf**(1.0d0/GAMMA)
            pinf  = (pinf1 + pinf2)**(GAMMA/GAMMA1)
            rinf  = r_inf*(pinf/p_inf)**(1.0d0/GAMMA)
      else
            uinf  = u_inf
            vinf  = v_inf
            pinf  = p_inf
            rinf  = r_inf
      endif

      dr = dsqrt( dsb(1,pt)**2 + dsb(2,pt)**2 )
      nx = dsb(1,pt)/dr
      ny = dsb(2,pt)/dr
      un = prim(2,pt)*nx + prim(3,pt)*ny
      a  = dsqrt(GAMMA*prim(4,pt)/prim(1,pt))
      S_int = prim(4,pt)/prim(1,pt)**GAMMA

      uninf = uinf*nx + vinf*ny
      utinf =-uinf*ny + vinf*nx
      ainf  = dsqrt(GAMMA*pinf/rinf)
      S_inf = pinf/rinf**GAMMA

      l1 = un - a
      l2 = un + a

      R1_int = 0.5d0*un - a/GAMMA1
      R2_int = 0.5d0*un + a/GAMMA1

      R1_inf = 0.5d0*uninf - ainf/GAMMA1
      R2_inf = 0.5d0*uninf + ainf/GAMMA1

      if( l1 .ge. 0.0d0)then
            R1 = R1_int
      else
            R1 = R1_inf
      endif

      if( l2 .ge. 0.0d0)then
            R2 = R2_int
      else
            R2 = R2_inf
      endif

      un = R1 + R2
      a  = 0.5d0*(R2 - R1)*GAMMA1

      if( un .ge. 0.0d0)then
            S  = S_int
            ut =-prim(2,pt)*ny + prim(3,pt)*nx
      else
            S  = S_inf
            ut = utinf
      endif

      r       = (a**2/(GAMMA*S))**(1.0d0/GAMMA1)
      p       = S*r**GAMMA
      u       = un*nx - ut*ny
      v       = un*ny + ut*nx
      E       = p/GAMMA1 + 0.5d0*r*(u**2 + v**2)

      flux(1) = dr*r*un
      flux(2) = dr*(p*nx + r*u*un)
      flux(3) = dr*(p*ny + r*v*un)
      flux(4) = dr*(E + p)*un

      return
      end
