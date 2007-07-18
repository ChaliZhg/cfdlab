C Compute gradients for reconstruction and viscous terms in NS equation
      subroutine green_gauss(elem, coord, elarea, cvareamc, prim, 
     &                             divf, mul, mu, qx, qy)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), prim(nvar,npmax), 
     &                 divf(nvar,npmax), mul(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax), elarea(ntmax),
     &                 cvareamc(npmax), mu(npmax)

      integer          is, iv

      do is=1,np
         do iv=1,nvar-1
            divf(iv,is) = 0.0d0
         enddo
      enddo

      if(iflow .eq. inviscid)then
         call green_gauss_euler(elem, coord, cvareamc, prim, 
     &                          divf, qx, qy)
      else
         call green_gauss_visc(elem, coord, elarea, cvareamc, prim, 
     &                          divf, mul, mu, qx, qy)
      endif

      return
      end

C Inviscid flow, compute only derivatives
      subroutine green_gauss_euler(elem, coord, cvareamc, prim, 
     &                             divf, qx, qy)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), prim(nvar,npmax), 
     &                 divf(nvar,npmax), qx(nvar,npmax), qy(nvar,npmax),
     &                 cvareamc(npmax)

      integer          is, iv, jt, k, n1, n2, n3
      double precision dbxx(3), dbyy(3), dxt(nvar), dyt(nvar), ais, n1b6

      n1b6 = 1.0d0/6.0d0

      do is=1,np
         do iv=1,nvar-1
            qx(iv,is)   = 0.0d0
            qy(iv,is)   = 0.0d0
         enddo
      enddo

c     loop on global list of triangles
      do jt=1,nt
         n1          = elem(1,jt)
         n2          = elem(2,jt)
         n3          = elem(3,jt)

         dbxx(1)     = coord(2,n2) - coord(2,n3)
         dbxx(2)     = coord(2,n3) - coord(2,n1)
         dbxx(3)     = coord(2,n1) - coord(2,n2)
         dbyy(1)     = coord(1,n3) - coord(1,n2)
         dbyy(2)     = coord(1,n1) - coord(1,n3)
         dbyy(3)     = coord(1,n2) - coord(1,n1)

         do k=1,nvar-1
            dxt(k)   = prim(k,n1)*dbxx(1) +
     &                 prim(k,n2)*dbxx(2) +
     &                 prim(k,n3)*dbxx(3)
            dyt(k)   = prim(k,n1)*dbyy(1) +
     &                 prim(k,n2)*dbyy(2) +
     &                 prim(k,n3)*dbyy(3)
 
            qx(k,n1) = qx(k,n1) + dxt(k)
            qx(k,n2) = qx(k,n2) + dxt(k)
            qx(k,n3) = qx(k,n3) + dxt(k)

            qy(k,n1) = qy(k,n1) + dyt(k)
            qy(k,n2) = qy(k,n2) + dyt(k)
            qy(k,n3) = qy(k,n3) + dyt(k)
         enddo

      enddo

      do is=1,np
         ais         = n1b6/cvareamc(is)
         do k=1,nvar-1
            qx(k,is) = qx(k,is)*ais
            qy(k,is) = qy(k,is)*ais
         enddo
      enddo

      return
      end

C Version for viscous flow, computes viscous terms also
      subroutine green_gauss_visc(elem, coord, elarea, cvareamc, prim, 
     &                            divf, mul, mu, qx, qy)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), prim(nvar,npmax), 
     &                 divf(nvar,npmax), mul(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax), elarea(ntmax),
     &                 cvareamc(npmax), mu(npmax)

      integer          is, iv, jt, k, n1, n2, n3
      double precision en(3), dbxx(3), dbyy(3), dxt(nvar), dyt(nvar), 
     &                 efx, efy, txx, txy, tyy, gpr, um, vm, ext, eyt, 
     &                 ais, areat, mut, prim5t, fv1, chi, chi3, mult, 
     &                 rhot, n1b6, n2b3

      n1b6 = 1.0d0/6.0d0
      n2b3 = 2.0d0/3.0d0

      do is=1,np
         do iv=1,nvar-1
            qx(iv,is)   = 0.0d0
            qy(iv,is)   = 0.0d0
         enddo
      enddo

c     loop on global list of triangles
      do jt=1,nt
         n1          = elem(1,jt)
         n2          = elem(2,jt)
         n3          = elem(3,jt)

         dbxx(1)     = coord(2,n2) - coord(2,n3)
         dbxx(2)     = coord(2,n3) - coord(2,n1)
         dbxx(3)     = coord(2,n1) - coord(2,n2)
         dbyy(1)     = coord(1,n3) - coord(1,n2)
         dbyy(2)     = coord(1,n1) - coord(1,n3)
         dbyy(3)     = coord(1,n2) - coord(1,n1)

         do k=1,nvar-1
            dxt(k)   = prim(k,n1)*dbxx(1) +
     &                 prim(k,n2)*dbxx(2) +
     &                 prim(k,n3)*dbxx(3)
            dyt(k)   = prim(k,n1)*dbyy(1) +
     &                 prim(k,n2)*dbyy(2) +
     &                 prim(k,n3)*dbyy(3)
 
            qx(k,n1) = qx(k,n1) + dxt(k)
            qx(k,n2) = qx(k,n2) + dxt(k)
            qx(k,n3) = qx(k,n3) + dxt(k)

            qy(k,n1) = qy(k,n1) + dyt(k)
            qy(k,n2) = qy(k,n2) + dyt(k)
            qy(k,n3) = qy(k,n3) + dyt(k)
         enddo

C Viscous flux contribution
C Gradient of internal energy
         en(1) = prim(4,n1)/prim(1,n1)
         en(2) = prim(4,n2)/prim(1,n2)
         en(3) = prim(4,n3)/prim(1,n3)
         ext   = en(1)*dbxx(1) + en(2)*dbxx(2) +  en(3)*dbxx(3)
         eyt   = en(1)*dbyy(1) + en(2)*dbyy(2) +  en(3)*dbyy(3)

c Average velocity on triangle
         um  = (prim(2,n1) + prim(2,n2) + prim(2,n3))/3.0d0
         vm  = (prim(3,n1) + prim(3,n2) + prim(3,n3))/3.0d0

c Average Reynolds number on triangle
         mut    = (mu(n1)     + mu(n2)     + mu(n3))/3.0d0
         mult   = (mul(n1)    + mul(n2)    + mul(n3))/3.0d0
         prim5t = (prim(5,n1) + prim(5,n2) + prim(5,n3))/3.0d0
         rhot   = (prim(1,n1) + prim(1,n2) + prim(1,n3))/3.0d0

         chi    = prim5t*rhot/mult
         chi3   = chi**3
         fv1    = chi3/(chi3 + Cv11)

c Viscous stresses
         txx = n2b3*mut*( 2.0d0*dxt(2) - dyt(3) )
         txy =      mut*( dyt(2) + dxt(3) )
         tyy = n2b3*mut*( 2.0d0*dyt(3) - dxt(2) )

C Viscous Energy flux
         gpr = gamma*(mult/prandtl + 
     &                rhot*prim5t*fv1/prandtl_turb)/gamma1
         efx = txx*um + txy*vm + gpr*ext
         efy = txy*um + tyy*vm + gpr*eyt

c Add the viscous flux
         areat      = 0.25d0/elarea(jt)
         divf(2,n1) = divf(2,n1) + (txx*dbxx(1) + txy*dbyy(1))*areat
         divf(3,n1) = divf(3,n1) + (txy*dbxx(1) + tyy*dbyy(1))*areat
         divf(4,n1) = divf(4,n1) + (efx*dbxx(1) + efy*dbyy(1))*areat

         divf(2,n2) = divf(2,n2) + (txx*dbxx(2) + txy*dbyy(2))*areat
         divf(3,n2) = divf(3,n2) + (txy*dbxx(2) + tyy*dbyy(2))*areat
         divf(4,n2) = divf(4,n2) + (efx*dbxx(2) + efy*dbyy(2))*areat

         divf(2,n3) = divf(2,n3) + (txx*dbxx(3) + txy*dbyy(3))*areat
         divf(3,n3) = divf(3,n3) + (txy*dbxx(3) + tyy*dbyy(3))*areat
         divf(4,n3) = divf(4,n3) + (efx*dbxx(3) + efy*dbyy(3))*areat

      enddo   

      do is=1,np
         ais         = n1b6/cvareamc(is)
         do k=1,nvar-1
            qx(k,is) = qx(k,is)*ais
            qy(k,is) = qy(k,is)*ais
         enddo
      enddo

      return
      end
