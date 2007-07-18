C.....Spalart-Allmaras turbulence model
C.....Implicit treatment of destruction term
      subroutine sa_model(elem, edge, spts, fpts, opts, coord, elarea,
     &                    cvarea, prim, prim_old, wd, mul, ds, 
     &                    dsb, dt, qx, qy)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax), edge(2,nemax), spts(nspmax),
     &                 fpts(nfpmax), opts(nopmax)
      double precision coord(2,npmax), elarea(ntmax), cvarea(npmax),
     &                 prim(nvar,npmax),
     &                 wd(npmax), mul(npmax), ds(2,nemax),
     &                 dsb(2,npmax), prim_old(nvar,npmax), dt(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax)

      double precision divf(npmax)


      call viscous_sa(elem, coord, elarea, prim, divf, mul)
      
      call divergence_sa(edge, fpts, opts, ds, dsb, prim, divf)

      call update_sa(divf, prim_old, prim, cvarea, spts, dt, mul, wd, 
     &               qx, qy)


      return
      end

c Viscous for SA model
      subroutine viscous_sa(elem, coord, elarea, prim, divf, mul)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), prim(nvar,npmax), 
     &                 divf(npmax), mul(npmax), elarea(ntmax)

      integer          is, jt, n1, n2, n3
      double precision dbxx(3), dbyy(3), nux, nuy, areat, prim5t,
     &                 divsa1, divsa2, divsa3, mult, mutb, rhot

      do is=1,np
            divf(is) = 0.0d0
      enddo

      do jt=1,nt
         areat    = 0.25d0/elarea(jt)

         n1       = elem(1,jt)
         n2       = elem(2,jt)
         n3       = elem(3,jt)

         dbxx(1)  = coord(2,n2) - coord(2,n3)
         dbxx(2)  = coord(2,n3) - coord(2,n1)
         dbxx(3)  = coord(2,n1) - coord(2,n2)
         dbyy(1)  = coord(1,n3) - coord(1,n2)
         dbyy(2)  = coord(1,n1) - coord(1,n3)
         dbyy(3)  = coord(1,n2) - coord(1,n1)

         nux      = prim(5,n1)*dbxx(1) + 
     &              prim(5,n2)*dbxx(2) +
     &              prim(5,n3)*dbxx(3)
         nuy      = prim(5,n1)*dbyy(1) +
     &              prim(5,n2)*dbyy(2) +
     &              prim(5,n3)*dbyy(3)

         mult     = (mul(n1)    + mul(n2)    + mul(n3))/3.0d0
         prim5t   = (prim(5,n1) + prim(5,n2) + prim(5,n3))/3.0d0
         rhot     = (prim(1,n1) + prim(1,n2) + prim(1,n3))/3.0d0
         mutb     = mult/rhot + prim5t

         divsa1   = (Cb2Sig1*mutb - Cb2Sig2*(mul(n1)/prim(1,n1) +
     &              prim(5,n1)))*(nux*dbxx(1) + nuy*dbyy(1))*areat

         divsa2   = (Cb2Sig1*mutb - Cb2Sig2*(mul(n2)/prim(1,n2) +
     &              prim(5,n2)))*(nux*dbxx(2) + nuy*dbyy(2))*areat

         divsa3   = (Cb2Sig1*mutb - Cb2Sig2*(mul(n3)/prim(1,n3) +
     &              prim(5,n3)))*(nux*dbxx(3) + nuy*dbyy(3))*areat

         divf(n1) = divf(n1) + divsa1
         divf(n2) = divf(n2) + divsa2
         divf(n3) = divf(n3) + divsa3
     
      enddo

      return
      end

c Convective flux for SA model
      subroutine divergence_sa(edge, fpts, opts, ds, dsb, prim, divf)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), fpts(nfpmax), opts(nopmax)
      double precision ds(2,nemax), dsb(2,npmax), prim(nvar,npmax), 
     &                 divf(npmax)

      integer          i, n1, n2
      double precision unl, unr, una, flux

      do i=1,ne
            n1   = edge(1,i)
            n2   = edge(2,i)

            unl  = prim(2,n1)*ds(1,i) + prim(3,n1)*ds(2,i)
            unr  = prim(2,n2)*ds(1,i) + prim(3,n2)*ds(2,i)
            una  = 0.5d0*(unl + unr)

            if(una .gt. 0.0d0)then
                  flux = una*prim(5,n1)
            else
                  flux = una*prim(5,n2)
            endif

            divf(n1) = divf(n1) + flux
            divf(n2) = divf(n2) - flux
      enddo

      call outer_flux_sa(fpts, opts, prim, dsb, divf)

      return
      end

C     Flux for outer boundary edges - far field bc
      subroutine outer_flux_sa(fpts, opts, prim, dsb, divf)
      implicit none
      include 'param.h'
      integer          fpts(nfpmax), opts(nopmax)
      double precision prim(nvar,npmax), dsb(2,npmax), divf(npmax)

      integer          i, ip
      double precision un, flux

      do i=1,nfp
            ip = fpts(i)
            un = prim(2,ip)*dsb(1,ip) + prim(3,ip)*dsb(2,ip)
            if(un .gt. 0.0d0)then
                  flux = prim(5,ip)*un
            else
                  flux = prim_inf(5)*un
            endif

            divf(ip) = divf(ip) + flux
      enddo

      do i=1,nop
            ip = opts(i)
            un = prim(2,ip)*dsb(1,ip) + prim(3,ip)*dsb(2,ip)
            if(un .gt. 0.0d0)then
                  flux = prim(5,ip)*un
            else
                  flux = prim_inf(5)*un
            endif

            divf(ip) = divf(ip) + flux
      enddo

      return
      end

C.....Update the turbulent viscosity
      subroutine update_sa(divf, prim_old, prim, cvarea, spts, dt,
     &                     mul, wd, qx, qy)
      implicit none
      include 'param.h'
      integer          spts(nspmax)
      double precision divf(npmax), prim_old(nvar,npmax), 
     &                 prim(nvar,npmax), cvarea(npmax), dt(npmax),
     &                 mul(npmax), wd(npmax), qx(nvar,npmax),
     &                 qy(nvar,npmax)

      integer          i, j
      double precision chi, chi3, fv1, fv2, fv3, r, Omega, S, g, fw, 
     &                 fact, n1b6, source

      n1b6 = 1.0d0/6.0d0
      do i=1,np
            chi    = prim(5,i)*prim(1,i)/mul(i)
            chi3   = chi**3
            fv1    = chi3/(chi3 + Cv11)
            fv2    = 1.0d0/(1.0d0 + chi/Cv2)**3
            r      = prim(5,i)/wd(i)/(wd(i)*kolm2)
            Omega  = qy(2,i) - qx(3,i)
            fv3    = (1.0d0 + chi*fv1)*(1.0d0 - fv2)/dmax1(chi, 0.001d0)
            S      = fv3*dabs(Omega) + r*fv2
            r      = r/S
            r      = dmin1(r,2.0d0)
            g      = r + Cw2*(r**6 - r)
            fw     = g*(Cw31/(g**6 + Cw32))**n1b6
            source = Cb1*S*prim(5,i)
            fact   = 1.0d0 + Cw1*fw*(prim(5,i)/wd(i))*(dt(i)/wd(i))
            prim(5,i) = (prim(5,i) -  (dt(i)/cvarea(i))*divf(i) +
     &                   dt(i)*source)/fact
      enddo
                  
      do i=1,nsp
            j         = spts(i)
            prim(5,j) = 0.0d0
      enddo

      return
      end
