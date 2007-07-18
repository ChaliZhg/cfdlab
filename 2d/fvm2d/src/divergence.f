C Calculates the flux divergence, which is also the steady state residual
      subroutine divergence(edge, spts, fpts, opts, ds, dsb, 
     &                      prim, coord, qx, qy, divf)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), spts(nspmax), fpts(nfpmax),
     &                 opts(nopmax)
      double precision ds(2,nemax), dsb(2,npmax), prim(nvar,npmax), 
     &                 coord(2,npmax), qx(nvar,npmax), qy(nvar,npmax), 
     &                 divf(nvar,npmax)

      integer          i, ip, iv
      double precision flux(nvar)

      if(iflux .eq. iroe)then
         call divergence_roe(edge, ds, prim, coord, qx, qy, divf)
      elseif(iflux .eq. ikfvs)then
         call divergence_kfvs(edge, ds, prim, coord, qx, qy, divf)
      elseif(iflux .eq. ihllc)then
         call divergence_hllc(edge, ds, prim, coord, qx, qy, divf)
      endif

C     Flux for solid boundary edges
      do i=1,nsp
         ip = spts(i)
         call solid_flux(ip, prim, dsb, flux)
         divf(2,ip) = divf(2,ip) + flux(2)
         divf(3,ip) = divf(3,ip) + flux(3)
      enddo

C     Flux for far field points - far field bc
      do i=1,nfp
         ip = fpts(i)
         call farfield_flux(ip, coord, prim, dsb, flux)
         do iv=1,nvar-1
            divf(iv,ip) = divf(iv,ip) + flux(iv)
         enddo
      enddo

C     Flux for outflow points
      do i=1,nop
         ip = opts(i)
         call outflow_flux(ip, prim, dsb, flux)
         do iv=1,nvar-1
            divf(iv,ip) = divf(iv,ip) + flux(iv)
         enddo
      enddo

      return
      end
