C Calculates finite volume residual
      subroutine fvm_residual(ptype, elem, edge, spts, fpts, opts, bpts,
     &                        coord, elarea, cvarea, ds, dsb, dt, divf,
     &                        prim, prim_old, qx, qy, drmin, mu, mul, 
     &                        cvareamc, wd)
      implicit none
      include 'param.h'

      integer          ptype(npmax), elem(nvemax,ntmax), edge(2,nemax),
     &                 spts(nspmax), fpts(nfpmax), opts(nopmax),
     &                 bpts(nbpmax)
      double precision coord(2, npmax), elarea(ntmax), 
     &                 cvarea(npmax), ds(2,nemax), dsb(2,npmax), 
     &                 dt(npmax), divf(nvar,npmax), prim(nvar,npmax), 
     &                 prim_old(nvar,npmax), qx(nvar,npmax), 
     &                 qy(nvar,npmax), drmin(npmax), mu(npmax), 
     &                 mul(npmax), cvareamc(npmax), wd(npmax)

c     Sutherland viscosity
      if(iflow .ne. inviscid)then
         call sutherland(prim, mul)
         call viscosity(prim, mul, mu)
      endif

c     Find gradients for reconstruction
      call green_gauss(elem, coord, elarea, cvareamc, prim, divf, 
     &                 mul, mu, qx, qy)
      
c     Calculate the flux divergence
      call divergence(edge, spts, fpts, opts, ds, dsb, prim, 
     &                coord, qx, qy, divf)

      return
      end
