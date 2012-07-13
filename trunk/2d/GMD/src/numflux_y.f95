subroutine numflux_y(conjm1, conj, conjp1, conjp2, flux, dflux)
   use comvar
   implicit none

   real :: conjm1(4), conj(4), conjp1(4), conjp2(4), flux(4), dflux(4)

   real :: conl(4), conr(4)
   real :: rhol, vexl, veyl, prel, sonl, hl, laml
   real :: rhor, vexr, veyr, prer, sonr, hr, lamr
   real :: lam

   ! reconstructed states
   call reconstruct(conjm1, conj, conjp1, conjp2, conl, conr)

   if(fluxtype == iroe)then
      call roe_flux(0.0, 1.0, conl, conr, flux, dflux)
   else if(fluxtype == irusanov)then
      call rusanov_flux(0.0, 1.0, conl, conr, flux, dflux)
   else if(fluxtype == iadv)then
      call adv_flux(0.0, 1.0, conl, conr, flux, dflux)
   else if(fluxtype == icusp)then
      call cusp_flux(0.0, 1.0, conl, conr, flux, dflux)
   else if(fluxtype == ikep)then
      call kep_flux(0.0, 1.0, conl, conr, flux, dflux)
   else if(fluxtype == ikepes)then
      call kepes_flux(0.0, 1.0, conl, conr, flux, dflux)
   else
      write(*,*)'Uknown flux type fluxtype =', fluxtype
      stop
   endif

   if(ikepes_diss == yes)then
      call kepes_diss(0.0, 1.0, conjm1, conj, conjp1, conjp2, dflux)
      flux = flux - dflux
   endif

end subroutine numflux_y
