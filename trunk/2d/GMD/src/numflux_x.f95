subroutine numflux_x(conjm1, conj, conjp1, conjp2, flux)
   use comvar
   implicit none

   real :: conjm1(4), conj(4), conjp1(4), conjp2(4), flux(4)

   real :: conl(4), conr(4)
   real :: rhol, vexl, veyl, prel, sonl, hl, laml
   real :: rhor, vexr, veyr, prer, sonr, hr, lamr
   real :: lam

   ! reconstructed states
   call reconstruct(conjm1, conj, conjp1, conjp2, conl, conr)

   if(fluxtype == iroe)then
      call roe_flux(1.0, 0.0, conl, conr, flux)
   else if(fluxtype == irusanov)then
      call rusanov_flux(1.0, 0.0, conl, conr, flux)
   else
      write(*,*)'Uknown flux type fluxtype =', fluxtype
      stop
   endif

end subroutine numflux_x
