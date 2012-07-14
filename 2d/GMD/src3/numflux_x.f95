subroutine numflux_x(conjm1, conj, conjp1, conjp2, flux)
   use comvar
   implicit none

   real :: conjm1(nvar), conj(nvar), conjp1(nvar), conjp2(nvar), flux(nvar)

   real :: conl(nvar), conr(nvar), dflux(nvar)

   ! reconstructed states
   call reconstruct(conjm1, conj, conjp1, conjp2, conl, conr)

   if(fluxtype == ikep)then
      call kep_flux(1.0, 0.0, 0.0, conl, conr, flux)
   else if(fluxtype == ikepes)then
      call kepes_flux(1.0, 0.0, 0.0, conl, conr, flux)
   else
      write(*,*)'Uknown flux type fluxtype =', fluxtype
      stop
   endif

   if(ikepes_diss == yes)then
      call kepes_diss(1.0, 0.0, 0.0, conjm1, conj, conjp1, conjp2, dflux)
      flux = flux - dflux
   endif

end subroutine numflux_x
