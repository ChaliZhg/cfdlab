subroutine numflux(lx, ly, prijm1, prij, prijp1, prijp2, flux)
   use comvar
   implicit none

   real :: lx, ly
   real :: prijm1(nvar), prij(nvar), prijp1(nvar), prijp2(nvar), &
           flux(nvar)

   real :: pril(nvar), prir(nvar)

   if(fluxtype == ient)then ! muscl scheme
      ! reconstructed states
      call reconstruct(prijm1, prij, prijp1, prijp2, pril, prir)
      call ent_flux(lx, ly, pril, prir, pril, prir, flux)
   elseif(fluxtype == ifent)then ! tecno scheme
      call fent_flux(lx, ly, prijm1, prij, prijp1, prijp2, flux)
   else
      write(*,*)'Uknown flux type fluxtype =', fluxtype
      stop
   endif

end subroutine numflux
