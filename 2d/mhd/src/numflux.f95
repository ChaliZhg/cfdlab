subroutine numflux(lx, ly, prijm1, prij, prijp1, prijp2, flux)
   use comvar
   implicit none

   real :: lx, ly
   real :: prijm1(nvar), prij(nvar), prijp1(nvar), prijp2(nvar), &
           flux(nvar)

   real :: pril(nvar), prir(nvar)

   ! reconstructed states
   call reconstruct(prijm1, prij, prijp1, prijp2, pril, prir)

   if(fluxtype == ient)then
      call ent_flux(lx, ly, prij, prijp1, pril, prir, flux)
   else
      write(*,*)'Uknown flux type fluxtype =', fluxtype
      stop
   endif

end subroutine numflux
