subroutine reconstruct(conjm1, conj, conjp1, conjp2, conl, conr)
   use comvar
   implicit none

   real    :: conjm1(nvar), conj(nvar), conjp1(nvar), conjp2(nvar)
   real    :: conl(nvar), conr(nvar)

   integer :: i
   real    :: kkk=1.0/3.0
   real    :: minmod

   ! reconstructed states
   if(limtype == ford)then
   ! first order
      do i=1,nvar
         conl(i) = conj(i)
         conr(i) = conjp1(i)
      enddo
   else if(limtype == muscl3)then
   !muscl scheme
      do i=1,nvar
         conl(i) = conj(i)   + 0.25*( (1.0-kkk)*(conj(i) - conjm1(i)) &
                                    + (1.0+kkk)*(conjp1(i) - conj(i)) )
         conr(i) = conjp1(i) - 0.25*( (1.0+kkk)*(conjp1(i) - conj(i)) &
                                    + (1.0-kkk)*(conjp2(i) - conjp1(i)) )
      enddo
   else if(limtype == mmod)then
   ! minmod limiter
      do i=1,nvar
         conl(i) = conj(i) + 0.5*minmod( conj(i)-conjm1(i), &
                                          0.5*(conjp1(i)-conjm1(i)), &
                                          conjp1(i)-conj(i) )
         conr(i) = conjp1(i) - 0.5*minmod( conjp1(i)-conj(i), &
                                          0.5*(conjp2(i)-conj(i)), &
                                          conjp2(i)-conjp1(i) )
      enddo
   endif

end subroutine reconstruct
