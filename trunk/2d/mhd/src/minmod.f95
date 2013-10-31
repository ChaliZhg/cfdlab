real function minmod(a, b, c)
   implicit none

   real :: a, b, c

   if(a*b > 0.0 .and. b*c > 0.0)then
      minmod = sign(1.0, a) * min(abs(a), abs(b), abs(c))
   else
      minmod = 0.0
   endif

end function minmod

real function minmod2(a, b, c)
   implicit none

   real :: a, b, c

   if(a*c > 0.0)then
      minmod2 = b
   else
      minmod2 = 0.0
   endif

end function minmod2
