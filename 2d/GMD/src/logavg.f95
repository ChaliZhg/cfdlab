real function logavg(a, b)
implicit none
real :: a, b

real :: xi, f, u, u2, u3, FF

   xi = b/a
   f = (xi - 1.0) / (xi + 1.0)
   u = f * f

   if (u < 1.0e-2)then
      u2 = u * u
      u3 = u2 * u
      FF = 1.0 + u/3.0 + u2/5.0 + u3/7.0
   else
      FF = log(xi)/2.0/f
   endif

   logavg = 0.5*(a+b)/FF

end function logavg
