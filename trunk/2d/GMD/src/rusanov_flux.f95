subroutine rusanov_flux(sx, sy, conl, conr, flux)
   use comvar
   implicit none

   real :: sx, sy, conl(4), conr(4), flux(4)

   real :: rhol, vexl, veyl, prel, sonl, hl, vnl, laml
   real :: rhor, vexr, veyr, prer, sonr, hr, vnr, lamr
   real :: lam

   ! left state
   rhol = conl(1)
   vexl = conl(2)/conl(1)
   veyl = conl(3)/conl(1)
   prel = (gamma-1.0)*(conl(4) - 0.5*(conl(2)**2 + conl(3)**2)/conl(1))
   sonl = sqrt(gamma*prel/rhol)
   hl   = sonl**2/(gamma-1.0) + 0.5*(vexl**2 + veyl**2)
   vnl  = vexl*sx + veyl*sy
   laml = max( abs(vnl - sonl), abs(vnl), abs(vnl + sonl) )

   ! right state
   rhor = conr(1)
   vexr = conr(2)/conr(1)
   veyr = conr(3)/conr(1)
   prer = (gamma-1.0)*(conr(4) - 0.5*(conr(2)**2 + conr(3)**2)/conr(1))
   sonr = sqrt(gamma*prer/rhor)
   hr   = sonr**2/(gamma-1.0) + 0.5*(vexr**2 + veyr**2)
   vnr  = vexr*sx + veyr*sy
   lamr = max( abs(vnr - sonr), abs(vnr), abs(vnr + sonr) )

   lam  = max(laml, lamr)

   flux(1) = 0.5*(rhol*vnl + rhor*vnr - lam*(conr(1)-conl(1)))
   flux(2) = 0.5*(prel*sx + rhol*vexl*vnl + prer*sx + rhor*vexr*vnr - &
             lam*(conr(2)-conl(2)))
   flux(3) = 0.5*(prel*sy + rhol*veyl*vnl + prer*sy + rhor*veyr*vnr - &
             lam*(conr(3)-conl(3)))
   flux(4) = 0.5*(rhol*vnl*hl + rhor*vnr*hr - lam*(conr(4)-conl(4)))

end subroutine rusanov_flux
