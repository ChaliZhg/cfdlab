subroutine divflux(lx, ly, ll, left, right, rr, flux_l, flux_r)
      use comvar
      implicit none
      real :: lx, ly
      real,dimension(nvar),intent(in) :: ll, left, right, rr
      real,dimension(nvar),intent(out) :: flux_l, flux_r

      real :: un, phi_ll, phi_l, phi_r, phi_rr
      real :: uB_ll, uB_l, uB_r, uB_rr
      real :: minmod, dl, dr, dphi

      un = 0.5 * ( (left(2) + right(2))*lx + (left(3)+right(3))*ly )

      call compute_phi(ll,    phi_ll, uB_ll)
      call compute_phi(left,  phi_l , uB_l )
      call compute_phi(right, phi_r , uB_r )
      call compute_phi(rr,    phi_rr, uB_rr)

      dl = minmod(phi_l-phi_ll, 0.5*(phi_r-phi_ll), phi_r-phi_l)
      dr = minmod(phi_r-phi_l,  0.5*(phi_rr-phi_l), phi_rr-phi_r)

      dphi = (phi_r - 0.5*dr) - (phi_l + 0.5*dl)
      dphi = -0.125*abs(un)*dphi

      flux_l(1)   = 0.0
      flux_l(2:4) = left(6:8) * dphi
      flux_l(5)   = uB_l * dphi
      flux_l(6:8) = left(2:4) * dphi

      flux_r(1)   = 0.0
      flux_r(2:4) = right(6:8) * dphi
      flux_r(5)   = uB_r * dphi
      flux_r(6:8) = right(2:4) * dphi

end subroutine divflux

subroutine compute_phi(prim, phi, uB)
      use comvar
      implicit none
      real :: prim(nvar), phi, uB

      uB =  prim(2)*prim(6) + prim(3)*prim(7) + prim(4)*prim(8)
      phi = prim(1)/prim(5) * uB

end subroutine compute_phi
