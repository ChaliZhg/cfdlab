subroutine solveGMD(rho, vex, vey, pre, omg, co0, co1, phi, psi, &
                    phid, psid, res)

   use comvar

   implicit none

   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: pre(-1:nx+2, -1:ny+2)
   real :: omg( 1:nx+1,  1:ny+1)
   real :: co0(4, -1:nx+2, -1:ny+2)
   real :: co1(4, -1:nx+2, -1:ny+2)
   real :: phi(4,nx+1,ny+1)
   real :: psi(4,nx+1,ny+1)
   real :: phid(4,nx+1,ny+1)
   real :: psid(4,nx+1,ny+1)
   real :: res(4,0:nx+1,0:ny+1)

   integer :: it, i, j, rks
   real    :: lambda
   real    :: xflux(4), yflux(4)
   real    :: dxflux(4), dyflux(4)
   real    :: time
   real    :: resid(4), resid1(4)
   real    :: etax, etay, neta, kx, ky, omg0, xd, yd, ep
   real    :: ane, ase, anw, asw
   real    :: bne, bse, bnw, bsw


   ! set initial condition
   call init_cond(rho, vex, vey, pre)
   call prim2cons(rho, vex, vey, pre, co1)
   call periodic(co1)
   call cons2prim(co1, rho, vex, vey, pre)
   call saveprim(0.0, rho, vex, vey, pre)
   call vorticity(rho, vex, vey, pre, omg)
   call savevort(0.0, omg)
   call timestep(rho, vex, vey, pre)

   lambda = dt/dx/dy
   time   = 0.0

   do it=1,itmax

      co0(:,:,:) = co1(:,:,:)

      do rks=1,3

         call cons2prim(co1,rho,vex,vey,pre)
         call vorticity(rho, vex, vey, pre, omg)

         ! symmetric potential phi
         phi  = 0.0
         phid = 0.0

         do i=0,nx
            j = 0
            call numflux_x(co1(:,i-1,j), co1(:,i,j), &
                           co1(:,i+1,j), co1(:,i+2,j), xflux, dxflux)
            phi (:,i+1,j+1) = xflux(:)
            phid(:,i+1,j+1) = dxflux(:)
            j = ny + 1
            call numflux_x(co1(:,i-1,j), co1(:,i,j), &
                           co1(:,i+1,j), co1(:,i+2,j), xflux, dxflux)
            phi (:,i+1,j) = xflux(:)
            phid(:,i+1,j) = dxflux(:)
         enddo

         do i=0,nx
            do j=1,ny
               call numflux_x(co1(:,i-1,j), co1(:,i,j), &
                              co1(:,i+1,j), co1(:,i+2,j), xflux, dxflux)
               phi (:,i+1,j)   = phi (:,i+1,j)   + xflux(:)
               phi (:,i+1,j+1) = phi (:,i+1,j+1) + xflux(:)
               phid(:,i+1,j)   = phid(:,i+1,j)   + dxflux(:)
               phid(:,i+1,j+1) = phid(:,i+1,j+1) + dxflux(:)
            enddo
         enddo

         ! symmetric potential psi
         psi  = 0.0
         psid = 0.0

         do j=0,ny
            i = 0
            call numflux_y(co1(:,i,j-1), co1(:,i,j), &
                           co1(:,i,j+1), co1(:,i,j+2), yflux, dyflux)
            psi (:,i+1,j+1) = yflux(:)
            psid(:,i+1,j+1) = dyflux(:)

            i = nx + 1
            call numflux_y(co1(:,i,j-1), co1(:,i,j), &
                           co1(:,i,j+1), co1(:,i,j+2), yflux, dyflux)
            psi (:,i,j+1) = yflux(:)
            psid(:,i,j+1) = dyflux(:)
         enddo

         do i=1,nx
            do j=0,ny
               call numflux_y(co1(:,i,j-1), co1(:,i,j), &
                              co1(:,i,j+1), co1(:,i,j+2), yflux, dyflux)
               psi (:,i,  j+1) = psi (:,i,  j+1) + yflux(:)
               psi (:,i+1,j+1) = psi (:,i+1,j+1) + yflux(:)
               psid(:,i,  j+1) = psid(:,i,  j+1) + dyflux(:)
               psid(:,i+1,j+1) = psid(:,i+1,j+1) + dyflux(:)
            enddo
         enddo

         ! compute residual
         resid = 0.0
         do i=1,nx
            do j=1,ny

               res(:,i,j) = ( ( phi(:,i+1,j) + phi(:,i+1,j+1) ) - &
                              ( phi(:,i  ,j) + phi(:,i  ,j+1) ) )*dy + &
                            ( ( psi(:,i,j+1) + psi(:,i+1,j+1) ) - &
                              ( psi(:,i  ,j) + psi(:,i+1,j  ) ) )*dx
               res(:,i,j) = 0.25*res(:,i,j)

               ! add vorticity confinement terms
               if(vconf==yes)then
                  etax = -0.5*( (abs(omg(i+1,j)) + abs(omg(i+1,j+1))) - &
                                (abs(omg(i  ,j)) + abs(omg(i  ,j+1))) )/dx
                  etay = -0.5*( (abs(omg(i,j+1)) + abs(omg(i+1,j+1))) - &
                                (abs(omg(i  ,j)) + abs(omg(i+1,j  ))) )/dy
                  neta = sqrt(etax**2 + etay**2)
                  if(neta > 0.0)then
                     etax = etax/neta
                     etay = etay/neta
                  else
                     etax = 0.0
                     etay = 0.0
                  endif

                  ! x velocity at vertices
                  ane = 0.25*(vex(i,j) + vex(i+1,j) + vex(i,j+1) + vex(i+1,j+1))
                  ase = 0.25*(vex(i,j) + vex(i+1,j) + vex(i,j-1) + vex(i+1,j-1))
                  anw = 0.25*(vex(i,j) + vex(i-1,j) + vex(i,j+1) + vex(i-1,j+1))
                  asw = 0.25*(vex(i,j) + vex(i-1,j) + vex(i,j-1) + vex(i-1,j-1))

                  ! y velocity at vertices
                  bne = 0.25*(vey(i,j) + vey(i+1,j) + vey(i,j+1) + vey(i+1,j+1))
                  bse = 0.25*(vey(i,j) + vey(i+1,j) + vey(i,j-1) + vey(i+1,j-1))
                  bnw = 0.25*(vey(i,j) + vey(i-1,j) + vey(i,j+1) + vey(i-1,j+1))
                  bsw = 0.25*(vey(i,j) + vey(i-1,j) + vey(i,j-1) + vey(i-1,j-1))

                  ! vorticity at cell center (i,j)
                  omg0= 0.5*((bne + bse) - (bnw + bsw))/dx - &
                        0.5*((anw + ane) - (asw + ase))/dy

                  kx = -etay*omg0
                  ky = +etax*omg0

                  xd =-0.25*( ( psid(3,i+1,j) + psid(3,i+1,j+1) ) - &
                              ( psid(3,i  ,j) + psid(3,i  ,j+1) ) )/dx + &
                       0.25*( ( psid(2,i,j+1) + psid(2,i+1,j+1) ) - &
                              ( psid(2,i  ,j) + psid(2,i+1,j  ) ) )/dy
                  yd = 0.25*( ( phid(3,i+1,j) + phid(3,i+1,j+1) ) - &
                              ( phid(3,i  ,j) + phid(3,i  ,j+1) ) )/dx - &
                       0.25*( ( phid(2,i,j+1) + phid(2,i+1,j+1) ) - &
                              ( phid(2,i  ,j) + phid(2,i+1,j  ) ) )/dy

                  if(kx**2 + ky**2 > 0.0)then
                     ep = (xd*kx + yd*ky)/(kx**2 + ky**2)
                  else
                     ep=0.0
                  endif
                  ep = max(0.0, ep)

                  res(2,i,j) = res(2,i,j) - ep*kx*(dx*dy)
                  res(3,i,j) = res(3,i,j) - ep*ky*(dx*dy)
                  res(4,i,j) = res(4,i,j) - &
                               ep*(vex(i,j)*kx + vey(i,j)*ky)*(dx*dy)
               endif
               resid = resid + res(:,i,j)**2
            enddo
         enddo

         resid = sqrt(resid)

         do i=1,nx
            do j=1,ny
               co1(:,i,j) = ark(rks)*co0(:,i,j) + &
                            (1.0-ark(rks))*(co1(:,i,j) - lambda*res(:,i,j))
            enddo
         enddo

         call periodic(co1)

      enddo ! Rk stage loop

      if(it==1)then
         resid1 = resid
      endif

      time = time + dt
      write(*,'(I6,F10.2,4E12.4)')it,time,resid(:)/resid1(:)

      if(mod(it,itsave)==0)then
         call cons2prim(co1,rho,vex,vey,pre)
         call saveprim(time, rho, vex, vey, pre)
         call vorticity(rho, vex, vey, pre, omg)
         call savevort(time, omg)
      endif

   enddo ! time iteration loop


end subroutine solveGMD
