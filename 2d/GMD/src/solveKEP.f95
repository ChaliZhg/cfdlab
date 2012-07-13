!------------------------------------------------------------------------------
! Kinetic energy preserving scheme for NS
!------------------------------------------------------------------------------
subroutine solveKEP(rho, vex, vey, pre, omg, co0, co1, phi, psi, &
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

   real :: ent(-1:nx+2, -1:ny+2)
   real :: tmp(-1:nx+2, -1:ny+2)

   ! Pressure at vertices
   real :: p(nx+1,ny+1)
   real :: Tx(nx+1,ny+1)
   real :: Ty(nx+1,ny+1)
   ! Shear stress at vertices
   real :: sxx(nx+1,ny+1)
   real :: sxy(nx+1,ny+1)
   real :: syy(nx+1,ny+1)

   integer :: it, i, j, rks
   real    :: xflux(4), yflux(4)
   real    :: lambda, time
   real    :: D, ux, uy, vx, vy, q
   real    :: p_avg, sxx_avg, sxy_avg, syy_avg, u_avg, v_avg, mx_avg, my_avg, H_avg
   real    :: Tx_avg, Ty_avg
   real    :: resid(4), resid1(4)
   real    :: ke, entropy
   logical :: tostop

   ! set initial condition
   call init_cond(rho, vex, vey, pre)
   call prim2cons(rho, vex, vey, pre, co1)
   call periodic(co1)
   call cons2prim(co1, rho, vex, vey, pre)
   call saveprim(0.0, rho, vex, vey, pre)
   call vorticity(rho, vex, vey, pre, omg)
   call savevort(0.0, omg)
   call timestep(rho, vex, vey, pre)

   time   = 0.0
   it     = 0

   do while(time < final_time .and. it < itmax)

      ! Exactly match final time
      tostop = .false.
      if(time + dt > final_time)then
         dt = final_time - time
         tostop = .true.
      endif

      lambda = dt/dx/dy

      co0(:,:,:) = co1(:,:,:)

      ! RK loop
      do rks=1,nrk

         res = 0

         call cons2prim(co1,rho,vex,vey,pre)
         ! Enthalpy
         ent = gamma * pre / ((gamma-1.0) * rho) + 0.5 * (vex**2 + vey**2)
         tmp = pre / (gas_const * rho)

         ! Pressure and shear stress at vertices
         do i=1,nx+1
            do j=1,ny+1
               p(i,j) = (pre(i-1,j-1) + pre(i,j-1) + pre(i,j) + pre(i-1,j))/4.0
               Tx(i,j) = 0.5*( tmp(i,j)   - tmp(i-1,j)   + tmp(i,j-1) - tmp(i-1,j-1) ) / dx
               Ty(i,j) = 0.5*( tmp(i-1,j) - tmp(i-1,j-1) + tmp(i,j)   - tmp(i,j-1)   ) / dy
               ux = 0.5*( vex(i,j)   - vex(i-1,j)   + vex(i,j-1) - vex(i-1,j-1) ) / dx
               vx = 0.5*( vey(i,j)   - vey(i-1,j)   + vey(i,j-1) - vey(i-1,j-1) ) / dx
               uy = 0.5*( vex(i-1,j) - vex(i-1,j-1) + vex(i,j)   - vex(i,j-1)   ) / dy
               vy = 0.5*( vey(i-1,j) - vey(i-1,j-1) + vey(i,j)   - vey(i,j-1)   ) / dy
               D  = ux + vy
               sxx(i,j) = 2.0 * mu * ux - (2.0/3.0) * mu * D
               syy(i,j) = 2.0 * mu * vy - (2.0/3.0) * mu * D
               sxy(i,j) = mu * (uy + vx)
            enddo
         enddo

         ! x flux
         do i=1,nx+1
            do j=1,ny
               p_avg    = 0.5*( p(i,j) + p(i,j+1) )
               Tx_avg   = 0.5*( Tx(i,j)  + Tx(i,j+1) )
               sxx_avg  = 0.5*( sxx(i,j) + sxx(i,j+1) )
               sxy_avg  = 0.5*( sxy(i,j) + sxy(i,j+1) )
               syy_avg  = 0.5*( syy(i,j) + syy(i,j+1) )
               u_avg    = 0.5*( vex(i-1,j) + vex(i,j) )
               v_avg    = 0.5*( vey(i-1,j) + vey(i,j) )
               mx_avg   = 0.5*( rho(i-1,j)*vex(i-1,j) + rho(i,j)*vex(i,j) )
               H_avg    = 0.5*( ent(i-1,j) + ent(i,j) )
               q        = -kth * Tx_avg
               xflux(1) =                 mx_avg
               xflux(2) = p_avg + u_avg * mx_avg - sxx_avg
               xflux(3) =         v_avg * mx_avg - sxy_avg
               xflux(4) =         H_avg * mx_avg - (sxx_avg * u_avg + sxy_avg * v_avg) + q
               res(:,i-1,j) = res(:,i-1,j) + xflux(:) * dy
               res(:,i  ,j) = res(:,i  ,j) - xflux(:) * dy
            enddo
         enddo

         ! y flux
         do j=1,ny+1
            do i=1,nx
               p_avg    = 0.5*( p(i,j) + p(i+1,j) )
               Ty_avg   = 0.5*( Ty(i,j)  + Ty(i+1,j) )
               sxx_avg  = 0.5*( sxx(i,j) + sxx(i+1,j) )
               sxy_avg  = 0.5*( sxy(i,j) + sxy(i+1,j) )
               syy_avg  = 0.5*( syy(i,j) + syy(i+1,j) )
               u_avg    = 0.5*( vex(i,j-1) + vex(i,j) )
               v_avg    = 0.5*( vey(i,j-1) + vey(i,j) )
               my_avg   = 0.5*( rho(i,j-1)*vey(i,j-1) + rho(i,j)*vey(i,j) )
               H_avg    = 0.5*( ent(i,j-1) + ent(i,j) )
               q        = -kth * Ty_avg
               yflux(1) =                 my_avg
               yflux(2) =         u_avg * my_avg - sxy_avg
               yflux(3) = p_avg + v_avg * my_avg - syy_avg
               yflux(4) =         H_avg * my_avg - (sxy_avg * u_avg + syy_avg * v_avg) + q
               res(:,i,j-1) = res(:,i,j-1) + yflux(:) * dx
               res(:,i,  j) = res(:,i  ,j) - yflux(:) * dx
            enddo
         enddo

         resid = 0

         ! Update solution
         do i=1,nx
            do j=1,ny
               co1(:,i,j) = ark(rks)*co0(:,i,j) + &
                            (1.0-ark(rks))*(co1(:,i,j) - lambda*res(:,i,j))
               resid = resid + res(:,i,j)**2
            enddo
         enddo

         resid = sqrt(resid)

         call periodic(co1)

      enddo ! Rk stage loop

      it = it + 1

      if(it==1)then
         resid1 = resid
      endif

      call cons2prim(co1,rho,vex,vey,pre)
      call global_quantities(rho,vex,vey,pre,ke,entropy)

      time = time + dt
      write(*,'(I6,F10.4,6E12.4)')it,time,resid(:)/resid1(:),ke,entropy

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call saveprim(time, rho, vex, vey, pre)
         call vorticity(rho, vex, vey, pre, omg)
         call savevort(time, omg)
      endif

   enddo ! time iteration loop


end subroutine solveKEP
