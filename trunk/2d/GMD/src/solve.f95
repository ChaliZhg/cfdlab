subroutine solve(rho, vex, vey, pre, co0, co1, phi, psi, res)

   use comvar

   implicit none

   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: pre(-1:nx+2, -1:ny+2)
   real :: co0(4, -1:nx+2, -1:ny+2)
   real :: co1(4, -1:nx+2, -1:ny+2)
   real :: phi(4,nx+1,ny+1)
   real :: psi(4,nx+1,ny+1)
   real :: res(4,nx,ny)

   integer :: it, i, j, rks
   real    :: lambda
   real    :: xflux(4), yflux(4)
   real    :: time
   real    :: resid(4), resid1(4)


   ! set initial condition
   call init_cond(rho, vex, vey, pre)
   call prim2cons(rho, vex, vey, pre, co1)
   call periodic(co1)
   !call saveprim(rho, vex, vey, pre)
   call savevort(rho, vex, vey, pre)
   call timestep(rho, vex, vey, pre)

   lambda = dt/dx/dy
   time   = 0.0

   do it=1,itmax

      co0(:,:,:) = co1(:,:,:)

      do rks=1,3

         ! symmetric potential phi
         phi = 0.0

         do i=0,nx
            j = 0
            call numflux_x(co1(:,i-1,j), co1(:,i,j), &
                           co1(:,i+1,j), co1(:,i+2,j), xflux)
            phi(:,i+1,j+1) = xflux(:)
            j = ny + 1
            call numflux_x(co1(:,i-1,j), co1(:,i,j), &
                           co1(:,i+1,j), co1(:,i+2,j), xflux)
            phi(:,i+1,j) = xflux(:)
         enddo

         do i=0,nx
            do j=1,ny
               call numflux_x(co1(:,i-1,j), co1(:,i,j), &
                              co1(:,i+1,j), co1(:,i+2,j), xflux)
               phi(:,i+1,j)   = phi(:,i+1,j)   + xflux(:)
               phi(:,i+1,j+1) = phi(:,i+1,j+1) + xflux(:)
            enddo
         enddo

         ! symmetric potential psi
         psi = 0.0

         do j=0,ny
            i = 0
            call numflux_y(co1(:,i,j-1), co1(:,i,j), &
                           co1(:,i,j+1), co1(:,i,j+2), yflux)
            psi(:,i+1,j+1) = yflux(:)

            i = nx + 1
            call numflux_y(co1(:,i,j-1), co1(:,i,j), &
                           co1(:,i,j+1), co1(:,i,j+2), yflux)
            psi(:,i,j+1) = yflux(:)
         enddo

         do i=1,nx
            do j=0,ny
               call numflux_y(co1(:,i,j-1), co1(:,i,j), &
                              co1(:,i,j+1), co1(:,i,j+2), yflux)
               psi(:,i,  j+1) = psi(:,i,  j+1) + yflux(:)
               psi(:,i+1,j+1) = psi(:,i+1,j+1) + yflux(:)
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
               resid = resid + res(:,i,j)**2
            enddo
         enddo
         res = 0.5 * 0.5 * res
         resid = sqrt(resid)

         do i=1,nx
            do j=1,ny
               co1(:,i,j) = ark(rks)*co0(:,i,j) + &
                            (1.0-ark(rks))*(co1(:,i,j) - lambda*res(:,i,j))
            enddo
         enddo

         call periodic(co1)

      enddo ! Rk stage loop

      if(mod(it,itsave)==0)then
         call cons2prim(co1,rho,vex,vey,pre)
         !call saveprim(rho, vex, vey, pre)
         call savevort(rho, vex, vey, pre)
      endif

      if(it==1)then
         resid1 = resid
      endif

      time = time + dt
      write(*,'(I6,F10.2,4E12.4)')it,time,resid(:)/resid1(:)

   enddo ! time iteration loop


end subroutine solve
