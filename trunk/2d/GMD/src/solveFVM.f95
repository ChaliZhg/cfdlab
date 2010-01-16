subroutine solveFVM(rho, vex, vey, pre, omg, co0, co1, res)

   use comvar

   implicit none

   real :: rho(-1:nx+2, -1:ny+2)
   real :: vex(-1:nx+2, -1:ny+2)
   real :: vey(-1:nx+2, -1:ny+2)
   real :: pre(-1:nx+2, -1:ny+2)
   real :: omg( 1:nx+1,  1:ny+1)
   real :: co0(4, -1:nx+2, -1:ny+2)
   real :: co1(4, -1:nx+2, -1:ny+2)
   real :: res(4,0:nx+1,0:ny+1)

   integer :: it, i, j, rks
   real    :: lambda
   real    :: xflux(4), yflux(4)
   real    :: dxflux(4), dyflux(4)
   real    :: time
   real    :: resid(4), resid1(4)


   ! set initial condition
   call init_cond(rho, vex, vey, pre)
   call prim2cons(rho, vex, vey, pre, co1)
   call periodic(co1)
   call saveprim(rho, vex, vey, pre)
   call vorticity(rho, vex, vey, pre, omg)
   call savevort(omg)
   call timestep(rho, vex, vey, pre)

   lambda = dt/dx/dy
   time   = 0.0

   do it=1,itmax

      co0(:,:,:) = co1(:,:,:)

      do rks=1,3

         res = 0.0

         ! x fluxes
         do i=0,nx
            do j=1,ny
               call numflux_x(co1(:,i-1,j), co1(:,i,j), co1(:,i+1,j), &
                              co1(:,i+2,j), xflux, dxflux)
               res(:,i,j)   = res(:,i,j)   + dy*xflux(:)
               res(:,i+1,j) = res(:,i+1,j) - dy*xflux(:)
            enddo
         enddo

         ! y fluxes
         do j=0,ny
            do i=1,nx
               call numflux_y(co1(:,i,j-1), co1(:,i,j), co1(:,i,j+1), &
                              co1(:,i,j+2), yflux, dyflux)
               res(:,i,j)   = res(:,i,j)   + dx*yflux(:)
               res(:,i,j+1) = res(:,i,j+1) - dx*yflux(:)
            enddo
         enddo

         ! update conserved variables
         resid = 0.0
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

      if(mod(it,itsave)==0)then
         call cons2prim(co1,rho,vex,vey,pre)
         call saveprim(rho, vex, vey, pre)
         call vorticity(rho, vex, vey, pre, omg)
         call savevort(omg)
      endif

      if(it==1)then
         resid1 = resid
      endif

      time = time + dt
      write(*,'(I6,F10.2,4E12.4)')it,time,resid(:)/resid1(:)

   enddo ! time iteration loop


end subroutine solveFVM
