subroutine solveFVM(rho, vex, vey, vez, pre, omg, co0, co1, res)

   use comvar

   implicit none

   real :: rho(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vex(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vey(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: vez(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: pre(-1:nx+2, -1:ny+2, -1:nz+2)
   real :: omg( 1:nx+1,  1:ny+1,  1:nz+1)
   real :: co0(nvar, -1:nx+2, -1:ny+2, -1:nz+2)
   real :: co1(nvar, -1:nx+2, -1:ny+2, -1:nz+2)
   real :: res(nvar,0:nx+1,0:ny+1,0:nz+1)

   integer :: it, i, j, k, rks
   real    :: lambda
   real    :: xflux(nvar), yflux(nvar), zflux(nvar)
   real    :: time
   real    :: resid(nvar), resid1(nvar)
   real    :: ke, entropy, ke0, entropy0
   logical :: tostop


   ! set initial condition
   call init_cond(rho, vex, vey, vez, pre)
   call prim2cons(rho, vex, vey, vez, pre, co1)
   call periodic(co1)
   call cons2prim(co1, rho, vex, vey, vez, pre)
   call saveprim(0.0, rho, vex, vey, vez, pre)
   call vorticity(rho, vex, vey, vez, pre, omg)
   call savevort(0.0, omg)
   call timestep(rho, vex, vey, vez, pre)

   time   = 0.0
   it     = 0

   do while(time < final_time .and. it < itmax)

      ! Exactly match final time
      tostop = .false.
      if(time + dt > final_time)then
         dt = final_time - time
         tostop = .true.
      endif

      lambda = dt/dx/dy/dz

      co0(:,:,:,:) = co1(:,:,:,:)

      do rks=1,nrk

         res = 0.0

         ! x fluxes
         do i=0,nx
            do j=1,ny
               do k=1,nz
                  call numflux_x(co1(:,i-1,j,k), co1(:,i,j,k), co1(:,i+1,j,k), &
                                 co1(:,i+2,j,k), xflux)
                  res(:,i,j,k)   = res(:,i,j,k)   + dy*dz*xflux(:)
                  res(:,i+1,j,k) = res(:,i+1,j,k) - dy*dz*xflux(:)
               enddo
            enddo
         enddo

         ! y fluxes
         do j=0,ny
            do k=1,nz
               do i=1,nx
                  call numflux_y(co1(:,i,j-1,k), co1(:,i,j,k), co1(:,i,j+1,k), &
                                 co1(:,i,j+2,k), yflux)
                  res(:,i,j,k)   = res(:,i,j,k)   + dx*dz*yflux(:)
                  res(:,i,j+1,k) = res(:,i,j+1,k) - dx*dz*yflux(:)
               enddo
            enddo
         enddo

         ! z fluxes
         do k=0,nz
            do i=1,nx
               do j=1,ny
                  call numflux_z(co1(:,i,j,k-1), co1(:,i,j,k), co1(:,i,j,k+1), &
                                 co1(:,i,j,k+2), zflux)
                  res(:,i,j,k)   = res(:,i,j,k)   + dx*dy*zflux(:)
                  res(:,i,j,k+1) = res(:,i,j,k+1) - dx*dy*zflux(:)
               enddo
            enddo
         enddo

         ! update conserved variables
         resid = 0.0
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  co1(:,i,j,k) = ark(rks)*co0(:,i,j,k) + &
                              (1.0-ark(rks))*(co1(:,i,j,k) - lambda*res(:,i,j,k))
                  resid = resid + res(:,i,j,k)**2
               enddo
            enddo
         enddo
         resid = sqrt(resid)

         call periodic(co1)

      enddo ! Rk stage loop

      call cons2prim(co1,rho,vex,vey,vez,pre)
      call global_quantities(rho,vex,vey,vez,pre,ke,entropy)

      it = it + 1
      if(it==1)then
         resid1 = resid
         ke0 = ke
         entropy0 = entropy
      endif

      time = time + dt
      !write(*,'(I6,F10.2,7E12.4)')it,time,resid(:)/resid1(:),ke,entropy
      !write(*,'(I6,F10.2,7E12.4)')it,time,resid(:),ke,entropy
      write(*,'(I6,4E18.8)')it,dt,time,ke/ke0,entropy/entropy0

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call saveprim(time, rho, vex, vey, vez, pre)
         call vorticity(rho, vex, vey, vez, pre, omg)
         call savevort(time, omg)
      endif

   enddo ! time iteration loop


end subroutine solveFVM
