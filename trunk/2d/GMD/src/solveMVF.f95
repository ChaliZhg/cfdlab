subroutine solveMVF(rho, vex, vey, pre, omg, co0, co1, res)

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
   real :: fd(4,nx+1,ny+1)
   real :: gd(4,nx+1,ny+1)

   real :: p(nx+1, ny+1)

   real :: vcx(-1:nx+2, -1:ny+2)
   real :: vcy(-1:nx+2, -1:ny+2)

   integer :: it, i, j, rks
   real    :: lambda
   real    :: xflux(4), yflux(4)
   real    :: dxflux(4), dyflux(4)
   real    :: time
   real    :: resid(4), resid1(4)
   real    :: etax, etay, neta, kx, ky, omg0, ep, cv, maxcv
   real    :: ane, ase, anw, asw
   real    :: bne, bse, bnw, bsw
   real    :: S2, S3, k2, p_avg
   logical :: tostop


   ! Flux function is only one choice
   fluxtype = iadv

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

      do rks=1,nrk

         maxcv = 0.0
         res = 0.0

         call cons2prim(co1,rho,vex,vey,pre)

         ! Pressure at vertices
         do i=1,nx+1
            do j=1,ny+1
               p(i,j) = (pre(i-1,j-1) + pre(i,j-1) + pre(i,j) + pre(i-1,j))/4.0
            enddo
         enddo

         ! x fluxes
         do i=0,nx
            do j=1,ny
               p_avg    = 0.5*( p(i,j) + p(i,j+1) )
               call numflux_x(co1(:,i-1,j), co1(:,i,j), co1(:,i+1,j), &
                              co1(:,i+2,j), xflux, dxflux)
               xflux(2)     = xflux(2) + p_avg
               res(:,i,j)   = res(:,i,j)   + dy*xflux(:)
               res(:,i+1,j) = res(:,i+1,j) - dy*xflux(:)
               fd(:,i+1,j)  = dxflux(:)
            enddo
         enddo

         ! y fluxes
         do j=0,ny
            do i=1,nx
               p_avg    = 0.5*( p(i,j) + p(i+1,j) )
               call numflux_y(co1(:,i,j-1), co1(:,i,j), co1(:,i,j+1), &
                              co1(:,i,j+2), yflux, dyflux)
               yflux(3)     = yflux(3) + p_avg
               res(:,i,j)   = res(:,i,j)   + dx*yflux(:)
               res(:,i,j+1) = res(:,i,j+1) - dx*yflux(:)
               gd(:,i,j+1)  = dyflux(:)
            enddo
         enddo

         ! add vorticity confinement terms
         if(vconf==yes)then
            call vorticity(rho, vex, vey, pre, omg)

            do i=1,nx
               do j=1,ny
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

                  !cv = 1.0
                  !kx = -cv*rho(i,j)*etay*omg0
                  !ky = +cv*rho(i,j)*etax*omg0

                  kx = -rho(i,j)*etay*omg0
                  ky = +rho(i,j)*etax*omg0
                  k2 = kx**2 + ky**2
                  ! Scheme I
                  !S2 = (fd(2,i+1,j) - fd(2,i,j))/dx + (gd(2,i,j+1) - gd(2,i,j))/dy
                  !S3 = (fd(3,i+1,j) - fd(3,i,j))/dx + (gd(3,i,j+1) - gd(3,i,j))/dy
                  ! Scheme II
                  !S2 = -(fd(3,i+1,j) - fd(3,i,j))/dx - (gd(3,i,j+1) - gd(3,i,j))/dy
                  !S3 = -(fd(2,i+1,j) - fd(2,i,j))/dx - (gd(2,i,j+1) - gd(2,i,j))/dy
                  ! Scheme III
                  !S2 =  (fd(2,i+1,j) - fd(2,i,j))/dx - (fd(3,i+1,j) - fd(3,i,j))/dx
                  !S3 = -(gd(2,i,j+1) - gd(2,i,j))/dy + (gd(3,i,j+1) - gd(3,i,j))/dy
                  if(k2 > 0.0)then
                     cv = (S2*kx + S3*ky)/k2
                     cv = max(0.0, cv)
                  !   !cv = min(0.5, cv)
                  else
                     cv = 0.0
                  endif
                  !kx = cv*kx
                  !ky = cv*ky
                  !maxcv = max(maxcv, cv)

                  ! For debugging
                  rho(i,j)=omg0
                  vcx(i,j)=kx
                  vcy(i,j)=ky
                  pre(i,j)=cv
                  kx=0
                  ky=0

                  res(2,i,j) = res(2,i,j) - kx*(dx*dy)
                  res(3,i,j) = res(3,i,j) - ky*(dx*dy)
                  res(4,i,j) = res(4,i,j) - &
                               (vex(i,j)*kx + vey(i,j)*ky)*(dx*dy)
               enddo
            enddo
         endif

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

      it = it + 1
      if(it==1)then
         resid1 = resid
      endif

      time = time + dt
      write(*,'(I6,F10.2,5E12.4)')it,time,resid(:)/resid1(:),maxcv

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call cons2prim(co1,rho,vex,vey,pre)
         call saveprim(time, rho, vex, vey, pre)
         !call saveprim(time, rho, vcx, vcy, pre)
         call vorticity(rho, vex, vey, pre, omg)
         call savevort(time, omg)
      endif

   enddo ! time iteration loop


end subroutine solveMVF
