subroutine solveFVM(pri, co0, co1, res, divB)

   use comvar
   use omp_lib

   implicit none

   real :: pri(nvar, -1:nx+2, -1:ny+2)
   real :: co0(nvar, -1:nx+2, -1:ny+2)
   real :: co1(nvar, -1:nx+2, -1:ny+2)
   real :: res(nvar,  0:nx+1,  0:ny+1)
   real :: divB(1:nx+1, 1:ny+1)

   integer :: it, i, j, rks
   real    :: lambda
   real    :: xflux(nvar), yflux(nvar)
   real    :: time
   real    :: resid(nvar), resid1(nvar)
   real    :: ke, entropy
   real    :: div, source(nvar), maxdivB, divf_l(nvar), divf_r(nvar)
   real    :: minmod, minmod2, Bx, By
   logical :: tostop

   call omp_set_num_threads(4)

   ! set initial condition
   call init_cond(pri, co1)
   call periodic(co1)
   call cons2prim(co1, pri)
   call saveprim(0.0, pri)

   time   = 0.0
   it     = 0

   do while(time < final_time .and. it < itmax)

      call timestep(pri)

      ! Exactly match final time
      tostop = .false.
      if(time + dt > final_time)then
         dt = final_time - time
         tostop = .true.
      endif

      lambda = dt/dx/dy

      co0(:,:,:) = co1(:,:,:)

      do rks=1,nrk

         res = 0.0

         ! x fluxes
         !$omp parallel do private(xflux) shared(res)
         do i=0,nx
            do j=1,ny
               call numflux(1.0, 0.0,     &
                            pri(:,i-1,j), &
                            pri(:,i,j),   &
                            pri(:,i+1,j), &
                            pri(:,i+2,j), &
                            xflux)
               !call divflux(1.0, 0.0,     &
               !             pri(:,i-1,j), &
               !             pri(:,i,j),   &
               !             pri(:,i+1,j), &
               !             pri(:,i+2,j), &
               !             divf_l, divf_r)
               !res(:,i,j)   = res(:,i,j)   + dy*(xflux(:) + divf_l(:))
               !res(:,i+1,j) = res(:,i+1,j) - dy*(xflux(:) + divf_r(:))
               res(:,i,j)   = res(:,i,j)   + dy*xflux(:)
               res(:,i+1,j) = res(:,i+1,j) - dy*xflux(:)
            enddo
         enddo
         !$omp end parallel do

         ! y fluxes
         !$omp parallel do private(yflux) shared(res)
         do j=0,ny
            do i=1,nx
               call numflux(0.0, 1.0,     &
                            pri(:,i,j-1), &
                            pri(:,i,j),   &
                            pri(:,i,j+1), &
                            pri(:,i,j+2), &
                            yflux)
               !call divflux(0.0, 1.0,     &
               !             pri(:,i,j-1), &
               !             pri(:,i,j),   &
               !             pri(:,i,j+1), &
               !             pri(:,i,j+2), &
               !             divf_l, divf_r)
               !res(:,i,j)   = res(:,i,j)   + dx*(yflux(:) + divf_l(:))
               !res(:,i,j+1) = res(:,i,j+1) - dx*(yflux(:) + divf_r(:))
               res(:,i,j)   = res(:,i,j)   + dx*yflux(:)
               res(:,i,j+1) = res(:,i,j+1) - dx*yflux(:)
            enddo
         enddo
         !$omp end parallel do

         ! update conserved variables
         resid = 0.0
         !$omp parallel do private(Bx,By,div,source) shared(resid)
         do i=1,nx
            do j=1,ny
               !source = 0.0
               ! Limited divergence
               !Bx = minmod2(2.0*(pri(6,i,j)-pri(6,i-1,j)), &
               !            0.5*(pri(6,i+1,j)-pri(6,i-1,j)), &
               !            2.0*(pri(6,i+1,j)-pri(6,i,j)))
               !By = minmod2(2.0*(pri(7,i,j)-pri(7,i,j-1)), &
               !            0.5*(pri(7,i,j+1)-pri(7,i,j-1)), &
               !            2.0*(pri(7,i,j+1)-pri(7,i,j)))
               !div = Bx/dx + By/dy
               ! Central divergence
               div = 0.5*( pri(6,i+1,j) - pri(6,i-1,j) ) / dx + &
                     0.5*( pri(7,i,j+1) - pri(7,i,j-1) ) / dy
               source(1) = 0.0
               source(2:4) = pri(6:8,i,j)
               source(5) = pri(2,i,j)*pri(6,i,j) + &
                           pri(3,i,j)*pri(7,i,j) + &
                           pri(4,i,j)*pri(8,i,j)
               source(6:8) = pri(2:4,i,j)
               source = source * div

               co1(:,i,j) = ark(rks)*co0(:,i,j) + &
                            brk(rks)*(co1(:,i,j) - lambda*res(:,i,j) - dt*source)
               resid = resid + res(:,i,j)**2
            enddo
         enddo
         !$omp end parallel do
         resid = sqrt(resid)

         call periodic(co1)
         call cons2prim(co1,pri)

      enddo ! Rk stage loop

      it = it + 1
      if(it==1)then
         resid1 = 1.0
         do i=1,nvar
            if(resid(i) > 0.0) resid1(i) = resid(i)
         enddo
      endif

      !call global_quantities(rho,vex,vey,pre,ke,entropy)
      call compute_divB(pri,divB,maxdivB)

      time = time + dt
      write(*,'(I6,E12.3,8E12.4)')it,time,resid(:)/resid1(:)
      write(*,'("     Div B = ", E12.3)') maxdivB
      call flush()

      if(mod(it,itsave)==0 .or. it==itmax .or. tostop)then
         call saveprim(time, pri)
      endif

   enddo ! time iteration loop


end subroutine solveFVM
