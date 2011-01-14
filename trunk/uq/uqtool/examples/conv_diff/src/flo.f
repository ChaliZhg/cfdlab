      program main
c ========== Program to solve u u_x = u_xx + S(x) =================
      implicit none
      include 'param.h'
      real, allocatable :: q(:),xc(:),qexact(:),qold(:),xv(:)
      real, allocatable :: res(:),dx(:)
      integer nc
      integer i, niter, iter,stage, mode
      real :: dt, maxspeed, CFL, factor, xi_1,tol_conv
      real :: frac,cost,exactcost,residue,alpha

      call read_input(mode, xi_1, CFL, niter,tol_conv)

      call bc()

c     Set up mesh
      open(10, file='grid.dat', status='old')
      read(10,*) nc
      allocate( q(nc) ,qexact(nc), qold(nc))
      allocate( res(nc))
      allocate( xv(nc+1),xc(nc),dx(nc))
      do i=1,nc+1
         read(10,*) xv(i)
      enddo
      close(10)

      do i=1,nc
         xc(i) = 0.5*(xv(i) + xv(i+1))
         dx(i) = xv(i+1) - xv(i)
      enddo

c     Read primal solution, evaluate residual, save to file
      if(mode.eq.2)then
         print*,'Reading primal.dat ...'
         open(10, file='primal.dat', status='old')
         read(10,*) (q(i),i=1,nc)
         close(10)
         res = 0.0
         call residu(nc, q, res, xc, xv, dx)
         call source(nc, q, res, xc, dx, xi_1)
         print*,'Saving primal residual ...'
         open(10, file='p_residual.dat')
         write(10,'(e24.14)') (res(i), i=1,nc)
         close(10)
         call costfun(nc, q, cost,dx)
         open(10,file='obj.dat')
         write(10,'(e24.14)') cost
         close(10)
         stop
      endif

c     initialize conditions
       
      do i=1,nc
         q(i) = 0.d0
      enddo

      open(10, file='init.dat')
      do i=1,nc
         write(10,*) xc(i), q(i)
      enddo
      close(10)


      do iter=1,niter
         qold(1:nc)=q(1:nc)
          dt=1.e20
          do i=1,nc
          dt       = min(min(CFL*dx(i)/1.5,CFL*dx(i)*dx(i)),dt)
          enddo
         do stage=1,4
            res = 0.0
            call residu(nc, q, res, xc, xv, dx)
            call source(nc, q, res, xc, xv, xi_1)
            residue = 0.d0
            do i=1,nc
               residue=residue+abs(res(i))
            enddo
            alpha=1./float(5-stage)
            do i=1,nc
            q(i) = qold(i) - alpha*(dt/dx(i))*res(i)
            enddo
         enddo
c        if(mod(iter,1).eq.0) print*,iter,residue
         if(abs(residue).lt.tol_conv) exit
      enddo
      print*,iter,residue

      open(11, file='flo.dat')
      open(12, file='exact.dat')
      open(13, file='primal.dat')
      do i=1,nc
         write(11,'(2e21.13)')xc(i), q(i)
         qexact(i)=10.0*xc(i)*(1.-xc(i))*sin(xi_1*xc(i));
         write(12,*) xc(i), qexact(i)
         write(13,'(e24.14)') q(i)
      enddo
      close(11)
      close(12)
      close(13)

      print*,'Computed'
      call costfun(nc, q, cost,dx)
      open(10,file='obj.dat')
      write(10,'(e24.14)') cost
      close(10)
      print*,'Exact (Numerical)'
      call costfun(nc, qexact, cost,dx)
      print*,'Exact (Analytical)'

      ExactCost= 5.0*(xi_1*(45.0+2.0*xi_1**4) + 45.0*xi_1*cos(2*xi_1) +
     1           15.0*(-3.0 + xi_1**2)*sin(2.0*xi_1))
      ExactCost = ExactCost/(6.0*xi_1**5)
 
      print*,'Cost= ',ExactCost

 1000  format(4(1X,E12.6))

      deallocate( q )
      deallocate( res )

      stop
      end
