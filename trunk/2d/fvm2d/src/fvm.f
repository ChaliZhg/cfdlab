C Cell-vertex finite volume code on triangles
C Performs inviscid, laminar and turbulent computations
C Spalart Allmaras turbulence model
C This is the main program for Finite Volume Solver
      program fvm
      implicit none
      include 'param.h'

      integer          ptype(npmax), elem(nvemax,ntmax), edge(2,nemax),
     &                 spts(nspmax), fpts(nfpmax), opts(nopmax),
     &                 bpts(nbpmax)
      double precision coord(2, npmax), elarea(ntmax), 
     &                 cvarea(npmax), ds(2,nemax), dsb(2,npmax), 
     &                 dt(npmax), divf(nvar,npmax), prim(nvar,npmax), 
     &                 prim_old(nvar,npmax), qx(nvar,npmax), 
     &                 qy(nvar,npmax), drmin(npmax), mu(npmax), 
     &                 mul(npmax), cvareamc(npmax), wd(npmax)

      integer          fres, irk

      call banner
      call math
      call read_input
      call read_grid(coord, elem, ptype, spts, fpts, opts, bpts)
      call prep_gnuplot
      call geometric(ptype, elem, edge, coord, elarea, cvarea, ds, dsb, 
     &               cvareamc, wd, drmin, spts, fpts, opts, bpts)
      call initialize(prim, mul, mu)

      iter = iterlast
      res  = MINRESIDUE + 1.0d0
      res1 = 0.0d0

      print*,'Beginning of iterations ...'

      fres = 11
      call system("rm -f RES.DAT")
      open(unit=fres, file='RES.DAT')
      rewind(fres)
      do while(res .gt. MINRESIDUE .and. iter .lt. iterlast+MAXITER)

         call save_old(prim, prim_old)
         call time_step(prim, mu, drmin, dt)

         do irk=1,NIRK
            call fvm_residual(ptype, elem, edge, spts, fpts, opts, bpts,
     &                        coord, elarea, cvarea, ds, dsb, dt, divf,
     &                        prim, prim_old, qx, qy, drmin, mu, mul, 
     &                        cvareamc, wd)

            call update(irk, divf, prim_old, prim, cvarea, spts, dt)
         enddo

         if(iflow .eq. turbulent)then
            call sa_model(elem, edge, spts, fpts, opts, coord, elarea,
     &                    cvarea, prim, prim_old, wd, mul, ds,
     &                    dsb, dt, qx, qy)
         endif

         call check_positivity(prim, coord, ptype)

         call clcd(spts, dsb, prim, mu, qx, qy)

         call residue(prim, prim_old, dt)
         write(fres,'(I8,E14.6,I8,3E14.6)')iter,res,iresi,resi,cl,cd

         iter = iter + 1

         if(mod(iter,saveinterval) .eq. 0)then
            call flush(fres)
            call write_result(prim, mu, qx, qy, coord, elem, spts,
     &                        dsb, ptype)
            if(iter .ge. MAXITER)then
               print*,'*** Stopping execution ***'
               print*,'\tMaximum number of iterations reached'
               goto 100
            endif
         endif
      enddo
      close(fres)

      call write_result(prim, mu, qx, qy, coord, elem, spts, dsb, ptype)

      print*,'*** Stopping execution ***'
      print*,'\tResidue has been reduced below MINRES =',MINRESIDUE

100   call finalize(prim, coord, elem, ptype)

      stop
      end
