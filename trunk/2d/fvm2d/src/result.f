      subroutine write_result(prim, mu, qx, qy, coord, elem, spts, dsb, 
     &                        ptype)
      implicit none
      include 'param.h'
      integer          ptype(npmax), elem(nvemax,ntmax), spts(nspmax)
      double precision prim(nvar,npmax), qx(nvar,npmax), qy(nvar,npmax),
     &                 coord(2,npmax), dsb(2,npmax), mu(npmax)

      integer          ifile, ii, is, it, i
      double precision ft1(3), ft2(3), ft3(3), coort(2,3), val1(100), 
     &                 val2(100), val3(100), delro1, delro2, delro3, 
     &                 deltat, xx0, yy0, xx1, yy1, q2, mach, cp, ent, 
     &                 nuref, ux, uy, vx, vy, nl, nx, ny, tx, ty, tw,
     &                 us, ys, cf
 
      pmin =  1.0d10
      pmax = -1.0d10
      mmin =  1.0d10
      mmax = -1.0d10
      rmin =  1.0d10
      rmax = -1.0d10
      umin =  1.0d10
      umax = -1.0d10
      vmin =  1.0d10
      vmax = -1.0d10
      emin =  1.0d10
      emax = -1.0d10
      nmin =  1.0d10
      nmax = -1.0d10
      do is=1,np
          rmin = dmin1(rmin, prim(1,is)) 
          rmax = dmax1(rmax, prim(1,is)) 
          umin = dmin1(umin, prim(2,is)) 
          umax = dmax1(umax, prim(2,is)) 
          vmin = dmin1(vmin, prim(3,is)) 
          vmax = dmax1(vmax, prim(3,is)) 
          pmin = dmin1(pmin, prim(4,is)) 
          pmax = dmax1(pmax, prim(4,is)) 
          q2   = prim(2,is)**2 + prim(3,is)**2
          mach = dsqrt(q2*prim(1,is)/(GAMMA*prim(4,is)))
          mmin = dmin1(mmin, mach)
          mmax = dmax1(mmax, mach)
          ent  = dlog10(prim(4,is)/prim(1,is)**GAMMA/ent_inf)
          emin = dmin1(emin, ent)
          emax = dmax1(emax, ent)
          nmin = dmin1(nmin, prim(5,is))
          nmax = dmax1(nmax, prim(5,is))
      enddo

      write(*,9)('-',i=1,70)
9     format(70a)
      write(*,10)mach_inf,aoa_deg,Rey,CFL
10    format(' Mach =',f6.3,', AOA =',f6.2, ', Rey = ',e10.4,
     &       ', CFL =',f5.2)
      write(*,11)iflux,ilimit,gridfile
11    format(' Flux =',i2, ',     Lim = ',i2,',    Grid= ',a30)
      write(*,9)('-',i=1,70)
      write(*,'(" Iterations        =",i12)')iter
      write(*,'(" Cl, Cd            =",2f12.6)')cl,cd
      write(*,'(" L2 residue        =",f12.6)')dlog10(res)
      write(*,'(" Linf residue      =",e16.6)')resi
      write(*,'(" Linf point        =",f12.6,f12.6,i8,i4)')
     &      coord(1,iresi), coord(2,iresi), iresi, ptype(iresi)
      write(*,'(" Minimum density   =",f12.6)')rmin
      write(*,'(" Maximum density   =",f12.6)')rmax
      write(*,'(" Minimum pressure  =",f12.6)')pmin
      write(*,'(" Maximum pressure  =",f12.6)')pmax
      write(*,'(" Minimum mach no   =",f12.6)')mmin
      write(*,'(" Maximum mach no   =",f12.6)')mmax
      write(*,'(" Minimum x vel     =",f12.6)')umin
      write(*,'(" Maximum x vel     =",f12.6)')umax
      write(*,'(" Minimum y vel     =",f12.6)')vmin
      write(*,'(" Maximum y vel     =",f12.6)')vmax
      write(*,'(" Minimum entropy   =",f12.6)')emin
      write(*,'(" Maximum entropy   =",f12.6)')emax
      write(*,'(" Minimum viscosity =",e16.6)')nmin
      write(*,'(" Maximum viscosity =",e16.6)')nmax
      call flush(6)

C Write pressure coefficient
      open(unit=10, file='WALL.DAT')
      do i=1,nsp
            is = spts(i)      
            cp = -(prim(4,is) - p_inf)/(0.5d0*r_inf*q_inf**2)

            if(iflow .ne. inviscid)then
                  ux = qx(2,is)
                  uy = qy(2,is)
                  vx = qx(3,is)
                  vy = qy(3,is)

                  nl = dsqrt(dsb(1,is)**2 + dsb(2,is)**2)
                  nx = dsb(1,is)/nl
                  ny = dsb(2,is)/nl
                  tx = ny
                  ty =-nx

                  tw = mu(is)*((ux*tx + vx*ty)*nx + (uy*tx + vy*ty)*ny)
                  cf = tw/(0.5d0*r_inf*q_inf**2)
                  us = sqrt( dabs(tw)/prim(1,is) )
                  ys = prim(1,is)*us*wd1(i)/mu(is)
                  write(10,*) coord(1,is), cp, cf, wd1(i), ys
            else
                  write(10,*) coord(1,is), cp
            endif
      enddo
      close(10)

c iso-pressure/mach/visc for gnuplot
      open(22,file='GNU.PRES',status='unknown')
      open(23,file='GNU.MACH',status='unknown')
      open(24,file='GNU.VISC',status='unknown')
      rewind(22)
      rewind(23)
      rewind(24)
         
c Normalizing value for turbulent viscosity
      if(nmax .ne. 0.0d0)then
            nuref = nmax
      else
            nuref = 1.0d0
      endif

      delro1 = (pmax-pmin)/niso
      delro2 = (mmax-mmin)/niso
      delro3 = (nmax-nmin)/nuref/niso
      do ii=1,niso+1
            val1(ii) = pmin       + (ii-1)*delro1
            val2(ii) = mmin       + (ii-1)*delro2
            val3(ii) = nmin/nuref + (ii-1)*delro3
      enddo

      do it=1,nt
            do i=1,3
                  is         = elem(i,it)
                  coort(1,i) = coord(1,is)
                  coort(2,i) = coord(2,is)
                  ft1(i)     = prim(4,is)
                  q2         = prim(2,is)**2 + prim(3,is)**2
                  ft2(i)     = dsqrt(q2*prim(1,is)/(GAMMA*prim(4,is)))
                  ft3(i)     = prim(5,is)/nuref
            enddo
            call isocont(22,ft1,coort,niso,val1)
            call isocont(23,ft2,coort,niso,val2)
            call isocont(24,ft3,coort,niso,val3)
      enddo
      close(22)
      close(23)
      close(24)

cc velocity vector for gnuplot
      ifile  = 26
      open(ifile,file='GNU.VECT',status='unknown')
      rewind(ifile)
      deltat = 1.0d-2
      do is=1,np
            xx0 = coord(1,is)
            yy0 = coord(2,is)
            xx1 = prim(2,is)*deltat
            yy1 = prim(3,is)*deltat
            write(ifile,*) xx0, yy0, xx1, yy1
      enddo
      close(ifile)
      
c Run gnuplot and Start xv if not already started
      if(xvstatus .eq. no .and. iterlast .eq. 0)then
            call system('gnuplot res.gnu')
            call system('xv -poll res.png &')
            xvstatus = yes
      else
            call system('gnuplot res.gnu &')
      endif

      return
      end
