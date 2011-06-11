C----------------------------------------------------------------------
C.....Definition of some constants
C----------------------------------------------------------------------
      subroutine math
      implicit none
      include 'param.h'

      GAMMA        = 1.4d0
      GAMMA1       = GAMMA-1.0d0
      GAS_CONST    = 1.0d0
      M_PI         = 4.0d0*datan2(1.0d0, 1.0d0)

      prandtl      = 0.72d0
      prandtl_turb = 0.9d0

      ALBADA11     = 2.0d0/3.0d0
      ALBADA12     = 1.0d0 - ALBADA11

      ALBADA21     = 4.0d0/3.0d0
      ALBADA22     = 1.0d0 - ALBADA21

      xvstatus     = no

c Constants for Spalart-Allmaras model
      Cb1          = 0.1355d0
      Cb2          = 0.622d0
      sigma_sa     = 2.0d0/3.0d0
      kolm         = 0.41d0
      Cw1          = Cb1/kolm**2 + (1.0d0 + Cb2)/sigma_sa
      Cw2          = 0.3d0
      Cw3          = 2.0d0
      Cv1          = 7.1d0
      Cv2          = 5.0d0

      Cv11         = Cv1**3
      Cw31         = 1.0d0 + Cw3**6
      Cw32         = Cw3**6
      kolm2        = kolm**2
      Cb2Sig1      = (1.0d0 + Cb2)/sigma_sa
      Cb2Sig2      = Cb2/sigma_sa

      return
      end

C----------------------------------------------------------------------
C.....Rotate a vector (u,v) by angle ang
C----------------------------------------------------------------------
      subroutine rotate(u, v, ang)
      implicit none
      double precision u, v, ang

      double precision ut, vt, ct, st

      ut = u
      vt = v

      ct = dcos(ang)
      st = dsin(ang)

      u  = ut*ct + vt*st
      v  =-ut*st + vt*ct

      return
      end

C----------------------------------------------------------------------
C.....Error function, from Abromovitz and Stegun
C----------------------------------------------------------------------
      double precision function ERRF(X)
      double precision X,ARG,E,VB,T,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6

      ARG = X*X
      if(ARG .lt. 20.0d0)then
            E = dexp(-ARG)
      else
            E = 0.0d0
      endif
      VB = dabs(X)
      T = 1.0d0/(1.0d0 + 0.3275911d0*VB)
      tmp1 = 1.061405429d0*T
      tmp2 = (tmp1 - 1.453152027d0)*T
      tmp3 = (tmp2 + 1.421413741d0)*T
      tmp4 = (tmp3 - 0.284496736d0)*T
      tmp5 = (tmp4 + 0.254829592d0)*T
      tmp6 = 1.0d0 - tmp5*E
      if(X .lt. 0.0d0)then
            ERRF = -tmp6
      else
            ERRF =  tmp6
      endif

      return
      end

      REAL*8 FUNCTION DRAND(iseed)
c     -----------------------------------------------------------------
c     Selection aleatoire d'un nombre entre 0 et 1 suivant une
c     valeur donnee iseed
c
c                  mod(iseed*7141 + 54773, 259200)
c         ran = -----------------------------------
c                            259200
c     -----------------------------------------------------------------
c
c     Parametres d'appel 
c
      INTEGER iseed
c
c     Variables locales 
c
      INTEGER ia, ic, im
      PARAMETER(ia = 7141, ic = 54773, im = 259200)
c
      iseed    = ABS(MOD(iseed*ia+ic, im))
c
      DRAND    = FLOAT(iseed)/FLOAT(im)
c
      RETURN
      END

C----------------------------------------------------------------------
C.....Read grid data from a file
C.....Currently supports only triangular elements
C----------------------------------------------------------------------
      subroutine read_grid(coord, elem, ptype, spts, fpts, opts, bpts)
      implicit none
      include 'param.h'
      double precision coord(2, npmax)
      integer          elem(nvemax, ntmax), ptype(npmax), spts(nspmax),
     &                 fpts(nfpmax), opts(nopmax), bpts(nbpmax)

      integer          ngrid, ip, i, j

      print*,'Reading grid from file ',gridfile

      ngrid = 10
      open(ngrid, file=gridfile, status="old")
      rewind(ngrid)
      read(ngrid,*) np, nt
      print*, 'Number of points    : ', np
      print*, 'Number of triangles : ', nt

      if(np.gt.npmax) then
            print*, 'Increase the size of npmax'
            stop
      endif

      if(nt.gt.ntmax) then
            print*, 'Increase the size of ntmax'
            stop
      endif

      do ip=1,np
            read(ngrid,*) i, coord(1,ip), coord(2,ip), ptype(ip)
      enddo

      do ip=1,nt
            read(ngrid,*) i, elem(1, ip), elem(2, ip), elem(3, ip), j
      enddo

      close(ngrid)

c     Find bounding box
      xmin = 1000000.0d0
      ymin = 1000000.0d0
      xmax =-1000000.0d0
      ymax =-1000000.0d0
      do ip=1,np
            xmin = dmin1(xmin, coord(1,ip))
            ymin = dmin1(ymin, coord(2,ip))

            xmax = dmax1(xmax, coord(1,ip))
            ymax = dmax1(ymax, coord(2,ip))
      enddo

      print*,'Bounding box:'
      print*, '\t\txmin = ', xmin
      print*, '\t\txmax = ', xmax
      print*, '\t\tymin = ', ymin
      print*, '\t\tymax = ', ymax

      nsp = 0
      nfp = 0
      nop = 0
      nbp = 0

      do ip=1,np
            if(ptype(ip) .eq. solid)then
                  nsp       = nsp + 1
                  spts(nsp) = ip
            endif
            if(ptype(ip) .eq. farfield)then
                  nfp       = nfp + 1
                  fpts(nfp) = ip
            endif
            if(ptype(ip) .eq. outflow)then
                  nop       = nop + 1
                  opts(nop) = ip
            endif
            if(ptype(ip) .ne. interior)then
                  nbp       = nbp + 1
                  bpts(nbp) = ip
            endif
      enddo

      print*,'Number of solid    points = ',nsp
      print*,'Number of farfield points = ',nfp
      print*,'Number of outflow  points = ',nop
      print*,'Number of boundary points = ',nbp

      if( nsp+nfp+nop .ne. nbp )then
            print*,'There seem to be some unrecognized point types'
            stop
      endif


      return
      end

C----------------------------------------------------------------------
C.....Variables stored are primitive - density, u, v, pressure
C.....Initialize primitive variables to free stream values
C----------------------------------------------------------------------
      subroutine initialize(prim, nut, mul, mu)
      implicit none
      include 'param.h'
      double precision prim(nvar, npmax), nut(npmax), mul(npmax), 
     +                 mu(npmax)
      
      integer          i, j
      double precision u1, u2, u3, u4, u5

      q_inf   = 1.0d0
      r_inf   = 1.0d0
      p_inf   = 1.0d0/(GAMMA*mach_inf**2)
      T_inf   = p_inf/(r_inf*GAS_CONST)
      aoa     = aoa_deg*M_PI/180.0d0
      u_inf   = q_inf*dcos(aoa)
      v_inf   = q_inf*dsin(aoa)
      ent_inf = p_inf/r_inf**GAMMA
      a_inf   = dsqrt(GAMMA*p_inf/r_inf)
      H_inf   = a_inf**2/GAMMA1 + 0.5d0*q_inf

c Required by Sutherland law
      T_infd  = 300.0d0
      SCONST  = 110.4d0*T_inf/T_infd

c Primitive variables in free-stream
      prim_inf(1) = r_inf
      prim_inf(2) = u_inf
      prim_inf(3) = v_inf
      prim_inf(4) = p_inf
      nut_inf     = 0.0d0

      if(iflow.eq.inviscid) print*,'Euler computation'
      if(iflow.eq.laminar)  print*,'Laminar Navier-Stokes computation'
      if(iflow.eq.turbulent)print*,'Turbulent Navier-Stokes computation'
      print*,'Free-stream values:'
      print*,'\t\t Mach number =',mach_inf
      print*,'\t\t AOA         =',aoa_deg
      print*,'\t\t u velocity  =',u_inf
      print*,'\t\t v velocity  =',v_inf
      print*,'\t\t Pressure    =',p_inf

      if(vortex .eq. yes)then
            print*,'Using point vortex correction for far-field points'
            print*,'\tVortex center = ',xref, yref
      endif

C Runge-Kutta time stepping
      NIRK    = 3
      airk(1) = 0.0d0
      airk(2) = 3.0d0/4.0d0
      airk(3) = 1.0d0/3.0d0
      birk(1) = 1.0d0
      birk(2) = 1.0d0/4.0d0
      birk(3) = 2.0d0/3.0d0

      if(istart .eq. scratch)then
            print*,'Initializing solution to free stream values'
            do j=1,np
                  do i=1,nvar
                     prim(i,j) = prim_inf(i)
                  enddo
            enddo

            if(iflow .eq. turbulent)then
                  call sutherland(prim, mul)
                  do i=1,np
                     nut(i) = 0.1d0*mul(i)/prim(1,i)
                  enddo
                  call viscosity(prim, nut, mul, mu)
                  nut_inf = nut(1)
            elseif(iflow .eq. laminar)then
                  call sutherland(prim, mul)
                  do i=1,np
                     nut(i) = 0.0d0
                  enddo
                  call viscosity(prim, nut, mul, mu)
                  nut_inf = 0.0d0
            else
                  do i=1,np
                        mul(i)    = 0.0d0
                        mu(i)     = 0.0d0
                        nut(i) = 0.0d0
                  enddo
                  nut_inf = 0.0d0
            endif

      else
            print*,'Initializing solution to old values from SOL'
            open(unit=20, file='SOL', status='old')
            do j=1,np
                  read(20,*) u1, u2, u3, u4, u5
                  prim(1,j) = u1
                  prim(2,j) = u2/u1
                  prim(3,j) = u3/u1
                  prim(4,j) = GAMMA1*( u4 - 0.5d0*(u2**2 + u3**2)/u1 )
                  nut(j)    = u5
            enddo
            close(20)
      endif

      return
      end

C----------------------------------------------------------------------
C.....Calculate lift and drag coefficients
C----------------------------------------------------------------------
      subroutine clcd(coord, prim, mu, qx, qy)
      implicit none
      include 'param.h'
      include 'gdata.h'
      double precision coord(2,npmax), prim(nvar,npmax), mu(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax)

      integer          i, ie, ip1, ip2
      double precision dx, dy, sx, sy, xf, yf, p, p1, p2, txx, txy, tyy,
     &                 txx1, txy1, tyy1, txx2, txy2, tyy2

      xf = 0.0d0
      yf = 0.0d0
      do i=nswe1,nswe2
         ie   = beindx(i)
         ip1  = edge(1,ie)
         ip2  = edge(2,ie)
         dx   = coord(1,ip2) - coord(1,ip1)
         dy   = coord(2,ip2) - coord(2,ip1)
         sx   = dy
         sy   =-dx

         p1   = prim(4,ip1)
         txx1 = 2.0d0/3.0d0*mu(ip1)*(2.0d0*qx(2,ip1) - qy(3,ip1))
         txy1 =             mu(ip1)*(      qy(2,ip1) + qx(3,ip1))
         tyy1 = 2.0d0/3.0d0*mu(ip1)*(2.0d0*qy(3,ip1) - qx(2,ip1))

         p2   = prim(4,ip2)
         txx2 = 2.0d0/3.0d0*mu(ip2)*(2.0d0*qx(2,ip2) - qy(3,ip2))
         txy2 =             mu(ip2)*(      qy(2,ip2) + qx(3,ip2))
         tyy2 = 2.0d0/3.0d0*mu(ip2)*(2.0d0*qy(3,ip2) - qx(2,ip2))

         p    = 0.5d0*(p1 + p2)
         txx  = 0.5d0*(txx1 + txx2)
         txy  = 0.5d0*(txy1 + txy2)
         tyy  = 0.5d0*(tyy1 + tyy2)

         xf = xf + p*sx - txx*sx - txy*sy
         yf = yf + p*sy - txy*sx - tyy*sy
      enddo

      cl =(-dsin(aoa)*xf + dcos(aoa)*yf)/(0.5d0*r_inf*q_inf**2)
      cd =( dcos(aoa)*xf + dsin(aoa)*yf)/(0.5d0*r_inf*q_inf**2)

      return
      end

C----------------------------------------------------------------------
C.....Compute local time step using cfl condition
C----------------------------------------------------------------------
      subroutine time_step(prim, mu, drmin, dt)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), mu(npmax), drmin(npmax), 
     &                 dt(npmax)

      integer          i
      double precision q, a, dtv

      if(iflow .ne. inviscid)then
c        Time-step for viscous flow
         do i=1,np
            q        = dsqrt(prim(2,i)**2 + prim(3,i)**2)
            a        = dsqrt(GAMMA*prim(4,i)/prim(1,i))
            dtv      = 2.0d0*GAMMA*mu(i)/(prim(1,i)*prandtl)
            dt(i)    = CFL*drmin(i)**2/(drmin(i)*(q + a) + dtv)
         enddo
      else
c        Time-step for inviscid flow
         do i=1,np
            q        = dsqrt(prim(2,i)**2 + prim(3,i)**2)
            a        = dsqrt(GAMMA*prim(4,i)/prim(1,i))
            dt(i)    = CFL*drmin(i)/(q + a)
         enddo
      endif

c     Find global time-step
      dtglobal = 1.0d10
      do i=1,np
         dtglobal = dmin1(dtglobal, dt(i))
      enddo

      return
      end

C----------------------------------------------------------------------
C.....Read parameters from an input file and set freestream values
C----------------------------------------------------------------------
      subroutine read_input
      implicit none
      include 'param.h'
      integer inp, iargc, n, inpstatus
      character sdummy*32

      n = iargc()
      if(n .eq. 0)then
            print*,'You must specify an input file.'
            stop
      endif

      call getarg(1,inpfile)

      inp = 11
      open(unit=inp, file=inpfile, status='old')
      print*,'Reading parameters from ',inpfile
      read(inp,*)sdummy, istart
      read(inp,*)sdummy, iflow
      read(inp,*)sdummy, mach_inf
      read(inp,*)sdummy, aoa_deg
      read(inp,*)sdummy, Rey
      read(inp,*)sdummy, cfl
      read(inp,*)sdummy, iterlast
      read(inp,*)sdummy, maxiter
      read(inp,*)sdummy, minresidue
      read(inp,*)sdummy, saveinterval
      read(inp,*)sdummy, niso
      read(inp,*)sdummy, iflux
      read(inp,*)sdummy, ILIMIT
      read(inp,*)sdummy, vortex, xref, yref
      read(inp,*)sdummy, cell_type
      read(inp,*)sdummy, gridfile
      close(inp)

      inpstatus = yes

      if(istart .ne. scratch .and. istart .ne. restart)then
            print*,'Unknown start option',istart
            print*,'Possible values: 1=scratch or 2=restart'
            inpstatus = no
      endif

      if(iflow .ne. inviscid .and. iflow .ne. laminar .and.
     &   iflow .ne. turbulent)then
            print*,'Unknown flow type',iflow
            print*,'Possible values: 1=inviscid, 2=laminar, 3=turbulent'
            inpstatus = no
      endif

      if(iflux .ne. iroe .and. iflux .ne. ikfvs .and. 
     &   iflux .ne. ihllc)then
            print*,'Unknown flux',iflux
            print*,'Possible values: 1=roe, 2=kfvs'
            inpstatus = no
      endif

      if(ilimit .ne. no .and. ilimit .ne. yes)then
            print*,'Unknown limiter option',ilimit
            print*,'Possible values: 0=no, 1=yes'
            inpstatus = no
      endif

      if(vortex .ne. yes .and. vortex .ne. no)then
            print*,'Unknown vortex option',vortex
            print*,'Possible values: 0=no, 1=yes'
            inpstatus = no
      endif

      if(cell_type .ne. median .and. cell_type .ne. barth)then
            print*,'Unknown cell type',cell_type
            print*,'Possible values: 1=median cell, 2=barth cell'
            inpstatus = no
      endif

      if(inpstatus .eq. no) stop

      return
      end

C----------------------------------------------------------------------
C.....Calculate length used for time step
C.....For each point find the minimum altitude of all triangles
C.....surrounding that point
C----------------------------------------------------------------------
      subroutine dtlength(coord, elarea, elem, drmin)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), elarea(ntmax), drmin(npmax)

      integer          ip, it, n1, n2, n3
      double precision dx1, dx2, dx3, dy1, dy2, dy3, dr1, dr2, dr3, 
     &                 h1, h2, h3


      do ip=1,np
         drmin(ip) = 1.0d15
      enddo

      do it=1,nt
         n1  = elem(1,it)
         n2  = elem(2,it)
         n3  = elem(3,it)

         dx1 = coord(1,n2) - coord(1,n3)
         dy1 = coord(2,n2) - coord(2,n3)
         dr1 = dsqrt(dx1**2 + dy1**2)
         h1  = 2.0d0*elarea(it)/dr1

         dx2 = coord(1,n3) - coord(1,n1)
         dy2 = coord(2,n3) - coord(2,n1)
         dr2 = dsqrt(dx2**2 + dy2**2)
         h2  = 2.0d0*elarea(it)/dr2

         dx3 = coord(1,n1) - coord(1,n2)
         dy3 = coord(2,n1) - coord(2,n2)
         dr3 = dsqrt(dx3**2 + dy3**2)
         h3  = 2.0d0*elarea(it)/dr3

         drmin(n1) = dmin1( h1, drmin(n1) )
         drmin(n2) = dmin1( h2, drmin(n2) )
         drmin(n3) = dmin1( h3, drmin(n3) )
      enddo

      return
      end

C----------------------------------------------------------------------
C.....Save the current solution into another array
C----------------------------------------------------------------------
      subroutine save_old(prim, prim_old)
      implicit none
      include 'param.h'
      double precision prim(nvar,npmax), prim_old(nvar,npmax)

      integer          i, j

      do i=1,np
         do j=1,nvar
            prim_old(j,i) = prim(j,i)
         enddo
      enddo

      return
      end
