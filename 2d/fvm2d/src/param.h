      integer no, yes
      parameter(no=0, yes=1)

      integer istart, scratch, restart
      parameter(scratch=1, restart=2)
      common/starttype/istart

c npmax = number of points
c ntmax = number of triangles (elements in general)
c nemax = number of edges
c nvemax= max number of vertices in an element
      integer npmax, ntmax, nemax, nspmax, nfpmax, nopmax, nbpmax, 
     &        nbemax, nvemax
      parameter(npmax = 24000,
     &          ntmax = 2*npmax,
     &          nemax = 3*npmax,
     &          nspmax= 2000,
     &          nfpmax= 2000,
     &          nopmax= 2000,
     &          nbpmax= 2000,
     &          nbemax= 2000,
     &          nvemax=3)

c nvar = number of variables, fixed at 5, which includes turbulent
c viscosity
      integer nvar
      parameter(nvar=5)

C Actual number of points, elements, edges and boundary edges
C     np = number of points
C     nt = number of triangles/elements
C     ne = number of edges
C     nsp= number of solid boundary points
C     nfp= number of farfield boundary points
C     nop= number of outflow boundary points
C     nbp= number of boundary points, must equal nsp+nfp+nop
      integer np,nt,ne,nsp,nfp,nop,nbp,nbe
      common/dims/np,nt,ne,nsp,nfp,nop,nbp,nbe

C Range of bounding box
      double precision xmin, xmax, ymin, ymax
      common/range/xmin, xmax, ymin, ymax

c Type of grid
c any other value implies hybrid grid
      character gridfile*32, inpfile*32
      common/files/gridfile,inpfile

c Maximum elements surrounding a point
      integer mesup
      parameter(mesup=10)

c Maximum points surrounding a point
      integer mpsup
      parameter(mpsup=10)

      double precision CFL, MINRESIDUE, dtglobal
      integer iter, ITERLAST, MAXITER, saveinterval
      common/itparam/CFL,MINRESIDUE,dtglobal,iter,ITERLAST,
     &               MAXITER,saveinterval

      double precision airk(3), birk(3)
      integer NIRK
      common/timeintegration/airk,birk,NIRK

C     Number of contours
      integer niso
      common/contour/niso

      double precision mach_inf, aoa, aoa_deg, q_inf, u_inf, v_inf, 
     &                 r_inf, p_inf, T_inf, T_infd, ent_inf, 
     &                 prim_inf(nvar), a_inf, H_inf
      common/inf/mach_inf,aoa,aoa_deg,q_inf,u_inf,v_inf,r_inf,p_inf,
     &           T_inf,T_infd,ent_inf,prim_inf,a_inf,H_inf

C     GAMMA = ratio of specific heats
C     GAMMA1= GAMMA - 1
C     GAS_CONST = gas constant, this can be set to 1.0
C     M_PI = value of pi
      double precision GAMMA, GAMMA1, GAS_CONST, M_PI
      common/const/GAMMA, GAMMA1, GAS_CONST, M_PI

      double precision minelarea, maxelarea, mincvarea, maxcvarea
      double precision minflen, maxflen
      common/minmaxarea/minelarea, maxelarea, mincvarea, maxcvarea,
     &                  minflen, maxflen

C Range of some variables
C rmin,rmax = density
C pmin,pmax = pressure
C mmin,mmax = mach number
      double precision rmin, rmax, umin, umax, vmin, vmax, pmin, pmax, 
     &                 mmin, mmax, emin, emax, nmin, nmax
      common/minmaxprim/rmin, rmax, umin, umax, vmin, vmax, pmin, pmax, 
     &                  mmin, mmax, emin, emax, nmin, nmax

      double precision res, res1, resi
      integer          iresi
      common/resparam/res,res1,resi,iresi

      double precision wd1(nspmax), wdmin, wdmax
      common/wdparam/wd1, wdmin, wdmax

C     Define tags for point types
      integer interior, solid, farfield, outflow
      parameter(interior=0)
      parameter(solid=3)
      parameter(outflow=4)
      parameter(farfield=5)

C Small number
      double precision EPSILON
      parameter(EPSILON=1.0d-16)

C Limiter factor for MUSCL
      integer          GRADTYPE, ILIMIT
      double precision LFACT, ALBADA11, ALBADA12, ALBADA21, ALBADA22
      common/lim/LFACT, ALBADA11, ALBADA12, ALBADA21, ALBADA22, 
     &           GRADTYPE, ILIMIT

      double precision cl, cd
      common/liftdrag/cl,cd

C Size of connectivity list; required by mayavi
      integer lsize
      common/maya/lsize

      integer iflux, iroe, ikfvs, ihllc
      parameter(iroe=1, ikfvs=2, ihllc=3)
      common/flux/iflux

C laminar paramters
C Rey = Reynolds number
C SCONST = Constant in Sutherland Law
      integer          iflow
      double precision Rey, prandtl, prandtl_turb, SCONST
      common/viscparam/Rey, prandtl, prandtl_turb, SCONST, iflow

      integer inviscid, laminar, turbulent
      parameter(inviscid=1, laminar=2, turbulent=3)

      integer xvstatus
      common/xv/xvstatus

      integer cell_type, median, barth
      parameter(median=1, barth=2)
      common/celltypes/cell_type

      double precision Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1, 
     &                 Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2
      common/samodel/Cb1, Cb2, sigma_sa, kolm, Cw1, Cw2, Cw3, Cv1,
     &               Cv2, Cv11, Cw31, Cw32, kolm2, Cb2Sig1, Cb2Sig2

      integer          vortex
      double precision xref, yref
      common/farfieldbc/xref, yref, vortex
