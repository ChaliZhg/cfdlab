C=======================================================================
c Alternate logic to identify boundary faces given by Paul Fischer
C=======================================================================
      subroutine usrdat2

      include 'SIZE'
      include 'TOTAL'

      character*3 cbv
      integer e,f,fmid(6)
      real xmin, xmax, ymax

      common /scrns/ one(lx1,ly1,lz1,lelt)

      n = nx1*ny1*nz1*nelt
      call rone (one,n)

      ifield = 1
      call dssum(one,nx1,ny1,nz1)  ! Identify shared faces

      im = (1+nx1)/2               ! Face midpoints
      jm = (1+ny1)/2
      km = (1+nz1)/2

      fmid(4) =   1 + nx1*( jm-1) + nx1*ny1*( km-1)   !  x = -1
      fmid(2) = nx1 + nx1*( jm-1) + nx1*ny1*( km-1)   !  x =  1
      fmid(1) =  im + nx1*(  1-1) + nx1*ny1*( km-1)   !  y = -1
      fmid(3) =  im + nx1*(ny1-1) + nx1*ny1*( km-1)   !  y =  1
      fmid(5) =  im + nx1*( jm-1) + nx1*ny1*(  1-1)   !  z = -1
      fmid(6) =  im + nx1*( jm-1) + nx1*ny1*(nz1-1)   !  z =  1


      xmin = -5.0
      xmax = +10.0
      ymax =  3.0
      tol  =  1.e-5

      iin = 0
      iout= 0
      iwall=0
      isym = 0

      do e = 1,nelt      ! set boundary conditions
      do f = 1,2*ndim
         count = one(fmid(f),1,1,e)
         xf    = xm1(fmid(f),1,1,e)
         yf    = ym1(fmid(f),1,1,e)
         zf    = zm1(fmid(f),1,1,e)

         if (count.lt.1.1) then       ! This is boundary edge
            if(abs(xf-xmin).lt.tol)then
                cbc(f,e,1) = 'v  '   ! inlet
                iin = iin + 1
            else if(abs(xf-xmax).lt.tol)then
                cbc(f,e,1) = 'O  '   ! outlet
                iout = iout + 1
            else if(abs(yf-ymax).lt.tol)then
                cbc(f,e,1) = 'SYM'   ! top wall
                isym = isym + 1
            else
                cbc(f,e,1) = 'W  '   ! bottom wall
                iwall = iwall + 1
            endif
         endif
      enddo
      enddo

      print*,iin,iout,iwall,isym
      call exitt

      return
      end
