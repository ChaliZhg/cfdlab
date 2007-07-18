C.....Calculate element and control volume areas for median cell
      subroutine areas(coord, tcoord, elem, elarea, cvarea, cvareamc)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), elarea(ntmax), cvarea(npmax),
     &                 tcoord(2,ntmax), cvareamc(npmax)

      integer          i

      if(cell_type .eq. median)then
         call areas_mc(coord, tcoord, elem, elarea, cvarea, cvareamc)
      else
         call areas_bc(coord, tcoord, elem, elarea, cvarea, cvareamc)
      endif

      if(minelarea .le. 0.0d0)then
         print*,'Fatal: Element area is zero/negative'
         stop
      endif

      print*,'\tMinimum element area     =',minelarea
      print*,'\tMaximum element area     =',maxelarea

      maxcvarea = 0.0d0
      mincvarea = 1.0d8

      do i=1,np
         maxcvarea = dmax1(maxcvarea, cvarea(i))
         mincvarea = dmin1(mincvarea, cvarea(i))
      enddo

      if(mincvarea .le. 0.0d0)then
         print*,'Fatal: Control volume area is zero/negative'
         stop
      endif

      print*,'\tMinimum cv area          =',mincvarea
      print*,'\tMaximum cv area          =',maxcvarea

      return
      end

C.....Calculate element and control volume areas for median cell
      subroutine areas_mc(coord, tcoord, elem, elarea, cvarea, cvareamc)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), elarea(ntmax), cvarea(npmax),
     &                 tcoord(2,ntmax), cvareamc(npmax)

      double precision dx1, dy1, dx2, dy2, area3
      integer          i, n1, n2, n3

      print*,'Finding element and control volume areas for MEDIAN cell'

      do i=1,np
         cvareamc(i) = 0.0d0
      enddo

      maxelarea = 0.0d0
      minelarea = 1.0d8
      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)

c        Triangle area
         dx1= coord(1,n2) - coord(1,n1)
         dy1= coord(2,n2) - coord(2,n1)

         dx2= coord(1,n3) - coord(1,n1)
         dy2= coord(2,n3) - coord(2,n1)

         elarea(i) = 0.5d0*( dx1*dy2 - dx2*dy1 )
         maxelarea = dmax1(maxelarea, elarea(i))
         minelarea = dmin1(minelarea, elarea(i))

         area3        = elarea(i)/3.0d0
         cvareamc(n1) = cvareamc(n1) + area3
         cvareamc(n2) = cvareamc(n2) + area3
         cvareamc(n3) = cvareamc(n3) + area3

         call centroid(i, elem, coord, tcoord)

      enddo

      do i=1,np
         cvarea(i) = cvareamc(i)
      enddo

      return
      end

C.....Calculate element and control volume areas
      subroutine centroid(el, elem, coord, tcoord)
      implicit none
      include 'param.h'
      integer          el, elem(nvemax,ntmax)
      double precision coord(2,npmax), tcoord(2,ntmax)

      integer          n1, n2, n3


      n1           = elem(1,el)
      n2           = elem(2,el)
      n3           = elem(3,el)
      tcoord(1,el) = ( coord(1,n1) + coord(1,n2) + coord(1,n3) )/3.0d0
      tcoord(2,el) = ( coord(2,n1) + coord(2,n2) + coord(2,n3) )/3.0d0

      return
      end

C.....Calculate element and control volume areas for Barth cell
      subroutine areas_bc(coord, tcoord, elem, elarea, cvarea, cvareamc)
      implicit none
      include 'param.h'
      integer          elem(nvemax,ntmax)
      double precision coord(2,npmax), elarea(ntmax), cvarea(npmax),
     &                 tcoord(2,ntmax), cvareamc(npmax)

      double precision dx1, dy1, dx2, dy2, x(5), y(5), area, area3
      integer          i, j, k, n1, n2, n3

      print*,'Finding element and control volume areas for BARTH cell'

      do i=1,np
         cvarea(i)   = 0.0d0
         cvareamc(i) = 0.0d0
      enddo

      maxelarea = 0.0d0
      minelarea = 1.0d8
      do i=1,nt
         n1 = elem(1,i)
         n2 = elem(2,i)
         n3 = elem(3,i)

c        Triangle area
         dx1          = coord(1,n2) - coord(1,n1)
         dy1          = coord(2,n2) - coord(2,n1)

         dx2          = coord(1,n3) - coord(1,n1)
         dy2          = coord(2,n3) - coord(2,n1)

         elarea(i)    = 0.5d0*( dx1*dy2 - dx2*dy1 )
         maxelarea    = dmax1(maxelarea, elarea(i))
         minelarea    = dmin1(minelarea, elarea(i))

         area3        = elarea(i)/3.0d0
         cvareamc(n1) = cvareamc(n1) + area3
         cvareamc(n2) = cvareamc(n2) + area3
         cvareamc(n3) = cvareamc(n3) + area3

         call circumcenter(i, elem, coord, tcoord)

C        Control volume area
         do j=1,3
            n1 = elem(j,i)

            if(j .eq. 1)then
               n2 = elem(3,i)
            else
               n2 = elem(j-1,i)
            endif

            if(j .eq. 3)then
               n3 = elem(1,i)
            else
               n3 = elem(j+1,i)
            endif

            x(1) = coord(1,n1)
            y(1) = coord(2,n1)
            x(2) = 0.5d0*(coord(1,n1) + coord(1,n3))
            y(2) = 0.5d0*(coord(2,n1) + coord(2,n3))
            x(3) = tcoord(1,i)
            y(3) = tcoord(2,i)
            x(4) = 0.5d0*(coord(1,n1) + coord(1,n2))
            y(4) = 0.5d0*(coord(2,n1) + coord(2,n2))
            x(5) = x(1)
            y(5) = y(1)

            area = 0.0d0
            do k=1,4
               area = area + x(k)*y(k+1) - x(k+1)*y(k)
            enddo
            if(area .lt. 0.0d0)then
               print*,'Area is non-positive'
               print*,'\t Area         =', area
               print*,'\t Triangle     =',i
               stop 'areas'
            endif
            area       = 0.5d0*area
            cvarea(n1) = cvarea(n1) + area
         enddo
      enddo

      return
      end

C.....Calculate element and control volume areas
      subroutine circumcenter(el, elem, coord, tcoord)
      implicit none
      include 'param.h'
      integer          el, elem(nvemax,ntmax)
      double precision coord(2,npmax), tcoord(2,ntmax)

      double precision dx1, dy1, dx2, dy2, dx3, dy3, l1, l2, l3, 
     &                 beta1, beta2, beta3, det, b1, b2, xc, yc, 
     &                 x1, x2, x3, y1, y2, y3, stat
      integer          n1, n2, n3

      n1    = elem(1,el)
      n2    = elem(2,el)
      n3    = elem(3,el)

      x1    = coord(1,n1)
      y1    = coord(2,n1)
      x2    = coord(1,n2)
      y2    = coord(2,n2)
      x3    = coord(1,n3)
      y3    = coord(2,n3)

      dx1   = x2 - x3
      dy1   = y2 - y3
      l1    = dx1**2 + dy1**2

      dx2   = x3 - x1
      dy2   = y3 - y1
      l2    = dx2**2 + dy2**2

      dx3   = x1 - x2
      dy3   = y1 - y2
      l3    = dx3**2 + dy3**2

      beta1 = dmax1(0.0d0, l2 + l3 - l1)
      beta2 = dmax1(0.0d0, l3 + l1 - l2)
      beta3 = dmax1(0.0d0, l1 + l2 - l3)

C This fix is supposed to remove very small cv faces.
C I am not totally happy with this one.
c     if(beta1.lt.beta2/2.0d0 .and. beta1.lt.beta3/2.0d0) beta1=0.0d0
c     if(beta2.lt.beta3/2.0d0 .and. beta2.lt.beta1/2.0d0) beta2=0.0d0
c     if(beta3.lt.beta1/2.0d0 .and. beta3.lt.beta2/2.0d0) beta3=0.0d0

C Find circumcenter
      det   = (x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)
      b1    = 0.5d0*( (x2 - x1)*(x2 + x1) + (y2 - y1)*(y2 + y1) )
      b2    = 0.5d0*( (x3 - x1)*(x3 + x1) + (y3 - y1)*(y3 + y1) )
      xc    = ( (y3-y1)*b1 - (y2-y1)*b2)/det
      yc    = (-(x3-x1)*b1 + (x2-x1)*b2)/det

      stat= beta1*beta2*beta3

      if(stat .eq. 0.0d0)then
            if    (beta1 .eq. 0.0d0)then
                  tcoord(1,el) = 0.5d0*(x2 + x3)
                  tcoord(2,el) = 0.5d0*(y2 + y3)
            elseif(beta2 .eq. 0.0d0)then
                  tcoord(1,el) = 0.5d0*(x3 + x1)
                  tcoord(2,el) = 0.5d0*(y3 + y1)
            elseif(beta3 .eq. 0.0d0)then
                  tcoord(1,el) = 0.5d0*(x1 + x2)
                  tcoord(2,el) = 0.5d0*(y1 + y2)
            else
                  print*,'areas: Fatal error'
                  print*,'       Atleast one if statement above must'
                  print*,'       evaluate to be true'
                  stop
            endif
      else
            tcoord(1,el) = xc
            tcoord(2,el) = yc
      endif

      return
      end
