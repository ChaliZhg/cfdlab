C.....Calculate face length vector. Also save dual grid into a file for
C.....visualization
      subroutine flength(coord, tcoord, edge, edneigh, ds, dsb)
      implicit none
      include 'param.h'
      double precision coord(2,npmax), ds(2,nemax), dsb(2,npmax),
     &                 tcoord(2,ntmax)
      integer          edge(2,nemax), edneigh(2,nemax)

      integer          i, p1, p2, n1, n2, idual
      double precision x1, y1, x2, y2, xm, ym, nx, ny, flen

      idual = 20
      open(unit=idual, file='DUAL.DAT')

      do i=1,ne
         n1 = edge(1,i)
         n2 = edge(2,i)

         p1 = edneigh(1,i)
         p2 = edneigh(2,i)

         x1 = 0.0d0
         y1 = 0.0d0
         x2 = 0.0d0
         y2 = 0.0d0

c        Mid-point of the edge
         xm = 0.5d0*(coord(1,n1) + coord(1,n2))
         ym = 0.5d0*(coord(2,n1) + coord(2,n2))

         if(p1 .ne. 0)then
            x1 = tcoord(1,p1)
            y1 = tcoord(2,p1)
         else
            print*,'flength: Fatal error at edge',i
            print*,'         p1 is zero'
            stop
         endif

         if(p2 .ne. 0)then
            x2 = tcoord(1,p2)
            y2 = tcoord(2,p2)
         else
c           Edge number i is a boundary edge
            x2 = xm
            y2 = ym
         endif

         ds(1,i) = -( y2 - y1 )
         ds(2,i) =  ( x2 - x1 )

         write(idual,*)x1, y1
         write(idual,*)xm, ym
         write(idual,*)
         write(idual,*)x2, y2
         write(idual,*)xm, ym
         write(idual,*)
      enddo
      close(idual)

c     Boundary edges
      do i=1,np
         dsb(1,i) = 0.0d0
         dsb(2,i) = 0.0d0
      enddo

      do i=1,ne
         n1 =  edge(1,i)
         n2 =  edge(2,i)
         if(edneigh(1,i)*edneigh(2,i) .eq. 0)then
            nx =  coord(2,n2) - coord(2,n1)
            ny =-(coord(1,n2) - coord(1,n1))

            dsb(1,n1) = dsb(1,n1) + 0.5d0*nx
            dsb(2,n1) = dsb(2,n1) + 0.5d0*ny
            dsb(1,n2) = dsb(1,n2) + 0.5d0*nx
            dsb(2,n2) = dsb(2,n2) + 0.5d0*ny
         endif
      enddo

      if(cell_type .eq. barth) call reorder_edges(edge, ds, edneigh)

c     Find minimum/maximum face lengths
      maxflen = 0.0d0
      minflen = 1.0d8
      do i=1,ne
         flen    = dsqrt( ds(1,i)**2 + ds(2,i)**2 )
         maxflen = dmax1(maxflen, flen)   
         minflen = dmin1(minflen, flen)
      enddo
      print*,'\tMinimum cv length        =',minflen
      print*,'\tMaximum cv length        =',maxflen

      return
      end

C.....For barth cell, some edges may not make any contribution, eg. when
C.....both the triangles sharing this edge are right-angled. We remove such
C.....edges from the list
      subroutine reorder_edges(edge, ds, edneigh)
      implicit none
      include 'param.h'
      integer          edge(2,nemax), edneigh(2,nemax)
      double precision ds(2,nemax)

      integer          i, ecount
      double precision length

      print*,'Checking if any cell faces have zero length...'

      ecount = 0
      do i=1,ne
         length = dsqrt( ds(1,i)**2 + ds(2,i)**2 )
         if(length .ne. 0.0d0)then
            ecount            = ecount + 1
            edge(1,ecount)    = edge(1,i)
            edge(2,ecount)    = edge(2,i)
            edneigh(1,ecount) = edneigh(1,i)
            edneigh(2,ecount) = edneigh(2,i)
            ds(1,ecount)      = ds(1,i)
            ds(2,ecount)      = ds(2,i)
         endif
      enddo

      print*,'\tRemoved ',ne-ecount,' edges which have zero length'
      print*,'\tFinal numbers of edges   = ',ecount

      ne = ecount

      return
      end
