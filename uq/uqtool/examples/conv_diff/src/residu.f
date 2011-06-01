c     Finite volume residual
      subroutine residu(nc, q, res, xc, xv, dx, visc)
      implicit none
      include 'param.h'
      integer :: nc, visc
      real    :: q(nc), res(nc), xc(nc), xv(nc+1), dx(nc)

      integer :: i, sndOrder
      real :: f,qleft,qright,sl,sr,viscflux
      real :: minmod
      real :: grads(nc)

c     For second order scheme, set to 1
      sndOrder=1

      call calc_grad(nc,xc,xv,q,grads)

      ! Left flux of first cell 
      qleft = ql 
      qright= q(1)
      if(sndOrder.eq.1) then
         qright=qright-0.5*dx(1)*grads(1)
      endif
      call RoeFlux(qleft, qright, f)
      res(1) = res(1) - f

      ! All interior cells
      do i=2,nc
        qleft = q(i-1)
        qright =q(i)

         if(sndOrder.eq.1) then
            qleft =qleft +0.5*dx(i-1)*grads(i-1)
            qright=qright-0.5*dx(i)*grads(i)
         endif

         call RoeFlux(qleft, qright, f)
         res(i-1) = res(i-1) + f
         res(i)   = res(i)   - f
      enddo

      ! right flux of last cell 
      qleft  = q(nc)
      qright = qr
      if(sndOrder.eq.1) then
         qleft = qleft + 0.5*dx(nc)*grads(nc)
      endif
      call RoeFlux(qleft, qright, f)
      res(nc) = res(nc) + f

c     Add viscous contributions
      if(visc.eq.1)then
         viscflux = (q(1)-ql)/(xc(1)-xv(1))
         res(1) = res(1) + viscflux
         do i=2,nc
            viscflux = 2.0*(q(i) - q(i-1))/(dx(i) + dx(i-1))
            res(i-1) = res(i-1) - viscflux
            res(i  ) = res(i  ) + viscflux
         enddo
         viscflux =  (qr-q(nc))/(xv(nc+1)-xc(nc))
         res(nc) = res(nc) - viscflux
      endif

      return
      end

c     Roe flux function
      subroutine RoeFlux(ql, qr, f)
      implicit none
      real :: ql, qr, f

      real  :: a

      a = 0.5*abs(ql+qr)
      f = 0.5*(ql*ql*0.5 + qr*qr*0.5)
      f = f - 0.5*a*(qr-ql)

      return
      end

!     minmod function
      real function minmod(sl, sr)
      implicit none
      real :: sl, sr
      real :: c1,c2
      integer :: limiting

      limiting =0
      if(limiting.eq.0) then
         minmod=0.5*(sl+sr)
      else 
         if(sl*sr.le.0.0 .or. abs(sl+sr).lt.1.e-12) then
            minmod=0.0
         else 
            c1 =0.50*(sl+sr)
            c2 =2*sl*sr/(sl+sr)
            if(abs(c2).ge.abs(c1)) then
               minmod=c1
            else
               minmod=c2
            endif
         endif
      endif

      end
