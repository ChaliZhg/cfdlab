c     Finite volume residual
      subroutine residu(nc, q, res, xc, xv, dx)
      implicit none
      include 'param.h'
      integer :: nc
      real    :: q(nc), res(nc), xc(nc), xv(nc+1), dx(nc)

      integer :: i, sndOrder
      real :: f,qleft,qright,sl,sr,viscflux
      real :: minmod
      real :: grads(nc)

      sndOrder=1

      call calc_grad(nc,xc,xv,q,grads)

      ! Left flux of first cell 
      qleft = ql 
      qright= ql

      call RoeFlux(qleft, qright, f)
      res(1) = res(1) - f

      viscflux = (q(1)-ql)/(xc(1)-xv(1))
      res(1) = res(1) + viscflux

      ! All cells

      do i=2,nc
        qleft = q(i-1)
        qright =q(i)

      if(sndOrder.eq.1) then
          qleft =qleft +0.5d0*dx(i-1)*grads(i-1)
          qright=qright-0.5d0*dx(i)*grads(i)
      endif

      call RoeFlux(qleft, qright, f)
      res(i-1) = res(i-1) + f
      res(i)   = res(i)   - f

      viscflux = 0.5d0*(grads(i-1)+grads(i))
      res(i-1) = res(i-1) - viscflux
      res(i  ) = res(i  ) + viscflux
      enddo

      ! right flux of last cell 

      qleft=qr
      qright = qr

      call RoeFlux(qleft, qright, f)
      res(nc) = res(nc) + f

      viscflux =  (qr-q(nc))/(xv(nc+1)-xc(nc))
      res(nc) = res(nc) - viscflux

      return
      end

c     Roe flux function
      subroutine RoeFlux(ql, qr, f)
      implicit none
      real :: ql, qr, f

      real  :: a

      f = 0.5*(ql*ql*0.5 + qr*qr*0.5)

      a = 0.5*abs(ql+qr)

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
         if(sl*sr.le.0.d0 .or. dabs(sl+sr).lt.1.e-12) then
            minmod=0.d0
         else 
            c1 =0.5d0*(sl+sr)
            c2 =2*sl*sr/(sl+sr)
            if(dabs(c2).ge.dabs(c1)) then
               minmod=c1
            else
               minmod=c2
            endif
         endif
      endif

      end

      subroutine calc_grad(nc,xc,xv,q,grads)
      implicit none
      include 'param.h'
      integer nc
      real :: xc(*), xv(*), q(*), grads(*)

      integer i

      call calc_grad_lsq(ql,q(1),q(2),xv(1),xc(1),xc(2),grads(1))
      !grads(1) =(q(1)-ql)/(xc(1)-xv(1))
      do i=2,nc-1
      call calc_grad_lsq(q(i-1),q(i),q(i+1),xc(i-1),xc(i),xc(i+1),
     1                   grads(i))
      enddo
      call calc_grad_lsq(q(nc-1),q(nc),qr,xc(nc-1),xc(nc),xv(nc+1),
     1                   grads(nc))
      !grads(nc) =(qr-q(nc))/(xv(nc+1)-xc(nc))

      end

      subroutine calc_grad_lsq(qm,q,qp,xm,x,xp,grads)
      implicit none
      real :: qm,q,qp,xm,x,xp,grads
      real :: w1,w2,epsil,du1,du2,dx1,dx2

      epsil = 1.e-12

      du1=q-qm
      du2=qp-q
      dx1=x-xm
      dx2=xp-x

      w1=1.d0
      w2=1.d0

      grads = (w1*dx1*du1+w2*dx2*du2)/(w1*dx1*dx1+w2*dx2*dx2)

      return
      end
