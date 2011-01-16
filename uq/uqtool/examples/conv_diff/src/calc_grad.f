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

c     Least squares gradient
      subroutine calc_grad_lsq(qm,q,qp,xm,x,xp,grads)
      implicit none
      real :: qm,q,qp,xm,x,xp,grads
      real :: w1,w2,epsil,du1,du2,dx1,dx2

      du1=q-qm
      du2=qp-q
      dx1=x-xm
      dx2=xp-x

      w1=1.0
      w2=1.0

      grads = (w1*dx1*du1+w2*dx2*du2)/(w1*dx1*dx1+w2*dx2*dx2)

      return
      end
