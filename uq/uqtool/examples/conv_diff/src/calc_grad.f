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
      real :: qm,q,qp,xm,x,xp,grads,qxp,qxm
      real :: w1,w2,epsil,du1,du2,dx1,dx2,xvm,xvp,qmax,qmin,phi

      epsil = 1.e-12

      du1=q-qm
      du2=qp-q
      dx1=x-xm
      dx2=xp-x

       w1=1.d0
       w2=1.d0

      grads = (w1*dx1*du1+w2*dx2*du2)/(w1*dx1*dx1+w2*dx2*dx2)

      xvm = 0.5*(x+xm)
      xvp = 0.5*(x+xp)

      qmax=max(max(q,qm),qp)
      qmin=min(min(q,qm),qp)

      qxp = q+(xvp-x)*grads
      qxm = q+(xvm-x)*grads

      phi = 1.d0

      if((qxp-q).gt.0.) phi = min(phi,(qmax-q)/(qxp-q))

      if((qxp-q).lt.0.) phi = min(phi,(qmin-q)/(qxp-q))

      if((qxm-q).gt.0.) phi = min(phi,(qmax-q)/(qxm-q))

      if((qxm-q).lt.0.) phi = min(phi,(qmin-q)/(qxm-q))


       !grads = grads*phi

      return
      end
