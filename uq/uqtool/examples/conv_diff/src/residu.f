c     Finite volume residual
      subroutine residu(nc, q, res, dx)
      implicit none
      include 'param.h'
      integer :: nc
      real :: q(nc), res(nc), dx

      integer :: i, sndOrder
      real :: f,qleft,qright,sl,sr,viscflux
      real :: minmod

      sndOrder=1

      ! Left flux of first cell 
      qleft = ql 
      qright= ql

      call RoeFlux(qleft, qright, f)
      res(1) = res(1) - f

      viscflux = -8.d0*ql+9.d0*q(1)-q(2)
      res(1) = res(1) + viscflux/dx

      ! Right flux of first cell 

      qleft = q(1)
      qright =q(2)

      if(sndOrder) then
          sl=2.0*(q(1)-ql)
          sr=q(2)-q(1)
          qleft =qleft+0.5*minmod(sl,sr)
          sl=q(2)-q(1)
          sr=q(3)-q(2)
          qright=qright-0.5*minmod(sl,sr)
      endif

      call RoeFlux(qleft, qright, f)
      res(1) = res(1) + f
      res(2) = res(2) - f

      viscflux = q(2)-q(1)
      res(1) = res(1) - viscflux/dx
      res(2) = res(2) + viscflux/dx
      ! All cells

      do i=3,nc-1
        qleft = q(i-1)
        qright =q(i)

      if(sndOrder) then
          sl=q(i-1)-q(i-2)
          sr=q(i)-q(i-1)
          qleft =qleft +0.5*minmod(sl,sr)
          sl=q(i)-q(i-1)
          sr=q(i+1)-q(i)
          qright=qright-0.5*minmod(sl,sr)
      endif

      call RoeFlux(qleft, qright, f)
      res(i-1) = res(i-1) + f
      res(i)   = res(i)   - f

      viscflux = q(i)-q(i-1)
      res(i-1) = res(i-1) - viscflux/dx
      res(i  ) = res(i  ) + viscflux/dx
      enddo

      ! left flux of last cell 

        qleft = q(nc-1)
        qright =q(nc)

      if(sndOrder) then
          sl=q(nc-1)-q(nc-2)
          sr=q(nc)-q(nc-1)
          qleft =qleft +0.5*minmod(sl,sr)
          sl=q(nc)-q(nc-1)
          sr=2.d0*(qr-q(nc))
          qright=qright-0.5*minmod(sl,sr)
      endif

      call RoeFlux(qleft, qright, f)
      res(nc-1) = res(nc-1) + f
      res(nc)   = res(nc)   - f

      viscflux = q(nc)-q(nc-1)
      res(nc-1) = res(nc-1) - viscflux/dx
      res(nc)   = res(nc)   + viscflux/dx

      ! right flux of last cell 

      qleft=qr
      qright = qr

      call RoeFlux(qleft, qright, f)
      res(nc) = res(nc) + f

      viscflux =  8.d0*qr-9.d0*q(nc)+q(nc-1)
      res(nc) = res(nc) - viscflux/dx

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
