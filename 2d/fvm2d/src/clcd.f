C.....Calculate lift and drag coefficients
      subroutine clcd(spts, dsb, prim, mu, qx, qy)
      implicit none
      include 'param.h'
      integer          spts(nspmax)
      double precision dsb(2,npmax), prim(nvar,npmax), mu(npmax),
     &                 qx(nvar,npmax), qy(nvar,npmax)

      integer          i, ip
      double precision xf, yf, txx, txy, tyy

      xf = 0.0d0
      yf = 0.0d0
      do i=1,nsp
         ip = spts(i)
         txx= 2.0d0/3.0d0*mu(ip)*(2.0d0*qx(2,ip) - qy(3,ip))
         txy=             mu(ip)*(      qy(2,ip) + qx(3,ip))
         tyy= 2.0d0/3.0d0*mu(ip)*(2.0d0*qy(3,ip) - qx(2,ip))
         xf = xf + prim(4,ip)*dsb(1,ip) - 
     &             txx*dsb(1,ip) - txy*dsb(2,ip)
         yf = yf + prim(4,ip)*dsb(2,ip) - 
     &             txy*dsb(1,ip) - tyy*dsb(2,ip)
      enddo

      cl =(-dsin(aoa)*xf + dcos(aoa)*yf)/(0.5d0*r_inf*q_inf**2)
      cd =( dcos(aoa)*xf + dsin(aoa)*yf)/(0.5d0*r_inf*q_inf**2)

      return
      end
