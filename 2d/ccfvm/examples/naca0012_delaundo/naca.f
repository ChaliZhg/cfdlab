      implicit double precision (a-h,o-z)
      dimension x(1000),y(1000),xo(1000),yo(1000),xw(1000)
      PRINT*,'Number of points on the profile'
      READ*,N
      PRINT*,'Number of points on the farfield boundary'
      READ*,M
      PRINT*,'Number of points on the wake'
      READ*,NW

      NoInnerBoundaries=1
      PI=4.0d0*DATAN(1.0d0)
      x(1)=1.00892d0
      y(1)=0.0d0
      DO i=1,N-1   !Finds points on the profile
        theta=(2.0d0*PI*(i-1))/(N-2)
        x(i+1)=0.5044d0*(dcos(theta)+1.0d0)
        y(i+1)=prof(x(i+1))
        IF(theta<PI)y(i+1)=-y(i+1)
      END DO

      DO i=1,N
        x(i)=x(i)/1.00892d0
        y(i)=y(i)/1.00892d0
      END DO

      router = 20.0d0
      DO i=1,M  !Calculates farfield points at a radius of 10 chord
         theta=(2.0d0*PI*(i-1))/M
         xo(i)=router*DCOS(theta)+0.5d0
         yo(i)=router*DSIN(theta)
      END DO

      iafl = 3
      ifar = 5

      print*,'Writing airfoil data into naca.pts'
      OPEN(1,FILE='naca.pts')
      WRITE(1,*)'NEWBND'
      WRITE(1,*)'NAMEBN'
      WRITE(1,*) iafl
      WRITE(1,*)'NFRSBN'
      WRITE(1,*) iafl
      WRITE(1,*)'NLSTBN'
      WRITE(1,*) iafl
      WRITE(1,*)'ITYPBN'
      WRITE(1,*) 1
      WRITE(1,*)'NRBNDE'
      WRITE(1,*) N-1-2
      WRITE(1,*)'BNDEXY'
      WRITE(1,5)x(1),y(1)
      !WRITE(1,5)(x(i),y(i),i=3,N-1)
      WRITE(1,5)(x(i),y(i),i=4,N-2)
      WRITE(1,5)x(1),y(1)

      WRITE(1,*)'NEWBND'
      WRITE(1,*)'NAMEBN'
      WRITE(1,*) ifar
      WRITE(1,*)'NFRSBN'
      WRITE(1,*) ifar
      WRITE(1,*)'NLSTBN'
      WRITE(1,*) ifar
      WRITE(1,*)'ITYPBN'
      WRITE(1,*) 2
      WRITE(1,*)'NRBNDE'
      WRITE(1,*) M+1
      WRITE(1,*)'BNDEXY'
      WRITE(1,5)(xo(i),yo(i),i=1,M)
      WRITE(1,5)xo(1),yo(1)

c Spacing on the wake
      if(NW .ne. 0)then
      dxw = x(1) - x(3)
      print*,'Spacing on wake = ',dxw
      xw(1) = x(1)
      fact = 2.00d0
      do i=1,NW-1
            dx      = 4.0d0*(20.0d0*dxw - dwx)*(i-1)*(NW-i)
            dx      = dx/NW/NW
            dx      = dxw + dx
            xw(i+1) = xw(i) + dx
      enddo
      if(xw(NW) .ge. xo(1))then
            print*,'Check the wake points'
      endif

      WRITE(1,*)'NEWBND'
      WRITE(1,*)'NAMEBN'
      WRITE(1,*) 3
      WRITE(1,*)'NFRSBN'
      WRITE(1,*) 1
      WRITE(1,*)'NLSTBN'
      WRITE(1,*) 0
      WRITE(1,*)'ITYPBN'
      WRITE(1,*) 4
      WRITE(1,*)'NRINDE'
      WRITE(1,*) NW
      WRITE(1,*)'BNDEXY'
      WRITE(1,5)(xw(i),0.0,i=1,NW)

      endif

2     FORMAT(A10)
3     FORMAT(I5)
4     FORMAT(3I5)
6     FORMAT(2I5)
5     FORMAT(2E18.7)

      write(1,*)'ENDDAT'
      CLOSE(1)

      stop
      end

C NACA00XX profile
      double precision function prof(x)
      implicit double precision (a-h,o-z)

      C1=0.2969d0
      C2=-0.1260d0
      C3=-0.3516d0
      C4=0.2843d0
      C5=-0.1015d0
      prof=5.0d0*0.12d0*(C1*DSQRT(x)+C2*x+C3*x**2+C4*x**3+C5*x**4)
c      RETURN prof
      end

