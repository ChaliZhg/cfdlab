      SUBROUTINE TRIDAG(A,B,C,R,U,N,CODE)
c*****************************************************************
c Solves for a vector U of length N the tridiagonal linear set
c M U = R, where A, B and C are the three main diagonals of matrix
c M(N,N), the other terms are 0. R is the right side vector.
c*****************************************************************
      PARAMETER(NMAX=10000)
      REAL*8 BET,GAM(NMAX),A(N),B(N),C(N),R(N),U(N)
      INTEGER CODE

      if(n.gt.nmax)then
         code = 100
         return
      endif

      IF(B(1).EQ.0.D0) THEN
         CODE=1
         RETURN
      END IF

      BET=B(1)
      U(1)=R(1)/BET
      DO J=2,N
         GAM(J)=C(J-1)/BET
         BET=B(J)-A(J)*GAM(J)
         IF(BET.EQ.0.D0) THEN
            CODE=2
            RETURN
         END IF
         U(J)=(R(J)-A(J)*U(J-1))/BET
      END DO

      DO J=N-1,1,-1
         U(J)=U(J)-GAM(J+1)*U(J+1)
      END DO
  
      CODE=0
      RETURN
      END
