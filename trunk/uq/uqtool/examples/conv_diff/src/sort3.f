      SUBROUTINE SSORT (X, IY, N, KFLAG)
      IMPLICIT NONE
c
c    Example of an Insertion Sort
c
C***BEGIN PROLOGUE  SSORT
C***PURPOSE  Sort an array and make the same interchanges in
C            an auxiliary array.  The array is sorted in
C            decreasing order.
C***TYPE      SINGLE PRECISION
C***KEYWORDS  SORT, SORTING
C
C   Description of Parameters
C      X - array of values to be sorted   (usually abscissas)
C      IY - array to be carried with X (all swaps of X elements are
C          matched in IY .  After the sort IY(J) contains the original
C          postition of the value X(J) in the unsorted X array.
C      N - number of values in array X to be sorted
C      KFLAG - Not used in this implementation
C
C***REVISION HISTORY  (YYMMDD)
C   950310  DATE WRITTEN
C   John Mahaffy
C***END PROLOGUE  SSORT
C     .. Scalar Arguments ..
      INTEGER KFLAG, N
C     .. Array Arguments ..
      REAL X(*)
      INTEGER IY(*)
C     .. Local Scalars ..
      REAL TEMP
      INTEGER I, J, K, ITEMP
C     .. External Subroutines ..
C     None
C     .. Intrinsic Functions ..
C     None
C
C***FIRST EXECUTABLE STATEMENT  SSORT
C
      DO 100 I=2,N
         IF ( X(I).GT.X(I-1) ) THEN
            DO 50 J=I-2,1,-1
              IF(X(I).LT.X(J)) go to 70
  50          CONTINUE
            J=0
  70        TEMP=X(I)
            ITEMP=IY(I)
            DO 90 K=I,J+2,-1
              IY(K)=IY(K-1)
  90          X(K)=X(K-1)
            X(J+1)=TEMP
            IY(J+1)=ITEMP
         ENDIF
  100 CONTINUE
      RETURN
      END
