 
C----------------------------------------------------------------------
C FUNCTION LEADF
C
C Called by: LAIAMP, LSLOOP, NEWLV, SEQ
C
C Uses global variables:
C      MULTI  - number of matrix elements with a given multipolarity
C      IFAC   -
C
C Formal parameters:
C      N1     - number of levels
C      N2     -
C      N3     - multipolarity
C
C Purpose: calculate an index used for some arrays.
 
      INTEGER*4 FUNCTION LEADF(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMDA , LAMMAX , LDNUM , LEAD , lsum , MULTI , N1 , 
     &          n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)

      lsum = 0
      n3m = N3 - 1
      IF ( n3m.NE.0 ) THEN
         DO k = 1 , n3m
            lsum = lsum + MULTI(k)
         ENDDO
      ENDIF
      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m
            lsum = lsum + LDNUM(N3,k)
         ENDDO
      ENDIF
      n1m = lsum + N2
      LEADF = LEAD(2,n1m)
      END
