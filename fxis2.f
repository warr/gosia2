 
C----------------------------------------------------------------------
C FUNCTION FXIS1
C
C Called by: NEWCAT
C
C Purpose: return -1 * sign(xi) for N = 2,3,5 and 6
C
C Uses global variables:
C      XI     - xi coupling coefficients
 
      REAL*8 FUNCTION FXIS2(I,N)
      IMPLICIT NONE
      INTEGER*4 I , N
      REAL*8 XI
      COMMON /CXI   / XI(500)

      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS2 = -SIGN(1.,XI(I))
         GOTO 99999
      ENDIF
      FXIS2 = 1.
99999 END
