 
C----------------------------------------------------------------------
C FUNCTION FXIS1
C
C Called by: NEWCAT
C
C Purpose: return -1 * sign(xi) except for N = 2,3,5 and 6
C
C Uses global variables:
C      XI     - xi coupling coefficients
C
C Formal parameters:
C      I      - index into XI array
C      N      - 
C
C Return value:
C      sign of xi
      
      REAL*8 FUNCTION FXIS1(I,N)
      IMPLICIT NONE
      INTEGER*4 I , N
      REAL*8 XI
      COMMON /CXI   / XI(1500)

      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS1 = 1.
         GOTO 99999
      ENDIF
      FXIS1 = -SIGN(1.D0,XI(I))
99999 END
