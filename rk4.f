 
C----------------------------------------------------------------------
C FUNCTION RK4
C
C Called by: KONTUR
C
C Purpose:
C
C Formal parameters:
C      Y      - 
C      H      -
C      F      - array of three coefficients
C
C Return value:
 
      REAL*8 FUNCTION RK4(Y,H,F)
      IMPLICIT NONE
      REAL*8 F , H , Y
      DIMENSION F(3)

      RK4 = Y + H*(F(1)+4.*F(2)+F(3))/6.
      END
