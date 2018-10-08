 
C----------------------------------------------------------------------
C FUNCTION TASIN
C
C Called by: CMLAB, COORD, TACOS
C
C Purpose: calculate an arcsine(x)
C
C Formal parameters:
C      X      - value for which we are to evaluate the arcsine
C
C Return value:
C      arcsine of X
C
C We take care of the special case of abs(x) = 1. Otherwise, we evaluate
C arctan(x / sqrt(1 - x^2).

      REAL*8 FUNCTION TASIN(X)
      IMPLICIT NONE
      REAL*8 dol , test , war , X
      INCLUDE 'fconst.inc'
      
      test = ABS(X) - 1.
      IF ( ABS(test).LT.1.D-9 ) THEN
         TASIN = pi/2.D0
         IF ( X.LT.0. ) TASIN = -pi/2.D0
         GOTO 99999
      ENDIF
      dol = SQRT(1.D0-X*X)
      war = X/dol
      TASIN = ATAN(war)
99999 END
