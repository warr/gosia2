 
C----------------------------------------------------------------------
C FUNCTION TACOS
C
C Called by: ARCCOS, CEGRY, COORD, GOSIA
C Calls:     TASIN
C
C Purpose: evaluate arccosine(x)
C
C Formal parameters:
C      X      - value for which we are to evaluate the arccosine
C
C Return value:
C      arccosine of X
C
C We use: arccos(x) = pi/2 - arcsin(x)
 
      REAL*8 FUNCTION TACOS(X)
      IMPLICIT NONE
      REAL*8 TASIN , X
      INCLUDE 'fconst.inc'
      
      TACOS = pi/2.D0 - TASIN(X)
      END
