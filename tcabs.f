 
C----------------------------------------------------------------------
C FUNCTION TCABS
C
C Called by: laiamp, pomnoz
C
C Purpose: evaluates the absolute value of a complex number
C
C Formal parameters:
C      Z      - argument for abs
C
C Return value:
C      absolute value of complex number Z.
 
      REAL*8 FUNCTION TCABS(Z)
      IMPLICIT NONE
      REAL*8 a , b
      COMPLEX*16 Z
      
      a = DBLE(Z)
      b = DIMAG(Z)
      IF ( ABS(a).LT.1.E-16 ) a = 0.
      IF ( ABS(b).LT.1.E-16 ) b = 0.
      TCABS = SQRT(a*a+b*b)
      END
