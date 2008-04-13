      REAL*8 FUNCTION TCABS(Z)
      IMPLICIT NONE
      REAL*8 a , b
      COMPLEX*16 Z
      a = DBLE(Z)
      b = IMAG(Z)
      IF ( ABS(a).LT.1.E-16 ) a = 0.
      IF ( ABS(b).LT.1.E-16 ) b = 0.
      TCABS = SQRT(a*a+b*b)
      END
