      COMPLEX*16 FUNCTION TCEXP(Z)
      IMPLICIT NONE
      REAL*8 a , b , c , d
      COMPLEX*16 Z
      a = DBLE(Z)
      b = IMAG(Z)
      a = EXP(a)
      c = a*COS(b)
      d = a*SIN(b)
      TCEXP = CMPLX(c,d)
      END
