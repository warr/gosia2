
C----------------------------------------------------------------------
C FUNCTION TCEXP
C
C Called by: EXPON
C
C Purpose: evaluates a complex exponential
C
C Formal parameters:
C      Z      - argument of exponential (complex)
C
C Return value:
C      complex exponential of Z.

      COMPLEX*16 FUNCTION TCEXP(Z)
      IMPLICIT NONE
      REAL*8 a , b , c , d
      COMPLEX*16 Z

      a = DBLE(Z)
      b = DIMAG(Z)
      a = EXP(a)
      c = a*COS(b)
      d = a*SIN(b)
      TCEXP = DCMPLX(c,d)
      END
