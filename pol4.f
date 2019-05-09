
C----------------------------------------------------------------------
C FUNCTION POL4
C
C Called by: TRINT
C
C Purpose: evaluate a 4th order polynomial
C
C Formal parameters:
C      C0     - coefficient of polynomial
C      C1     - coefficient of polynomial
C      C2     - coefficient of polynomial
C      C3     - coefficient of polynomial
C      C4     - coefficient of polynomial
C
C Evaluate C0 * A^4 + C1 * A^3 + C2 * A^2 + C3 * A + C4

      REAL*8 FUNCTION POL4(C0,C1,C2,C3,C4,A)
      IMPLICIT NONE
      REAL*8 A , C0 , C1 , C2 , C3 , C4

      POL4 = C4 + A*(C3+A*(C2+A*(C1+A*C0)))
      END
