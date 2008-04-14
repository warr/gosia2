      REAL*8 FUNCTION POL4(C0,C1,C2,C3,C4,A)
      IMPLICIT NONE
      REAL*8 A , C0 , C1 , C2 , C3 , C4
      POL4 = 1.
      IF ( ABS(A).GT.1.E+9 ) RETURN
      POL4 = C4 + A*(C3+A*(C2+A*(C1+A*C0)))
      END
