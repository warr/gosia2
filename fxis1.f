      REAL*8 FUNCTION FXIS1(I,N)
      IMPLICIT NONE
      INTEGER*4 I , N
      REAL*8 XI
      COMMON /CXI   / XI(500)
      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS1 = 1.
         GOTO 99999
      ENDIF
      FXIS1 = -SIGN(1.,XI(I))
99999 END
