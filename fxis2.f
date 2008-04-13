      FUNCTION FXIS2(I,N)
      IMPLICIT NONE
      REAL*8 FXIS2 , XI
      INTEGER*4 I , N
      COMMON /CXI   / XI(500)
      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS2 = -SIGN(1.,XI(I))
         GOTO 99999
      ENDIF
      FXIS2 = 1.
99999 END
