      REAL*8 FUNCTION FUNC(Y,I)
      IMPLICIT NONE
      REAL*8 Y
      INTEGER*4 I
      IF ( I.EQ.2 ) THEN
         IF ( Y.LT.1.E-12 ) Y = 1.E-12
         FUNC = LOG(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC = SQRT(Y)
         GOTO 99999
      ENDIF
      FUNC = Y
99999 END
