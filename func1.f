      REAL*8 FUNCTION FUNC1(Y,I)
      IMPLICIT NONE
      INTEGER*4 I
      REAL*8 Y
      IF ( I.EQ.2 ) THEN
         FUNC1 = EXP(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC1 = Y*Y
         GOTO 99999
      ENDIF
      FUNC1 = Y
99999 END
