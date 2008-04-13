      FUNCTION FUNC1(Y,I)
      IMPLICIT NONE
      REAL*8 FUNC1 , Y
      INTEGER*4 I
      IF ( I.EQ.2 ) THEN
         FUNC1 = EXP(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC1 = Y*Y
         GOTO 99999
      ENDIF
      FUNC1 = Y
99999 END
