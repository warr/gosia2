 
C----------------------------------------------------------------------
C FUNCTION FUNC1
C
C Called by: LAGRAN, SPLNER
C
C Purpose: evaluates f(y) = y, f(y) = e^y or f(y) = y^2, depending on the
C flag I
C
C Formal parameters:
C      Y      - argument to evaluate
C      I      - mode: 1 = linear, 2 = exp, 3 = square
C
C Return value:
C      returns either y, exp(y) or y^2 depending on mode
C
C Note that this is the inverse of FUNC

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
