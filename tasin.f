      REAL*8 FUNCTION TASIN(X)
      IMPLICIT NONE
      REAL*8 dol , test , war , X
      test = ABS(X) - 1.
      IF ( ABS(test).LT.1.E-9 ) THEN
         TASIN = 1.570796327
         IF ( X.LT.0. ) TASIN = -1.570796327
         GOTO 99999
      ENDIF
      dol = SQRT(1.-X*X)
      war = X/dol
      TASIN = ATAN(war)
99999 END
