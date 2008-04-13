      REAL*8 FUNCTION ARCTG(A,F,Pi)
      IMPLICIT NONE
      REAL*8 A , an , F , Pi , q , qa , qap
      INTEGER*4 ie , j , k
      q = ATAN(A)
      qa = q
      qap = q
      IF ( q.LE.F ) THEN
         DO j = 1 , 40
            an = j*Pi
            DO k = 1 , 2
               qap = qa
               ie = (-1)**k
               qa = an + ie*q
               IF ( qa.GT.F ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
 100  ARCTG = qa
      IF ( (qa-F).GT.Pi/4. ) ARCTG = qap
      END
