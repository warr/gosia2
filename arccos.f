 
C----------------------------------------------------------------------
C FUNCTION ARCOS
C
C Called by: GOSIA
C Calls:     TACOS
C
C Purpose: calculates an arccosine in a particular range
C
C Formal parameters:
C      A      - argument
C      F      - range
C      Pi     - Pi must be set to 3.14159... before calling ARCCOS
C
C Return value:
C      arccosine(A) within range of F
 
      REAL*8 FUNCTION ARCCOS(A,F,Pi)
      IMPLICIT NONE
      REAL*8 A , an , F , Pi , q , qa , qap , TACOS
      INTEGER*4 ie , j , k

      q = TACOS(A)
      qa = q
      qap = q
      IF ( q.LE.F ) THEN
         DO j = 1 , 20
            an = 2*j*Pi
            DO k = 1 , 2
               qap = qa
               ie = (-1)**k
               qa = an + ie*q
               IF ( qa.GT.F ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
 100  ARCCOS = qa
      IF ( (qa-F).GT.Pi/2. ) ARCCOS = qap
      END
