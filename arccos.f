 
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
C
C Return value:
C      arccosine(A) within range of F
 
      REAL*8 FUNCTION ARCCOS(A,F)
      IMPLICIT NONE
      REAL*8 A , an , F , q , qa , qap , TACOS
      INTEGER*4 ie , j , k
      INCLUDE 'fconst.inc'

      q = TACOS(A)
      qa = q
      qap = q
      IF ( q.LE.F ) THEN
         DO j = 1 , 20
            an = 2*j*pi
            DO k = 1 , 2
               qap = qa
               ie = (-1)**k
               qa = an + ie*q
               IF ( qa.GT.F ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
 100  ARCCOS = qa
      IF ( (qa-F).GT.pi/2.D0 ) ARCCOS = qap
      END
