 
C----------------------------------------------------------------------
C SUBROUTINE QM
C
C Called by: SNAKE
C
C Purpose: calculate Qm values
C
C Formal parameters:
C      C     - cosh(omega) + epsilon
C      D     - sqrt(epsilon^2 - a) * epsilon / C
C      B2    - B^2
C      B4    - B^4
C      Ert   - sqrt(epsilon^2 -1)
C      Lmda  - lambda
C      Cq    - array where the results are returned
C
C We used different formulae depending on lambda (see the table of magnetic
C collision functions in the gosia manual).
C
C Lmda = lambda + 6. i.e. Lmda = 7 is M1, Lmda = 8 is M2.
 
 
      SUBROUTINE QM(C,D,B2,B4,Ert,Lmda,Cq)
      IMPLICIT NONE
      REAL*8 B2 , B4 , C , Cq , D , Ert
      INTEGER*4 Lmda
      DIMENSION Cq(7)

      IF ( Lmda.EQ.8 ) THEN
         Cq(1) = -.9185586536*C*Ert/B4
         Cq(2) = -Cq(1)*D/C
         GOTO 99999
      ENDIF
      Cq(1) = -.3535533905*Ert/B2
99999 END
