
C----------------------------------------------------------------------
C FUNCTION FAZA
C
C Called by: LAISUM
C
C Purpose: calculate phase. QE and QM calculate only the magnitude of the
C          collision function for positive omega. Here we multiply by -i and
C          take into account whether the collision function is real or
C          imaginary. We also multiply by the sign of omega if the collision
C          function is asymmetric wrt. omega.
C
C Formal parameters:
C      La     - lambda 1...6 = E1...6, 7,8 = M1,M2
C      Mi     - mu + 1
C      Rmu    - mu
C      Rsg    - sign of omega
C Return value:
C      Phase

      COMPLEX*16 FUNCTION FAZA(La,Mi,Rmu,Rsg)
      IMPLICIT NONE
      INTEGER*4 ieven , La , Mi
      REAL*8 Rmu , Rsg
      COMPLEX*16 ci
      DATA ci/(0.,1.)/ ! sqrt(-1)

      IF ( La.GT.6 ) THEN ! M1, M2 multipolarity
         FAZA = -ci
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         IF ( La.EQ.7 ) RETURN
         IF ( Mi.EQ.2 ) RETURN
         FAZA = DCMPLX(Rsg,0.D0)
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         GOTO 99999
      ELSE                ! E1...6
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            FAZA = -ci ! mu is even
            RETURN
         ENDIF
      ENDIF
      FAZA = DCMPLX(Rsg,0.D0) ! mu is odd, so sign changes with sign of omega
99999 END
