
C----------------------------------------------------------------------
C SUBROUTINE FAZA1
C
C Called by: LAIAMP
C
C Purpose: calculate complex phase. Note it is only called for E1...6, M1 with
C          mu non zero and not at all for M2.
C
C Formal parameters:
C      La     - lambda 1...6 = E1...6, 7,8 = M1,M2
C      Mi     - mu + 1
C      Rmir   - m of first state
C      Rmis   - m of second state
C      Dis    - complex phase
C      Rmu    - mu

      SUBROUTINE FAZA1(La,Mi,Rmir,Rmis,Dis,Rmu)
      IMPLICIT NONE
      INTEGER*4 ieven , irs , La , Mi
      REAL*8 Rmir , Rmis , Rmu
      COMPLEX*16 Dis , ci
      DATA ci/(0.,1.)/ ! sqrt(-1)

      irs = (-1)**INT(Rmir+Rmis)
      IF ( La.EQ.7 ) THEN ! M1
         Dis = -ci*irs
         IF ( Rmu.LT.0. ) Dis = -Dis
         GOTO 99999
      ELSE ! E1...6
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            Dis = -ci*irs
            RETURN
         ENDIF
      ENDIF
      Dis = DCMPLX(-DBLE(irs),0.D0)
99999 END
