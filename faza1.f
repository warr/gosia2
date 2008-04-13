      SUBROUTINE FAZA1(La,Mi,Rmir,Rmis,Dis,Rmu)
      IMPLICIT NONE
      INTEGER*4 ieven , irs , La , Mi
      REAL*8 Rmir , Rmis , Rmu
      COMPLEX*16 Dis , ci
      DATA ci/(0.,1.)/
      irs = (-1)**INT(Rmir+Rmis)
      IF ( La.EQ.7 ) THEN
         Dis = -ci*irs
         IF ( Rmu.LT.0. ) Dis = -Dis
         GOTO 99999
      ELSE
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            Dis = -ci*irs
            RETURN
         ENDIF
      ENDIF
      Dis = CMPLX(-DBLE(irs),0.)
99999 END
