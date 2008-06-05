 
C----------------------------------------------------------------------
C SUBROUTINE MIXR
C
C Called by: FTBM, GOSIA
C
C Purpose: calculate theoretical mixing ratio and compare to experimental one.
C
C Uses global variables:
C      DMIX   - 0.8326 * gamma energy
C      DMIXE  - mixing ratio and its error
C      IMIX   - decay associated with known mixing ratio
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LNY    - use logs to calculate chi squared
C      NDL    - number of mixing ratios
C
C Formal parameters:
C      Nw     - number of data points used to calculate chi squared
C      Ipsw   - printing flag (0 means no print, 1 means print)
C      Chi    - chi squared
C      Chilo  - chi squared using logs
 
      SUBROUTINE MIXR(Nw,Ipsw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , dl , DMIX , DMIXE
      INTEGER*4 i , IMIX , INTR , inx , inx1 , IPS1 , Ipsw , it , 
     &          LNY , NDL , Nw
      INCLUDE 'lev.inc'
      INCLUDE 'comme.inc'
      COMMON /MIXD  / DMIXE(20,2) , DMIX(20) , IMIX(20) , NDL
      COMMON /LOGY  / LNY , INTR , IPS1

      IF ( NDL.EQ.0 ) RETURN
      Nw = Nw + NDL

      DO i = 1 , NDL ! For each mixing ratio
         it = IMIX(i) ! Decay for this mixing ratio
         inx = KSEQ(it,1) ! Index 1 of it'th decay
         inx1 = KSEQ(it,2) ! Index 2 of it'th decay
         IF ( ABS(ELM(inx1)).LT.1.E-5 ) ELM(inx1) = 1.E-5
         dl = DMIX(i)*ELM(inx)/ELM(inx1)
         IF ( Ipsw.EQ.1 ) DMIX(i) = dl
         Chi = Chi + (dl-DMIXE(i,1))**2/DMIXE(i,2)/DMIXE(i,2)
         IF ( LNY.EQ.1 ) Chilo = Chilo + 
     &                           (DMIXE(i,1)*LOG(ABS(dl/DMIXE(i,1)))
     &                           /DMIXE(i,2))**2
      ENDDO ! Loop on mixing ratios i

      IF ( Ipsw.EQ.0 ) RETURN
      
      WRITE (22,99001)
99001 FORMAT (1X//10X,'E2/M1 MIXING RATIOS'/10X,'TRANSITION',10X,
     &        'EXP.DELTA',10X,'CALC.DELTA',10X,'SIGMA'/)

      DO i = 1 , NDL ! For each mixing ratio
         dl = (DMIX(i)-DMIXE(i,1))/DMIXE(i,2) ! Relative error
         it = IMIX(i) ! Matrix element for this mixing ratio
         WRITE (22,99002) KSEQ(it,3) , KSEQ(it,4) , DMIXE(i,1) , DMIX(i)
     &                    , dl ! KSEQs are level numbers
99002    FORMAT (10X,1I2,'---',1I2,14X,1F7.2,12X,1F7.2,13X,1F5.2)
      ENDDO ! Loop on mixing ratios i

      END
