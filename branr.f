 
C----------------------------------------------------------------------
C SUBROUTINE BRANR
C
C Called by: FTBM
C Calls:     CONV
C
C Purpose: calculate the theoretical branching ratios and compare to the
C experimental ones.
C
C Uses global variables:
C      BRAT   - branching ratio and its error
C      DELTA  - \delta_\lambda: index 1 = electric^2, 2 = magnetic^2, 3 = cross term
C      ELM    - matrix elements
C      EN     - level energy
C      IBRC   - index branching ratios
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      MULTI  - number of matrix elements with given multipolarity
C      NBRA   - number of branching ratios
C
C Formal parameters:
C      Chisq  - chi squared
C      Nwyr   - number of data points contributing to chi squared
C      Chilo  - chi squared of logs

      SUBROUTINE BRANR(Chisq,Nwyr,Chilo)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , BRAT , ch1 , ch2 , Chilo , Chisq , CONV , 
     &       DELTA , DIPOL , ELM , ELML , ELMU , EN , ENDEC , eng1 , 
     &       eng2 , ENZ , SA , SPIN
      REAL*8 TAU , u , ZPOL
      INTEGER*4 i1 , i2 , IBRC , iflg , iout , IPRM , ISO , ITMA , itt , 
     &          j1 , j2 , k , KSEQ , lab1 , lab2 , LAMDA , LAMMAX , 
     &          LDNUM , LEAD , mul2
      INTEGER*4 MULTI , n1 , n2 , NBRA , Nwyr
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /PRT   / IPRM(20)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /LEV   / TAU(75) , KSEQ(500,4)

C     If no branching ratios were defined, return doing nothing
      IF ( NBRA.EQ.0 ) RETURN

C     If printing option is on, print something      
      IF ( IPRM(3).EQ.-1 ) WRITE (22,99001)
99001 FORMAT (1X,///10X,'EXP. AND CALCULATED BRANCHING RATIOS',//5X,
     &        'NS1',5X,'NF1',5X,'NS2',5X,'NF2',5X,'RATIO(1:2)',9X,
     &        'ERROR',7X,'CALC.RATIO',5X,'(EXP-CAL)/ERROR',//)

      Nwyr = Nwyr + NBRA ! Add to number of datapoints contributing to chisqr
      mul2 = MULTI(1) + MULTI(2)
      DO k = 1 , NBRA ! For each branching ratio
         ch1 = 0.
         ch2 = 0.
         iflg = 1
         itt = 1
         iout = 0
         n1 = IBRC(1,k) ! 1st matrix element
         n2 = IBRC(2,k) ! 2nd matrix element
         i1 = KSEQ(n1,1) ! Index of n1'th gamma
         i2 = KSEQ(n2,1) ! Index of n2'th gamma
         eng1 = EN(KSEQ(n1,3)) - EN(KSEQ(n1,4)) ! Energy of gamma 1
         eng2 = EN(KSEQ(n2,3)) - EN(KSEQ(n2,4)) ! Energy of gamma 2
         IF ( i1.NE.0 ) THEN
            IF ( i1.LE.MULTI(1) ) lab1 = 1
            IF ( i1.GT.MULTI(1) .AND. i1.LE.mul2 ) lab1 = 2
            IF ( i1.GT.mul2 ) lab1 = 3
         ENDIF
         IF ( i2.NE.0 ) THEN
            IF ( i2.LE.MULTI(1) ) lab2 = 1
            IF ( i2.GT.MULTI(1) .AND. i2.LE.mul2 ) lab2 = 2
            IF ( i2.GT.mul2 ) lab2 = 3
         ENDIF
         IF ( i1.NE.0 ) ch1 = ELM(i1)*ELM(i1)*DELTA(n1,1)
     &                        /(1.+CONV(eng1,lab1))
         IF ( i2.NE.0 ) ch2 = ELM(i2)*ELM(i2)*DELTA(n2,1)
     &                        /(1.+CONV(eng2,lab2))
         j1 = KSEQ(n1,2)
         IF ( j1.NE.0 ) THEN
            iflg = iflg + 1
            lab1 = lab1 + 2
            ch1 = ch1 + ELM(j1)*ELM(j1)*DELTA(n1,2)/(1.+CONV(eng1,lab1))
         ENDIF
         j2 = KSEQ(n2,2)
         IF ( j2.NE.0 ) THEN
            iflg = iflg + 1
            lab2 = lab2 + 2
            ch2 = ch2 + ELM(j2)*ELM(j2)*DELTA(n2,2)/(1.+CONV(eng2,lab2))
         ENDIF
         u = (ch1/ch2-BRAT(k,1))/BRAT(k,2)
         Chisq = Chisq + u*u
         Chilo = Chilo + (BRAT(k,1)*LOG(ch1/ch2/BRAT(k,1))/BRAT(k,2))**2
         IF ( IPRM(3).EQ.-1 ) WRITE (22,99002) KSEQ(n1,3) , KSEQ(n1,4) , 
     &                               KSEQ(n2,3) , KSEQ(n2,4) , BRAT(k,1)
     &                               , BRAT(k,2) , ch1/ch2 , -u
99002    FORMAT (5X,3(1I2,6X),1I2,5X,3(1F10.5,5X),5X,1F4.1)
      ENDDO ! Loop on branching ratios
       
      IF ( IPRM(3).EQ.-1 ) IPRM(3) = 0 ! Turn off printing option, so we don't print twice
      END
