 
C----------------------------------------------------------------------
C SUBROUTINE INTG
C
C Called by: FTBM, GOSIA
C Calls:     AMPDER, DOUBLE, HALF, RESET
C
C Purpose: the main integration routine.
C
C Uses global variables:
C      ACC50  - accuracy required for integration
C      ARM    - reduced matrix elements
C      CAT    -
C      D2W    - step in omega (= 0.03)
C      IFLG   -
C      INTERV -
C      IPATH  -
C      IRA    - limit of omega for integration for each multipolarity
C      ISG    -
C      ISMAX  -
C      ISO    -
C      KDIV   -
C      LAMR   -
C      MAXLA  -
C      NDIV   -
C      NMAX   - number of levels
C      NPT    -
C      NSTART -
C      NSW    -
C
C Formal parameters:
C      Ien    - experiment number
C
C Note that if it finds that the step size for the integral is too small, it
C calls DOUBLE to increase it by a factor of two, or if it finds that the
C step size is too big, it decreases it by a factor of two by calling HALF.
C The function RESET may also be called to reset the step size.
C
C We use the Adams-Moulton predictor-corrector model.
 
      SUBROUTINE INTG(Ien)
      IMPLICIT NONE
      REAL*8 ACC50 , ACCA , ACCUR , CAT , D2W , DIPOL , EN , f , rim , 
     &       rl , SPIN , srt , ZPOL
      INTEGER*4 i , i57 , Ien , IFAC , IFLG , ihold , intend , INTERV , 
     &          IPATH , ir , ir1 , IRA , ISG , ISG1 , ISMAX , ISO , k , 
     &          kast , KDIV , LAMR
      INTEGER*4 MAGA , MAXLA , mir , n , NDIM , NDIV , NMAX , NMAX1 , 
     &          NPT , NSTART , NSTOP , NSW
      COMPLEX*16 ARM , hold
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /AZ    / ARM(600,7)
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /A50   / ACC50
      COMMON /CLCOM0/ IFAC(75)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /FLA   / IFLG
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /CEXC9 / INTERV(50)
      
      intend = INTERV(Ien)
      D2W = .03 ! We use steps of 0.03 in omega
      NSW = 1
      kast = 0
      NDIV = 0
      KDIV = 0
 100  IF ( (NPT+NSW).GT.IRA(MAXLA) .AND. ISG.GT.0 ) RETURN
      DO i = 1 , 8
         LAMR(i) = 0
         IF ( (NPT+NSW).LT.IRA(i) ) LAMR(i) = 1
      ENDDO
C     Predictor 
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX
            ir = NSTART(n) - 1
 120        ir = ir + 1
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
            mir = CAT(ir,3)
            ir1 = ir - 2*mir
            ARM(ir1,7) = IFAC(n)*ARM(ir,7)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 120
         ENDDO
      ELSE
         DO ir = 1 , ISMAX
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
         ENDDO
      ENDIF
      NPT = NPT + NSW*ISG
      IF ( NPT.GT.0 ) THEN
         IF ( NDIV.EQ.0 ) GOTO 200
         KDIV = KDIV + 1
         IF ( KDIV.LT.NDIV ) GOTO 200
         KDIV = 0
         NPT = NPT + ISG
         IF ( NPT.GT.0 ) GOTO 200
      ENDIF
      NPT = -NPT + 2
      ISG = 1
 200  CALL RESET(ISO)
      IFLG = 1
      i57 = 7

C     Calculate derivatives of amplitudes
      CALL AMPDER(i57)
      
C     Corrector
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX
            ir = NSTART(n) - 1
 220        ir = ir + 1
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
            mir = CAT(ir,3)
            ir1 = ir - 2*mir
            ARM(ir1,5) = IFAC(n)*ARM(ir,5)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 220
         ENDDO
      ELSE
         DO ir = 1 , ISMAX
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
         ENDDO
      ENDIF
      kast = kast + 1
      IFLG = 0
      i57 = 5

C     Calculate derivatives of amplitudes
      CALL AMPDER(i57)
      IF ( (LAMR(2)+LAMR(3)).NE.0 ) THEN
         IF ( kast.GE.intend ) THEN
            kast = 0
            f = 0.
            DO k = 1 , NMAX
               ihold = IPATH(k)
               IF ( ihold.NE.0 ) THEN
                  hold = ARM(ihold,5) - ARM(ihold,7)
                  rl = DBLE(hold)
                  rim = IMAG(hold)
                  srt = rl*rl + rim*rim
                  f = MAX(f,srt)
               ENDIF
            ENDDO

C           Decide if we have appropriate accuracy
            f = SQRT(f)/14.
            IF ( f.GT.ACCUR .OR. f.LT.ACC50 ) THEN
               IF ( f.LT.ACC50 ) THEN
                  CALL DOUBLE(ISO) ! Double step size
                  D2W = 2.*D2W
                  NSW = 2*NSW
                  intend = (DBLE(intend)+.01)/2.
                  IF ( intend.EQ.0 ) intend = 1
                  IF ( NSW.LT.1 ) THEN
                     NDIV = (DBLE(NDIV)+.01)/2.
                     IF ( NDIV.LT.2 ) THEN
                        NDIV = 0
                        NSW = 1
                     ENDIF
                  ENDIF
               ELSE
                  CALL HALF(ISO) ! Halve step size
                  D2W = D2W/2.
                  NSW = (DBLE(NSW)+.01)/2.
                  intend = 2*intend
                  IF ( NSW.LT.1 ) THEN
                     NDIV = 2*NDIV
                     IF ( NDIV.EQ.0 ) NDIV = 2
                  ENDIF
               ENDIF
            ENDIF
             
         ENDIF
      ENDIF
      GOTO 100
      END