
C----------------------------------------------------------------------
C SUBROUTINE STING
C
C Called by: FTBM, GOSIA
C Calls:     LAIAMP, LAISUM, NEWLV
C
C Purpose: calculate and store excitation amplitudes,
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      ELM    - matrix elements
C      EXPO   - adiabatic exponential
C      IFLG   - flag to determine whether to calculate exponential (so we don't calculate twice)
C      IRA    - limit of omega for integration for each multipolarity
C      ISG    - sign of omega
C      ISMAX  - number of substates used
C      ISSTAR - first substate for given level
C      ISSTO  - last substate for given level
C      KDIV   - index of division
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LZETA  - index in ZETA to coupling coefficients for a given multipolarity
C      MAXLA  - multipolarity to calculate here
C      MSTORE - index of final level number and index of matrix element
C      NDIV   - number of divisions
C      NPT    - index in ADB array (this is omega / 0.03)
C
C Formal parameters:
C      Irld   - index into ARM array

      SUBROUTINE STING(Irld)
      IMPLICIT NONE
      REAL*8 rsg , w0
      INTEGER*4 i , i57 , ibg , iend , indx , Irld , is2 , j , j1 ,
     &          jj , lam , ld , maxh , mm , n , nz
      INCLUDE 'clcom.inc'
      INCLUDE 'az.inc'
      INCLUDE 'comme.inc'
      INCLUDE 'adbxi.inc'
      INCLUDE 'fla.inc'
      INCLUDE 'pint.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'caux.inc'
      INCLUDE 'clcom8.inc'
      INCLUDE 'rng.inc'

      maxh = MAXLA ! Save MAXLA, so we can restore it later
 100  ISG = -1
      n = 1
      rsg = -1.
      IFLG = 1
      w0 = IRA(MAXLA)*.03 + .03 ! Maximum omega to calculate for (steps of 0.03)

      DO j = 1 , ISMAX ! For substate used, zero ARM array
         DO jj = 1 , 6
            ARM(j,jj) = (0.,0.)
         ENDDO
      ENDDO
      ARM(Irld,5) = (1.,0.) ! Set y_n to 1

      DO j = 1 , 8
         LAMR(j) = 0 ! Initially mark that we shouldn't calculate any multipolarity
      ENDDO
      LAMR(MAXLA) = 1 ! Mark that we should calculate this multipolarity

      NPT = IRA(MAXLA) + 1 ! Number of omega values to calculate for this multipolarity

      IF ( MAXLA.EQ.7 .AND. IRA(2).NE.0 ) THEN ! Special case of M1
         LAMR(2) = 1
         NPT = NPT - 1
         w0 = w0 - .03
      ENDIF

      NDIV = 0
      KDIV = 0

      DO j = 1 , 4 ! Loop over terms for Adams-Moulton corrector-predictor
         NPT = NPT - 1
         DO j1 = 1 , LAMMAX ! Loop up to maximum multipolarity to calculate
            lam = LAMDA(j1) ! Get the multipolarity
            IF ( LAMR(lam).NE.0 ) THEN ! If this multipolarity should be calculated

C              Calculate and store exponentials in EXPO
               CALL NEWLV(n,ld,lam) ! n = 1, so for ground-state

               IF ( ld.NE.0 ) THEN ! If there are levels connected to g.s. by the right multipolarity
                  nz = LZETA(lam) ! Index into zeta array for this multipolarity
                  ld = LDNUM(lam,1) ! Number of levels connected to ground state for this multipolarity
                  i57 = 5 ! Use ARM(I,5) in LAISUM for excitation amplitudes
C                 Calculate sum over matrix elements
                  CALL LAISUM(Irld,n,rsg,lam,ld,nz,i57)
                  DO mm = 1 , ld ! Loop over levels
                     indx = MSTORE(2,mm) ! Index of matrix element in ELM
                     ibg = ISSTAR(mm) ! First substate for this level
                     iend = ISSTO(mm) ! Last substate for this level
                     DO is2 = ibg , iend ! Loop over substates
                        ARM(is2,4) = ARM(is2,4) + ARM(is2,6)*ELM(indx)
     &                               /EXPO(indx)
                        ARM(is2,6) = (0.,0.)
                     ENDDO ! Loop over substates
                  ENDDO ! Loop over matrix elements for ground state for this multipolarity
               ELSEIF ( j1.EQ.MAXLA ) THEN ! Else if it is the last multipolarity
                  IRA(MAXLA) = -IRA(MAXLA) ! Make IRA negative
                  DO jj = 1 , LAMMAX
                     lam = LAMDA(jj)
                     IF ( IRA(lam).GT.0 ) GOTO 105
                  ENDDO
 105              MAXLA = LAMDA(jj) ! Advance MAXLA to next multipolarity
                  GOTO 100 ! Back to start
               ENDIF

            ENDIF ! If we should calculate this multipolarity
         ENDDO ! Loop over multipolarities
         IF ( j.EQ.4 ) GOTO 200 ! We've set everything up, so finish

         DO i = 1 , ISMAX ! Shift terms up one
            ARM(i,j) = ARM(i,4)
            ARM(i,4) = (0.,0.)
         ENDDO

      ENDDO ! Loop over terms for Adams-Moulton corrector-predictor

C     Calculate amplitude
 200  CALL LAIAMP(Irld,w0)

      MAXLA = maxh ! Restore MAXLA

      DO jj = 1 , 8
         IRA(jj) = ABS(IRA(jj)) ! Make sure all the IRA are positive again
      ENDDO
      END
