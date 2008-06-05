 
C----------------------------------------------------------------------
C SUBROUTINE AMPDER
C
C Called by: INTG
C Calls:     LAISUM, NEWLV
C
C Purpose: to calculate the derivatives of the amplitudes needed for the
C Adams-Moulton predictor-corrector method.
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      EXPO   - exponents of adiabatic term
C      ISG    - phase
C      ISG1   - index of sigma
C      ISMAX  - number of substates used
C      ISSTAR - first substate for given level
C      ISSTO  - last substate for given level
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LAMR   - flag = 1 if we should calculate this multipolarity
C      LZETA  - index in ZETA to coupling coefficients for given multipolarity
C      MSTORE - index of final level number and index of matrix element
C      NMAX   - number of levels
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C
C Formal parameters:
C      I57    - switch which is either 5 or 7. This tells LAISUM to access either ARM(I,5) or ARM(I,7)

      SUBROUTINE AMPDER(I57)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELML , ELMU , rsg , SA , ZETA
      INTEGER*4 i1 , I57 , ibg , iend , iflg , indx , ir , is2 , ISG , 
     &          ISG1 , ISMAX , ISSTAR , ISSTO , k , KDIV , lam , 
     &          LAMR , lax
      INTEGER*4 ld , LZETA , m , mm , MSTORE , 
     &          n , NDIM , NDIV , nhold , NMAX , NMAX1 , NPT , NSTART , 
     &          NSTOP , NSW , nz
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(600,7)
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      INCLUDE 'clcom.inc'
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      INCLUDE 'adbxi.inc'
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)

C     Zero ARM(k,4) and ARM(k,6) for each substate used
      DO k = 1 , ISMAX ! ISMAX is number of substates used
         ARM(k,6) = (0.,0.)
         ARM(k,4) = (0.,0.)
      ENDDO

      ISG1 = ISG
      IF ( NPT.EQ.1 ) ISG1 = ABS(ISG1)
      rsg = DBLE(ISG)

      DO i1 = 1 , LAMMAX ! LAMMAX is number of multipolarities to calculate
         lam = LAMDA(i1) ! For each value of lambda, the user wants to calculate for
         lax = lam
         nz = LZETA(lam) ! Index into ZETA array for each multipolarity
         IF ( LAMR(lam).NE.0 ) THEN ! LAMR is flag to decide if we calculate for this multipolarity
            iflg = 1
            nhold = 1
 20         CALL NEWLV(nhold,ld,lam)
            IF ( ld.EQ.0 ) THEN ! If there are no decays
 30            nhold = nhold + 1
               IF ( NSTART(nhold).NE.0 ) GOTO 20
               GOTO 30
            ELSE
               ir = NSTART(nhold) - 1 ! Get first substate - 1 for this level
 40            ir = ir + 1 ! ir is a substate
               IF ( ir.LE.ISMAX ) THEN
                  n = CAT(ir,1) ! Level number of substate ir
                  IF ( n.NE.nhold ) THEN
                     DO mm = 1 , ld ! Loop over matrix elements
                        m = MSTORE(1,mm) ! Index of final level
                        IF ( m.NE.nhold ) THEN
                           indx = MSTORE(2,mm) ! Index of matrix element in ELM
                           ibg = ISSTAR(mm) ! First substate for this level
                           iend = ISSTO(mm) ! Last substate for this level
                           DO is2 = ibg , iend ! Loop over substates for level
                              ARM(is2,4) = ARM(is2,4) + ARM(is2,6)
     &                           *ELM(indx)/EXPO(indx)
                              ARM(is2,6) = (0.,0.)
                           ENDDO
                        ENDIF
                     ENDDO
 42                  CALL NEWLV(n,ld,lam)
                     IF ( ld.EQ.0 ) THEN ! if ld is zero, skip all the states for this level
                        ir = ir + NSTOP(n) - NSTART(n) + 1
                        n = n + 1
                        IF ( n.LE.NMAX ) GOTO 42
                        GOTO 100 ! IF this was the last level, loop back over lambda
                     ELSE
                        nhold = n
                     ENDIF
                  ENDIF ! If n .ne. nhold
                  CALL LAISUM(ir,n,rsg,lax,ld,nz,I57)
                  GOTO 40
               ENDIF ! If IR .le ISMAX
            ENDIF ! If LD .ne. 0
         ENDIF ! If LAMR(lam) .ne. 0
 100  ENDDO ! Loop over lambda
      END
