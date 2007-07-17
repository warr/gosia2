 
C----------------------------------------------------------------------
C SUBROUTINE LAISUM
C
C Called by: AMPDER, STING
C Calls:     FAZA
C
C Purpose: evaluate the sum  over matrix elements.
C
C Uses global variables:
C      ARM    - reduced matrix elements
C      CAT    -
C      ELM    - matrix elements
C      EXPO   - adiabatic exponential
C      ISG    -
C      ISG1   -
C      ISHA   -
C      ISO    -
C      ISSTAR -
C      ISSTO  -
C      KDIV   -
C      LOCQ   - location of collision functions in ZETA array
C      LP7    - start of collision functions in ZETA (45100)
C      MSTORE -
C      NDIV   -
C      NSTART -
C      NSTOP  -
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ir     -
C      N      -
C      Rsg    -
C      Lam    - multipolarity
C      Ld     -
C      Nz     -
C      I57    - switch which is either 5 or 7.
C
C   \sum_{lmn} \zeta^{lm}_{kn} . M^(1)_{kn} f_{lm}(\omega) a_n(\omega)
C where
C   f_{lm} = -i Q_{lm} exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C and
C   M^(1)_{kn} = <k||E(M)\lambda||n>
C
C EXPO is exp(i * xi * sinh(w) + w) calculated in function EXPON.
C ARM is the reduced matrix element.
C q is the Qe or Qm calculated by the functions QE and QM, respectively and
C stored in ZETA array in the function SNAKE.
C z is the coupling parameter zeta, calculated in the function LSLOOP.
      
      SUBROUTINE LAISUM(Ir,N,Rsg,Lam,Ld,Nz,I57)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , CAT , D2W , DIPOL , ELM , ELML , ELMU , EN , 
     &       q , rmir , rmis , rmu , Rsg , SA , SPIN , z , ZETA , ZPOL
      INTEGER*4 i2 , i3 , I57 , iii , indq , indx , Ir , irs , is , 
     &          is1 , is2 , ISG , ISG1 , ISHA , ISMAX , ismin , ISO , 
     &          isplus , ISSTAR , ISSTO
      INTEGER*4 KDIV , la , Lam , LAMR , Ld , LOCQ , LP1 , LP10 , LP11 , 
     &          LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , 
     &          LP9 , LZETA
      INTEGER*4 m , mrange , MSTORE , mua , N , NDIV , NPT , NSTART , 
     &          NSTOP , NSW , Nz
      COMPLEX*16 ARM , FAZA , pamp , EXPO , pamp1
      COMMON /PSPIN / ISHA
      COMMON /AZ    / ARM(600,7)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      COMMON /ADBXI / EXPO(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      
      rmir = CAT(Ir,3)
      iii = 0
      IF ( Lam.GT.6 ) iii = 1
      la = Lam
      IF ( Lam.GT.6 ) Lam = Lam - 6
      DO i2 = 1 , Ld
         pamp = (0.,0.)
         m = MSTORE(1,i2)
         indx = MSTORE(2,i2)
         ismin = 0
         is1 = NSTART(m)
         IF ( is1.NE.0 ) THEN
            isplus = INT(rmir-CAT(is1,3)) - Lam
            IF ( isplus.LT.0 ) THEN
               ismin = isplus
               isplus = 0
            ENDIF
            is2 = is1 + isplus - 1
            mrange = 2*Lam + 1 + ismin
            IF ( is2+mrange.GT.NSTOP(m) ) mrange = NSTOP(m) - is2
            IF ( mrange.GT.0 ) THEN
               DO i3 = 1 , mrange
                  is = is2 + i3
                  rmis = CAT(is,3)
                  IF ( ISO.NE.0 .OR. rmir.LE..1 .OR. rmis.LE..1 ) THEN
                     rmu = rmis - rmir
                     mua = ABS(rmu) + 1.1
                     IF ( la.LE.6 .OR. mua.NE.1 ) THEN
                        indq = LOCQ(Lam,mua) + NPT
                        Nz = Nz + 1
                        z = ZETA(Nz)         ! Zeta
                        q = ZETA(indq+LP7)   ! Q-function
                        IF ( NDIV.NE.0 ) q = ZETA(indq+LP7) + DBLE(KDIV)
     &                       *(ZETA(indq+LP7+ISG1)-ZETA(indq+LP7))
     &                       /DBLE(NDIV)
                        pamp1 = FAZA(la,mua,rmu,Rsg)*q*z
                        IF ( ISO.NE.0 .OR. rmir.LE..1 ) THEN
                           pamp = pamp1*ARM(is,I57) + pamp
                           IF ( ISO.EQ.0 .AND. rmis.GT..1 ) GOTO 10
                        ENDIF
                        IF ( N.NE.m ) THEN
                           irs = (-1)**(INT(rmir+rmis)-ISHA+iii)
                           ARM(is,6) = ARM(is,6) + irs*pamp1*ARM(Ir,I57)
                           ISSTAR(i2) = MIN(is,ISSTAR(i2))
                           ISSTO(i2) = MAX(is,ISSTO(i2))
                        ENDIF
                     ENDIF
                  ENDIF
 10            ENDDO
               IF ( N.EQ.m ) THEN
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)
               ELSE
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)*EXPO(indx)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      Lam = la
      END