 
C----------------------------------------------------------------------
C SUBROUTINE LSLOOP
C
C Called by: LOAD
C Calls:     CODE7, LEADF, WTHREJ
C
C Purpose: calculates the coupling parameter zeta and stores it in the
C          ZETA array starting at the beginning of this array (note that
C          this array has other things in it as well as zeta).
C
C Uses global variables:
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      IAPR   - index of initial and final levels for each matrix element
C      IFAC   - spin/parity phase factor
C      ISO    - isotropic flag
C      LP7    - maximum number of zeta coefficients (45100)
C      MAGA   - number of magnetic substates in approximate calculation
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      PSI    - psi coefficients
C      QAPR   - approximate Coulomb amplitudes
C      SPIN   - spin of level
C      ZETA   -
C
C Formal parameters:
C      Ir     - index of first substate of level
C      N      - index of level
C      Nz     - index into ZETA array for this multipolarity
C      Ld     - number of matrix elements with this multipolarity
C      Lam    - lambda
C      La     - 1...6 for E1...6 or 7,8 for M1,2
C      Ssqrt  - sqrt(2 * lambda + 1)
C      Icg    - (read only)
C      Iexp   - experiment number
C
C \zeta_{kn}^{(\lambda n)} = \sqrt{2 \lambda + 1} *
C                            (-1)^{I_n - M_n} *
C                            \three_j{I_n \lambda I_k -M_n \mu M_k} *
C                            \psi_{kn}
C
C For the evaluation of the 3-j symbol, ins = 2 I_n, lam2 = 2 \lambda,
C inr = 2 I_k, jg1 = -2 M_n, jg2 = 2 * \mu, jrmir = 2 * M_k. Note that the
C parameters to WTHREJ are all doubled, so that this routine can cope with
C half-integers.
 
      SUBROUTINE LSLOOP(Ir,N,Nz,Ld,Lam,La,Ssqrt,Icg,Iexp)
      IMPLICIT NONE
      REAL*8 ELM , ELML , ELMU , phz , 
     &       PSI , QAPR , rmir , rmis , SA , Ssqrt , WTHREJ
      INTEGER*4 i2 , i3 , IAPR , Icg , Iexp , IFAC , iiex , indx , 
     &          inqa , inr , ins , IPATH , Ir , is , is1 , is2 , ISEX , 
     &          ismin
      INTEGER*4 isplus , jg1 , jg2 , jrmir , La , Lam , lam2 , Ld , 
     &          LEADF , m , MAGA , MEM , mrange , mt , N , 
     &          NSTART , NSTOP , Nz
      INCLUDE 'coex.inc'
      COMMON /PCOM  / PSI(1500)
      INCLUDE 'ccoup.inc'
      INCLUDE 'clcom8.inc'
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      INCLUDE 'mgn.inc'
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      COMMON /CLCOM0/ IFAC(75)
      
      lam2 = 2*Lam
      inr = CAT(Ir,2)*2. ! 2 * Spin of substate Ir
      rmir = CAT(Ir,3)   ! m quantum number of substate Ir
      jrmir = 2.*rmir
      DO i2 = 1 , Ld
         m = LEADF(N,i2,La) ! Index of final level
         indx = MEM(N,m,La) ! Index of matrix element
         IAPR(indx,1) = N   ! Index of initial level
         IAPR(indx,2) = m   ! Index of final level
         ismin = 0
         ins = SPIN(m)*2.
         is1 = NSTART(m) ! Index of first substate of level m
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
                  rmis = CAT(is,3) ! m quantum number of substate is
                  IF ( ISO.NE.0 .OR. rmis.LE..1 .OR. rmir.LE..1 ) THEN
                     jg1 = -rmis*2.
                     jg2 = (rmis-rmir)*2.
                     IF ( Icg.NE.2 .OR. ABS(jg2).LE.2*MAGA(Iexp) ) THEN
                        IF ( La.LE.6 .OR. jg2.NE.0 ) THEN
                           Nz = Nz + 1
                           IF ( Nz.LE.LP7 ) THEN
                              iiex = (ins+jg1)/2
                              phz = (-1.0)**iiex
                              ZETA(Nz) = phz*PSI(indx) ! This is really zeta
     &                           *Ssqrt*WTHREJ(ins,lam2,inr,jg1,jg2,
     &                           jrmir)
                              IF ( Icg.NE.1 ) THEN
                                 mt = CAT(is,1) ! level number of substate is
                                 CALL CODE7(Ir,is,N,mt,inqa,indx)
                                 IF ( ABS(ELM(indx)).LT.1.E-6 )
     &                                ELM(indx) = 1.E-6
                                 IF ( inqa.NE.-1 ) THEN
                                    QAPR(indx,1,inqa) = ZETA(Nz)
     &                                 *ELM(indx)
                                    IF ( ISO.EQ.0 .AND. inqa.EQ.1 )
     &                                 QAPR(indx,1,7) = QAPR(indx,1,1)
     &                                 *IFAC(m)
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF ! If isotropic or rmis < 1 or rmir < 1
               ENDDO ! Loop on substates
            ENDIF ! If range of substates is greater than 0
         ENDIF ! If there are substates
      ENDDO ! Loop on matrix elements
      END
