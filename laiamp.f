 
C----------------------------------------------------------------------
C SUBROUTINE LAIAMP
C
C Called by: STING
C Calls:     FAZA1, LEADF, STAMP, TCABS
C
C Purpose: calculate reduced matrix element
C
C Uses global variables:
C      ARM    - reduced matrix elements
C      CAT    - substates of levels (n_level, J, m)
C      ELM    - matrix elements
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 - 1)
C      IEXP   - number of experiment
C      IFAC   -
C      ISG    -
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LAMR   -
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LZETA  - index in ZETA to coupling coefficients for given multipolarity
C      NSTART -
C      NSTOP  -
C      XI     - xi coupling coefficients
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ir     - index into ARM array
C      W0     - omega
      
      SUBROUTINE LAIAMP(Ir,W0)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELML , ELMU , EPS , epsi , EROOT , errt , 
     &       FIEX , pm , ppp , rmir , rmis , rmu , SA , TCABS , W0 , 
     &       XI , xiv
      REAL*8 z , ZETA
      INTEGER*4 i1 , i2 , i3 , IAXS , IEXP , indx , Ir , is , is1 , 
     &          is2 , ISG , ISG1 , ISMAX , ismin , isplus , KDIV , la , 
     &          lam , LAMDA , LAMMAX
      INTEGER*4 LAMR , ld , LDNUM , LEAD , LEADF , LZETA , m , MEM , 
     &          mrange , mua , MULTI , NDIV , NPT , NSTART , NSTOP , 
     &          NSW , nz
      COMPLEX*16 ARM , STAMP , dis , uhuj
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /AZ    / ARM(600,7)
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CXI   / XI(500)

      ppp = 0.
      epsi = EPS(IEXP)
      errt = EROOT(IEXP)
      rmir = CAT(Ir,3) ! m quantum number of substate Ir
      DO i1 = 1 , LAMMAX
         lam = LAMDA(i1)
         nz = LZETA(lam)
         IF ( LAMR(lam).NE.0 ) THEN
            la = lam
            IF ( lam.GT.6 ) lam = lam - 6
            ld = LDNUM(la,1)
            IF ( ld.NE.0 ) THEN
               DO i2 = 1 , ld
                  m = LEADF(1,i2,la)
                  indx = MEM(1,m,la)
                  xiv = XI(indx)
                  ismin = 0
                  is1 = NSTART(m)
                  IF ( NSTART(m).NE.0 ) THEN
                     isplus = INT(rmir-CAT(is1,3)) - lam
                     IF ( isplus.LT.0 ) THEN
                        ismin = isplus
                        isplus = 0
                     ENDIF
                     is2 = is1 + isplus - 1
                     mrange = 2*lam + 1 + ismin
                     IF ( is2+mrange.GT.NSTOP(m) ) mrange = NSTOP(m)
     &                    - is2
                     IF ( mrange.GT.0 ) THEN
                        DO i3 = 1 , mrange
                           is = is2 + i3
                           nz = nz + 1
                           z = ZETA(nz)
                           rmis = CAT(is,3) ! m quantum number of substate is
                           rmu = rmis - rmir
                           mua = ABS(rmu) + 1.1
                           IF ( lam.LE.6 .OR. mua.NE.1 ) THEN
                              CALL FAZA1(la,mua,rmir,rmis,dis,rmu)
                              pm = ELM(indx)*z
                              uhuj = STAMP(epsi,errt,xiv,.03D0,W0,lam,
     &                               mua)
                              ARM(is,5) = dis*pm*uhuj
                              ppp = ppp + TCABS(ARM(is,5))
     &                              *TCABS(ARM(is,5))
                           ENDIF
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      ARM(Ir,5) = CMPLX(SQRT(1.-ppp),0.)
      END
