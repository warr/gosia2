      SUBROUTINE LSLOOP(Ir,N,Nz,Ld,Lam,La,Ssqrt,Icg,Iexp)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , CAT , DIPOL , ELM , ELML , ELMU , EN , phz , 
     &       PSI , QAPR , rmir , rmis , SA , SPIN , Ssqrt , WTHREJ , 
     &       ZETA , ZPOL
      INTEGER*4 i2 , i3 , IAPR , Icg , Iexp , IFAC , iiex , indx , 
     &          inqa , inr , ins , IPATH , Ir , is , is1 , is2 , ISEX , 
     &          ISMAX , ismin , ISO
      INTEGER*4 isplus , jg1 , jg2 , jrmir , La , Lam , lam2 , Ld , 
     &          LEADF , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7
      INTEGER*4 LP8 , LP9 , LZETA , m , MAGA , MEM , mrange , mt , N , 
     &          NSTART , NSTOP , Nz
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /PCOM  / PSI(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CLCOM0/ IFAC(75)
      lam2 = 2*Lam
      inr = CAT(Ir,2)*2.
      rmir = CAT(Ir,3)
      jrmir = 2.*rmir
      DO i2 = 1 , Ld
         m = LEADF(N,i2,La)
         indx = MEM(N,m,La)
         IAPR(indx,1) = N
         IAPR(indx,2) = m
         ismin = 0
         ins = SPIN(m)*2.
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
                  IF ( ISO.NE.0 .OR. rmis.LE..1 .OR. rmir.LE..1 ) THEN
                     jg1 = -rmis*2.
                     jg2 = (rmis-rmir)*2.
                     IF ( Icg.NE.2 .OR. ABS(jg2).LE.2*MAGA(Iexp) ) THEN
                        IF ( La.LE.6 .OR. jg2.NE.0 ) THEN
                           Nz = Nz + 1
                           IF ( Nz.LE.LP7 ) THEN
                              iiex = (ins+jg1)/2
                              phz = (-1.0)**iiex
                              ZETA(Nz) = phz*PSI(indx)
     &                           *Ssqrt*WTHREJ(ins,lam2,inr,jg1,jg2,
     &                           jrmir)
                              IF ( Icg.NE.1 ) THEN
                                 mt = CAT(is,1)
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
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      END
