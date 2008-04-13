      SUBROUTINE AMPDER(I57)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELML , ELMU , rsg , SA , ZETA
      INTEGER*4 i1 , I57 , ibg , iend , iflg , indx , ir , is2 , ISG , 
     &          ISG1 , ISMAX , ISSTAR , ISSTO , k , KDIV , lam , LAMDA , 
     &          LAMMAX , LAMR , lax
      INTEGER*4 ld , LDNUM , LEAD , LZETA , m , mm , MSTORE , MULTI , 
     &          n , NDIM , NDIV , nhold , NMAX , NMAX1 , NPT , NSTART , 
     &          NSTOP , NSW , nz
      COMPLEX*16 ARM , EXPO
      COMMON /AZ    / ARM(600,7)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      COMMON /ADBXI / EXPO(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      DO k = 1 , ISMAX
         ARM(k,6) = (0.,0.)
         ARM(k,4) = (0.,0.)
      ENDDO
      ISG1 = ISG
      IF ( NPT.EQ.1 ) ISG1 = ABS(ISG1)
      rsg = DBLE(ISG)
      DO i1 = 1 , LAMMAX
         lam = LAMDA(i1)
         lax = lam
         nz = LZETA(lam)
         IF ( LAMR(lam).NE.0 ) THEN
            iflg = 1
            nhold = 1
 20         CALL NEWLV(nhold,ld,lam)
            IF ( ld.EQ.0 ) THEN
 30            nhold = nhold + 1
               IF ( NSTART(nhold).EQ.0 ) GOTO 30
               GOTO 20
            ELSE
               ir = NSTART(nhold) - 1
 40            ir = ir + 1
               IF ( ir.LE.ISMAX ) THEN
                  n = CAT(ir,1)
                  IF ( n.NE.nhold ) THEN
                     DO mm = 1 , ld
                        m = MSTORE(1,mm)
                        IF ( m.NE.nhold ) THEN
                           indx = MSTORE(2,mm)
                           ibg = ISSTAR(mm)
                           iend = ISSTO(mm)
                           DO is2 = ibg , iend
                              ARM(is2,4) = ARM(is2,4) + ARM(is2,6)
     &                           *ELM(indx)/EXPO(indx)
                              ARM(is2,6) = (0.,0.)
                           ENDDO
                        ENDIF
                     ENDDO
 42                  CALL NEWLV(n,ld,lam)
                     IF ( ld.EQ.0 ) THEN
                        ir = ir + NSTOP(n) - NSTART(n) + 1
                        n = n + 1
                        IF ( n.GT.NMAX ) GOTO 100
                        GOTO 42
                     ELSE
                        nhold = n
                     ENDIF
                  ENDIF
                  CALL LAISUM(ir,n,rsg,lax,ld,nz,I57)
                  GOTO 40
               ENDIF
            ENDIF
         ENDIF
 100  ENDDO
      END
