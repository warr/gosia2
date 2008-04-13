      SUBROUTINE DECAY(Chisq,Nlift,Chilo)
      IMPLICIT NONE
      REAL*8 AKS , bsum , Chilo , Chisq , DELLA , DELTA , df , DQ , 
     &       el1 , ELM , ELML , ELMU , emt , emt1 , ENDEC , ENZ , EPS , 
     &       EROOT , FIEX , FP
      REAL*8 gk , GKP , QCEN , SA , TAU , TIMEL , VACDP , vcd , XNOR , 
     &       ZETA
      INTEGER*4 i , IAXS , ibra , IBYP , idr , idrh , IEXP , ifn , il , 
     &          inx , inx1 , ITMA , iu , j , jlt , k , kl , KLEC , kq , 
     &          KSEQ
      INTEGER*4 l , l1 , lc1 , lc2 , LIFCT , LZETA , n1 , n2 , NDIM , 
     &          Nlift , NMAX , NMAX1
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /LIFE1 / LIFCT(50) , TIMEL(2,50)
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CATLF / FP(4,500,3) , GKP(4,500,2) , KLEC(75)
      COMMON /LCDL  / DELLA(500,3)
      DIMENSION gk(4)
      idr = 1
      DO il = 1 , NMAX1
         l = KSEQ(idr,3)
         n1 = 28*(l-1)
         ibra = KLEC(l)
         bsum = 0.
         idrh = idr
         DO j = 1 , ibra
            inx = KSEQ(idr,1)
            inx1 = KSEQ(idr,2)
            el1 = 0.
            IF ( inx.NE.0 ) el1 = ELM(inx)
            emt = el1*el1
            DELLA(idr,1) = emt
            IF ( inx1.NE.0 ) emt1 = ELM(inx1)*ELM(inx1)
            bsum = bsum + DELTA(idr,1)*emt
            IF ( inx1.NE.0 ) THEN
               DELLA(idr,3) = el1*ELM(inx1)
               DELLA(idr,2) = emt1
               bsum = bsum + DELTA(idr,2)*emt1
            ENDIF
            idr = idr + 1
         ENDDO
         idr = idrh
         TAU(l) = 1./bsum
         CALL GKVAC(l)
         DO j = 1 , ibra
            l1 = KSEQ(idr,4)
            n2 = 28*(l1-1)
            inx1 = KSEQ(idr,2)
            DO i = 1 , 4
               gk(i) = GKP(i,idr,1)*DELLA(idr,1)
            ENDDO
            IF ( inx1.NE.0 ) THEN
               DO i = 1 , 4
                  gk(i) = gk(i) + GKP(i,idr,2)*DELLA(idr,2)
               ENDDO
            ENDIF
            DO i = 1 , 4
               vcd = 1.
               IF ( i.NE.1 ) vcd = VACDP(i-1,l)
               gk(i) = gk(i)*TAU(l)
               ifn = 2*i - 1
               iu = (i-1)*7
               IF ( IAXS(IEXP).EQ.0 ) ifn = 1
               DO kq = 1 , ifn
                  lc1 = n1 + iu + kq
                  lc2 = n2 + iu + kq
                  ZETA(lc2) = ZETA(lc2) + gk(i)*vcd*ZETA(lc1)
               ENDDO
            ENDDO
            idr = idr + 1
         ENDDO
      ENDDO
      IBYP = 1
      IF ( Nlift.NE.0 .AND. IEXP.EQ.1 ) THEN
         DO jlt = 1 , Nlift
            kl = LIFCT(jlt)
            df = (TAU(kl)-TIMEL(1,jlt))/TIMEL(2,jlt)
            Chilo = Chilo + (LOG(TAU(kl)/TIMEL(1,jlt))*TIMEL(1,jlt)
     &              /TIMEL(2,jlt))**2
            Chisq = Chisq + df*df
         ENDDO
      ENDIF
      DO l = 2 , NMAX
         IF ( KLEC(l).NE.0 ) THEN
            n1 = 28*(l-1)
            DO j = 1 , 4
               vcd = 1.
               IF ( j.NE.1 ) vcd = VACDP(j-1,l)
               ifn = 2*j - 1
               iu = (j-1)*7
               DO k = 1 , ifn
                  lc1 = n1 + iu + k
                  ZETA(lc1) = ZETA(lc1)*vcd
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      END
