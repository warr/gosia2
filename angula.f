 
C----------------------------------------------------------------------
C SUBROUTINE ANGULA
C
C Called by: GOSIA, CEGRY
C Calls:     FIINT, FIINT1, RECOIL, YLM, YLM1
C
C Purpose: calculate angular distribution of emitted gamma rays
C
C Uses global variables:
C      BETAR  - recoil beta
C      DELLA  -
C      DSIGS  -
C      ENDEC  -
C      ENZ    -
C      FP     -
C      IAXS   - axial symmetry flag
C      IEXP   - experiment number
C      ITMA   - identify detectors according to OP,GDET
C      ITTE   - thick target experiment flag
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      TAU    -
C      TLBDG  - theta of particle detector
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ygn    -
C      Idr    - number of decays
C      Iful   -
C      Fi0    - phi_0
C      Fi1    - phi_1
C      Trec   -
C      Gth    -
C      Figl   -
C      Ngl    - detector number
C      Op2    - option (either 'INTG' or not).
      
      SUBROUTINE ANGULA(Ygn,Idr,Iful,Fi0,Fi1,Trec,Gth,Figl,Ngl,Op2)
      IMPLICIT NONE
      REAL*8 AGELI , alab , arg , at , attl , BETAR , bt , CC , DELLA , 
     &       DELTA , EG , ENDEC , ENZ , EPS , EROOT , f , Fi0 , fi01 , 
     &       Fi1 , fi11
      REAL*8 FIEX , Figl , FP , GKP , Gth , Q , qv , sm , TAU , Trec , 
     &       Ygn , ylmr , ZETA
      REAL*8 ttx , dsig , TETACM , TREP , DSIGS , XA , XA1 , EP , ! For gosia2
     &       TLBDG, VINF ! For gosia2
      INTEGER*4 IAXS , Idr , IEXP , ifn , Iful , ig , il , inat , inx1 , 
     &          ipd , is , ITMA , ITTE , iu , ixs , j , ji , jj , jm , k
      INTEGER*4 KLEC , kq , KSEQ , l , lf , lf1 , LZETA , mind , NANG , 
     &          Ngl , NICC , nlv
      INTEGER*4 NEXPT, IZ, IZ1 ! For gosia2
      CHARACTER*4 Op2 ! For gosia2
      DIMENSION f(4) , ylmr(9,9) , at(28) , alab(9,9) , attl(9,9) , 
     &          Ygn(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /LCDL  / DELLA(500,3)
      COMMON /CATLF / FP(4,500,3) , GKP(4,500,2) , KLEC(75)
      COMMON /BREC  / BETAR(50)
      COMMON /THTAR / ITTE(50)
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) ,
     &                TLBDG(50) , VINF(50)
      
      DO l = 1 , Idr
         nlv = KSEQ(l,3) ! Level number of l'th decay
         il = (nlv-1)*28
         inx1 = KSEQ(l,2) ! Index of l'th decay
         DO j = 1 , 4
            f(j) = FP(j,l,1)*DELLA(l,1)
         ENDDO
         IF ( inx1.NE.0 ) THEN
            DO j = 1 , 4
               f(j) = f(j) + 2.*FP(j,l,3)*DELLA(l,3) + FP(j,l,2)
     &                *DELLA(l,2)
            ENDDO
         ENDIF
         DO j = 1 , 4
            f(j) = f(j)*TAU(nlv)
            iu = (j-1)*7
            ifn = 2*j - 1
            IF ( IAXS(IEXP).EQ.0 ) ifn = 1
            DO kq = 1 , ifn
               is = iu + kq
               ig = is + il
               at(is) = ZETA(ig)*f(j)
            ENDDO
         ENDDO
         IF ( Iful.EQ.1 ) THEN
            DO j = 1 , 9
               DO k = 1 , 9
                  alab(j,k) = 0.
                  attl(j,k) = 0.
               ENDDO
            ENDDO
            DO j = 1 , 4
               lf = 2*j - 1
               lf1 = lf
               IF ( IAXS(IEXP).EQ.0 ) lf1 = 1
               DO k = 1 , lf1
                  inat = (j-1)*7 + k
                  alab(lf,k) = at(inat)
               ENDDO
            ENDDO
            bt = BETAR(IEXP)
            IF ( ITTE(IEXP).NE.1 ) CALL RECOIL(alab,attl,bt,Trec)
            IF ( l.EQ.1 ) CALL YLM1(Gth,ylmr)
            ixs = IAXS(IEXP)
            fi01 = Fi0 - Figl
            fi11 = Fi1 - Figl
            CALL FIINT1(fi01,fi11,alab,ixs)
            Ygn(l) = alab(1,1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
            DO j = 2 , 9
               sm = ylmr(j,1)*alab(j,1)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  DO k = 2 , j
                     sm = sm + 2.*ylmr(j,k)*alab(j,k)
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl)
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,j-1)*Q(2,ipd,j-1)+Q(1,ipd,j-1)*arg)
     &              /(Q(2,ipd,j-1)+arg)
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ELSE
            ixs = IAXS(IEXP)
            fi01 = Fi0 - Figl
            fi11 = Fi1 - Figl
            CALL FIINT(fi01,fi11,at,ixs)
            IF ( l.EQ.1 ) CALL YLM(Gth,ylmr)
            Ygn(l) = at(1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
            DO jj = 1 , 3
               ji = jj*7 + 1
               sm = ylmr(jj,1)*at(ji)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  mind = 2*jj + 1
                  DO jm = 2 , mind
                     ji = ji + 1
                     sm = ylmr(jj,jm)*at(ji)*2. + sm
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl) ! Detector ID
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,2*jj)*Q(2,ipd,2*jj)+Q(1,ipd,2*jj)*arg)
     &              /(Q(2,ipd,2*jj)+arg) ! solid angle attenuation coefficients
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ENDIF
      ENDDO
      
      IF ( Op2.EQ.'INTG') RETURN

      dsig = DSIGS(IEXP)
      ttx = TLBDG(IEXP) / 57.2957795
      DO j = 1, Idr
         YGN(j) = YGN(j) * dsig * SIN(ttx)
      ENDDO
      END
