 
C----------------------------------------------------------------------
C SUBROUTINE CEGRY
C
C Called by: FTBM
C Calls:     ANGULA, DECAY, EFFIX, SIXEL, TACOS
C
C Purpose: calculate the gamma-ray deexcitation.
C
C Uses global variables:
C      AGELI  - angles of the Ge detectors
C      BETAR  - recoil beta
C      CNOR   - normalization factors
C      CNOR1  - normalization factors (added for gosia2)
C      CORF   - internal correction factors
C      DEV    -
C      DYEX   - error on experimental yield
C      EMH    -
C      ENDEC  -
C      FIEX   - phi range of particle detector
C      ICLUST -
C      IDRN   -
C      IEXP   - number of experiment
C      IFMO
C      IGRD   -
C      ILE    -
C      IMIN   -
C      INM    -
C      INNR   -
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      IRAWEX -
C      ITMA   - identify detectors according to OP,GDET
C      ITS    -
C      IWF    -
C      IY     - index for yields
C      JSKIP  -
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      KVAR   -
C      LASTCL -
C      LFL    -
C      LNORM  - normalization constant control
C      LP2    - maximum number of matrix elements (500)
C      LP6    - 32
C      LP10   - 600
C      MCFIX  - fixing parameter (added for gosia2)
C      NANG   - number of gamma-ray detectors for each experiment
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NLIFT  - number of lifetimes
C      NMAX   - number of levels
C      NYLDE  - number of yields
C      ODL    - distance from target to front face of detector
C      PART   -
C      PARTL  -
C      SPIN   - spin of level
C      SUBCH1 -
C      SUBCH2 -
C      SUMCL  -
C      TAU    -
C      TREP   -
C      UPL    - upper limits for all gamma detectors
C      VACDP  -
C      YEXP   - experimental yield
C      YGN    -
C      YGP    -
C      YNRM   - relative normalization factors for gamma detectors
C
C Formal parameters:
C      Chisq  - chi squared
C      Itemp  -
C      Chilo  - chi squared of logs
C      Idr    - number of decays
C      Nwyr   - number of data points contributing to chi squared
C      Icall  -
C      Issp   -
C      Iredv  -
 
      SUBROUTINE CEGRY(Chisq,Itemp,Chilo,Idr,Nwyr,Icall,Issp,Iredv)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AGELI , AKS , BETAR , CC , ccc , ccd , 
     &       Chilo , Chisq , CNOR , cnr , cocos , CORF , d , decen , 
     &       DELTA , DEV , DIPOL , DIX
      REAL*8 dl , DQ , DSIGS , DYEX , effi , EG , EMH , EN , ENDEC , 
     &       ENZ , EP , EPS , EROOT , fi0 , fi1 , fic , FIEX , figl , 
     &       fm , g
      REAL*8 gth , ODL , part , partl , Q , QCEN , rik , rl , rx , ry , 
     &       rys , rz , sf , sgm , SGW , SPIN , SUBCH1 , SUBCH2 , sum3 , 
     &       SUMCL
      REAL*8 sumpr , TACOS , TAU , TETACM , tetrc , tfac , thc , TLBDG , 
     &       TREP , UPL , VACDP , VINF , wf , XA , XA1 , XNOR , YEXP , 
     &       YGN , YGP , YNRM
      REAL*8 ZPOL
      REAL*8 CNOR1 , CNOR2 ! For gosia2
      INTEGER*4 iabc , IAXS , IBYP , Icall , ICLUST , id , idc , Idr , 
     &          IDRN , IEXP , ifdu , IFMO , ifxd , IGRD , ii , ILE , 
     &          ile2 , IMIN , inclus , INM
      INTEGER*4 INNR , ipd , IPRM , IRAWEX , Iredv , ISO , Issp , 
     &          Itemp , ITMA , ITS , iva , iw , IWF , ixl , ixm , IY , 
     &          iyex , IZ , IZ1 , jj
      INTEGER*4 jj1 , jk , jpc , JSKIP , k , k9 , kc , kj , kk , KSEQ , 
     &          KVAR , l , l1 , LASTCL , LFL , LFL1 , LFL2 , lic , 
     &          licz , ll1
      INTEGER*4 LNORM , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , lth , lu , luu , 
     &          na , NANG , NDIM
      INTEGER*4 NDST , NEXPT , nf , nf1 , ni , ni1 , NICC , NLIFT , 
     &          NMAX , NMAX1 , Nwyr , NYLDE
      INTEGER*4 MCFIX ! For gosia2
      CHARACTER*4 wupl , war
      CHARACTER*4 op2 ! For gosia2
      DIMENSION part(32,50,2) , lic(32) , lth(500) , cnr(32,50) , 
     &          partl(32,50,2)
      COMMON /CLUST / ICLUST(50,200) , LASTCL(50,20) , SUMCL(20,500) , 
     &                IRAWEX(50)
      COMMON /ODCH  / DEV(500)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /BREC  / BETAR(50)
      COMMON /DIMX  / DIX(4) , ODL(200)
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      COMMON /CINIT / CNOR(32,75) , CNOR1(32,75) , CNOR2(32,75) , INNR ,
     &                MCFIX ! Changed in gosia2
      COMMON /PRT   / IPRM(20)
      COMMON /LIFE  / NLIFT
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /IGRAD / IGRD
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /MINNI / IMIN , LNORM(50)
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      COMMON /YTEOR / YGN(500) , YGP(500) , IFMO
      COMMON /SEL   / KVAR(500)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /WARN  / SGW , SUBCH1 , SUBCH2 , IWF
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /SKP   / JSKIP(50)
      COMMON /TRB   / ITS
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      COMMON /CCCDS / NDST(50)

      op2 = '    '
      ifxd = 0
      tetrc = TREP(IEXP)
      IF ( Icall.EQ.4 .AND. IPRM(13).EQ.-2 ) THEN
         IPRM(13) = 0
         WRITE (22,99001)
99001    FORMAT (1X//20X,'NORMALIZATION CONSTANTS'//2X,'EXPERIMENT',5X,
     &           'DETECTORS(1-32)')
         DO jpc = 1 , NEXPT
            k = NDST(jpc)
            WRITE (22,99012) jpc , (CNOR(l,jpc),l=1,k)
         ENDDO
         WRITE (22,99002)
99002    FORMAT (1X//20X,'RECOMMENDED RELATIVE GE(LI) EFFICIENCIES'//2X,
     &           'EXPERIMENT')
         DO jpc = 1 , NEXPT
            IF ( ABS(cnr(1,jpc)).LT.1.E-9 ) cnr(1,jpc) = 1.
            k = NDST(jpc)
            WRITE (22,99012) jpc , (cnr(l,jpc)/cnr(1,jpc),l=1,k)
         ENDDO
      ENDIF
      DO jpc = 1 , LP6
         lic(jpc) = 0
      ENDDO
      IF ( Icall.NE.7 ) THEN
         IF ( Itemp.EQ.0 ) THEN
            Nwyr = 0
            IF ( IGRD.NE.1 ) THEN
               IF ( IEXP.EQ.1 ) sumpr = 0.
               IF ( IEXP.EQ.1 ) sum3 = 0.
               DO jj = 1 , LP6
                  DO jk = 1 , 2
                     partl(jj,IEXP,jk) = 0.
                     part(jj,IEXP,jk) = 0.
                  ENDDO
               ENDDO
            ENDIF
            CALL DECAY(Chisq,NLIFT,Chilo)
            IF ( Icall.EQ.4 .AND. IPRM(14).EQ.-1 ) THEN
               IF ( IEXP.EQ.NEXPT ) IPRM(14) = 0
               WRITE (22,99003)
99003          FORMAT (1X//20X,'VACUUM DEPOLARIZATION COEFFICIENTS '//)
               WRITE (22,99004) IEXP
99004          FORMAT (5X,'EXPERIMENT',1X,1I2/5X,'LEVEL',10X,'G2',10X,
     &                 'G4',10X,'G6'/)
               DO iva = 2 , NMAX
                  WRITE (22,99005) iva , (VACDP(ii,iva),ii=1,3)
99005             FORMAT (7X,1I2,9X,3(1F6.4,6X))
               ENDDO
            ENDIF
            fi0 = FIEX(IEXP,1)
            fi1 = FIEX(IEXP,2)
            na = NANG(IEXP)
            DO k = 1 , LP2
               DO k9 = 1 , 20
                  SUMCL(k9,k) = 0.
               ENDDO
            ENDDO
            k9 = 0
            DO k = 1 , na
               gth = AGELI(IEXP,k,1)
               figl = AGELI(IEXP,k,2)
               ifxd = 0
               fm = (fi0+fi1)/2.
               IF ( Icall.EQ.4 ) ifxd = 1
               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k,op2)
               IF ( IFMO.NE.0 ) THEN
                  id = ITMA(IEXP,k)
                  d = ODL(id)
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k,op2)
                  DO ixl = 1 , Idr
                     ixm = KSEQ(ixl,3) ! Initial level of ixl'th decay
                     tfac = TAU(ixm)
                     YGN(ixl) = YGN(ixl) + .01199182*tfac*BETAR(IEXP)
     &                          *(sf*YGP(ixl)-YGN(ixl))
                  ENDDO
               ENDIF
               IF ( IRAWEX(IEXP).NE.0 ) THEN
                  ipd = ITMA(IEXP,k)
                  DO l = 1 , Idr
                     decen = ENDEC(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETAR(IEXP)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLUST(IEXP,k)
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr
                        SUMCL(inclus,l) = SUMCL(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LASTCL(IEXP,inclus) ) GOTO 20
                     DO l = 1 , Idr
                        YGN(l) = SUMCL(inclus,l)
                     ENDDO
                  ENDIF
               ENDIF
               k9 = k9 + 1
               IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                  WRITE (22,99006) IEXP , k9
99006             FORMAT (1X//5X,
     &                 'CALCULATED AND EXPERIMENTAL YIELDS   EXPERIMENT'
     &                 ,1X,1I2,1X,'DETECTOR',1X,1I2//6X,'NI',5X,'NF',7X,
     &                 'II',8X,'IF',9X,'ENERGY(MEV)',6X,'YCAL',8X,
     &                 'YEXP',7X,'PC. DIFF.',2X,'(YE-YC)/SIGMA')
               ENDIF
               lu = ILE(k9)
               DO iabc = 1 , LP2
                  lth(iabc) = 0
               ENDDO
               DO l = 1 , Idr
                  ni = KSEQ(l,3) ! Intial level of l'th decay
                  nf = KSEQ(l,4) ! Final level of l'th decay
                  IF ( l.EQ.IY(lu,k9) .OR. l.EQ.(IY(lu,k9)/1000) ) THEN
                     ifdu = 0
                     lic(k9) = lic(k9) + 1
                     licz = lic(k9)
                     Nwyr = Nwyr + 1
                     wf = CORF(lu,k9)
                     IF ( Icall.EQ.4 ) wf = 1.
                     IF ( Icall.EQ.1 .AND. Issp.EQ.1 ) wf = 1.
                     IF ( IY(lu,k9).GE.1000 ) THEN
                        ifdu = 1
                        l1 = IY(lu,k9)/1000
                        l1 = IY(lu,k9) - 1000*l1
                        YGN(l) = YGN(l) + YGN(l1)
                        lth(l1) = 1
                        IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                           war = '    '
                           sgm = (YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP))
     &                           /DYEX(k9,lu)
                           ni1 = KSEQ(l1,3) ! Initial level of l1'th decay
                           nf1 = KSEQ(l1,4) ! Final level of l1'th decay
                           WRITE (22,99007) ni , ni1 , nf , nf1 , 
     &                            SPIN(ni) , SPIN(ni1) , SPIN(nf) , 
     &                            SPIN(nf1) , ENDEC(l) , ENDEC(l1) , 
     &                            YGN(l)*CNOR(k9,IEXP) , YEXP(k9,lu) , 
     &                            100.*(YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP)
     &                            )/YEXP(k9,lu) , sgm , war
99007                      FORMAT (4X,1I2,'+',1I2,'--',1I2,'+',1I2,3X,
     &                             1F4.1,'+',1F4.1,'--',1F4.1,'+',1F4.1,
     &                             3X,1F6.4,'+',1F6.4,2X,1E9.4,6X,1E9.4,
     &                             3X,1F6.1,5X,1F4.1,10X,1A4)
                           SUBCH1 = SUBCH1 + sgm*sgm
                        ENDIF
                     ENDIF
                     ry = YGN(l)*wf*CNOR(k9,IEXP) - YEXP(k9,lu)
                     IF ( ifdu.NE.1 ) THEN
                        IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                           war = '    '
                           sgm = (YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP))
     &                           /DYEX(k9,lu)
                           WRITE (22,99013) ni , nf , SPIN(ni) , 
     &                            SPIN(nf) , ENDEC(l) , YGN(l)
     &                            *CNOR(k9,IEXP) , YEXP(k9,lu) , 
     &                            100.*(YEXP(k9,lu)-YGN(l)*CNOR(k9,IEXP)
     &                            )/YEXP(k9,lu) , sgm , war
                           SUBCH1 = SUBCH1 + sgm*sgm
                        ENDIF
                     ENDIF
                     rys = ry*ry
                     IF ( IGRD.EQ.1 ) Chisq = Chisq + rys/DYEX(k9,lu)
     &                    /DYEX(k9,lu)
                     IF ( k9.EQ.1 .AND. Iredv.EQ.1 ) DEV(licz) = ry
                     IF ( Iredv.NE.1 ) THEN
                        IF ( LFL.EQ.1 ) THEN
                           IF ( k9.EQ.1 ) THEN
                              luu = 6*licz - 5
                              jk = (luu-1)/LP10 + 1
                              kk = luu - LP10*(jk-1)
                              rik = DEV(licz) + YEXP(k9,lu)
                              sgm = -DEV(licz)/DYEX(k9,lu)
                              IF ( ITS.EQ.1 .AND. KVAR(INM).NE.0 )
     &                             WRITE (17,*) ni , nf , sgm , YGN(l)
     &                             *CNOR(k9,IEXP)/DYEX(k9,lu)
                              IF ( ITS.EQ.1 .AND. INM.EQ.1 )
     &                             WRITE (15,*) IEXP , rik/CNOR(1,IEXP)
     &                             , CNOR(1,IEXP) , DYEX(k9,lu) , 
     &                             YEXP(k9,lu)
                              CALL SIXEL(rik,ry,EMH,jk,kk,INM,licz)
                           ENDIF
                        ENDIF
                     ENDIF
                     IF ( IGRD.NE.1 ) THEN
                        IF ( JSKIP(IEXP).NE.0 ) THEN
                           dl = DYEX(k9,lu)*DYEX(k9,lu)
                           part(k9,IEXP,1) = part(k9,IEXP,1) + YGN(l)
     &                        *YGN(l)*wf*wf/dl
                           part(k9,IEXP,2) = part(k9,IEXP,2) - 2.*YGN(l)
     &                        *wf*YEXP(k9,lu)/dl
                           sumpr = sumpr + YEXP(k9,lu)*YEXP(k9,lu)/dl
                           partl(k9,IEXP,1) = partl(k9,IEXP,1)
     &                        + YEXP(k9,lu)*YEXP(k9,lu)/dl
                           partl(k9,IEXP,2) = partl(k9,IEXP,2)
     &                        + LOG(wf*YGN(l)/YEXP(k9,lu))*YEXP(k9,lu)
     &                        *YEXP(k9,lu)/dl
                           sum3 = sum3 + YEXP(k9,lu)*YEXP(k9,lu)
     &                            *LOG(wf*YGN(l)/YEXP(k9,lu))**2/dl
                        ENDIF
                     ENDIF
                     lu = lu + 1
                  ELSE
                     IF ( JSKIP(IEXP).EQ.0 ) YGN(IDRN) = 1.E+10
                     ry = YGN(l)/YGN(IDRN)
                     IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                        wupl = '    '
                        IF ( ry.GT.UPL(k9,IEXP) .AND. lth(l).EQ.0 )
     &                       wupl = 'UPL!'
                        IF ( IPRM(16).NE.0 .OR. wupl.NE.'    ' ) THEN
                           IF ( wupl.EQ.'    ' ) WRITE (22,99008) ni , 
     &                          nf , SPIN(ni) , SPIN(nf) , ENDEC(l) , 
     &                          YGN(l)*CNOR(k9,IEXP) , wupl
99008                      FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,
     &                             1F6.4,6X,1E9.4,10X,1A4)
                           IF ( wupl.NE.'    ' ) THEN
                              sgm = (ry-UPL(k9,IEXP))/UPL(k9,IEXP)
                              WRITE (22,99013) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , ENDEC(l) , YGN(l)
     &                               *CNOR(k9,IEXP) , UPL(k9,IEXP)
     &                               *CNOR(k9,IEXP)*YGN(IDRN) , 
     &                               100.*(1.-YGN(l)/UPL(k9,IEXP)
     &                               /YGN(IDRN)) , sgm , wupl
                              SUBCH1 = SUBCH1 + sgm*sgm
                           ENDIF
                        ENDIF
                     ENDIF
                     IF ( ry.GE.UPL(k9,IEXP) .AND. lth(l).NE.1 ) THEN
                        Chisq = Chisq + (ry-UPL(k9,IEXP))
     &                          *(ry-UPL(k9,IEXP))/UPL(k9,IEXP)
     &                          /UPL(k9,IEXP)
                        Chilo = Chilo + LOG(ry/UPL(k9,IEXP))**2
                        IF ( IWF.NE.0 ) THEN
                           WRITE (22,99009) IEXP , ni , nf , 
     &                            ry/UPL(k9,IEXP)
99009                      FORMAT (5X,'WARNINIG-EXP.',1I2,2X,'TRANS. ',
     &                             1I2,'--',1I2,5X,
     &                             'EXCEEDS UPPER LIMIT (RATIO=',1E14.6,
     &                             ')')
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO
               IF ( IEXP.EQ.NEXPT ) IWF = 0
               IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                  WRITE (22,99010) SUBCH1 - SUBCH2
99010             FORMAT (1X/50X,'CHISQ SUBTOTAL = ',E14.6)
                  SUBCH2 = SUBCH1
               ENDIF
 20         ENDDO
            IF ( IGRD.EQ.1 ) RETURN
            IF ( IEXP.NE.NEXPT ) RETURN
            IF ( Icall.EQ.1 ) RETURN
         ELSE
            ifxd = 1
            IF ( Itemp.NE.2 ) ifxd = 0
            Nwyr = 1
            CALL DECAY(ccd,0,ccc)
            fi0 = FIEX(IEXP,1)
            fi1 = FIEX(IEXP,2)
            na = NANG(IEXP)
            DO k = 1 , LP2
               DO kj = 1 , 20
                  SUMCL(kj,k) = 0
               ENDDO
            ENDDO
            k9 = 0
            DO k = 1 , na
               gth = AGELI(IEXP,k,1)
               figl = AGELI(IEXP,k,2)
               fm = (fi0+fi1)/2.
               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k, op2)
               IF ( IFMO.NE.0 ) THEN
                  id = ITMA(IEXP,k)
                  d = ODL(id)
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k,op2)
                  DO ixl = 1 , Idr
                     ixm = KSEQ(ixl,3) ! Initial level of ixl'th decay
                     tfac = TAU(ixm)
                     IF ( tfac.GT.1.E+4 ) GOTO 25
                     YGN(ixl) = YGN(ixl) + .01199182*tfac*BETAR(IEXP)
     &                          *(sf*YGP(ixl)-YGN(ixl))
                  ENDDO
 25               IFMO = 0
                  WRITE (22,99011)
99011             FORMAT (1X,/,2X,'DURING THE MINIMIZATION',1X,
     &    'IT WAS NECESSARY TO SWITCH OFF THE TIME-OF-FLIGHT CORRECTION'
     &    )
               ENDIF
               IF ( IRAWEX(IEXP).NE.0 ) THEN
                  ipd = ITMA(IEXP,k)
                  DO l = 1 , Idr
                     decen = ENDEC(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETAR(IEXP)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLUST(IEXP,k)
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr
                        SUMCL(inclus,l) = SUMCL(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LASTCL(IEXP,inclus) ) GOTO 40
                     DO l = 1 , Idr
                        YGN(l) = SUMCL(inclus,l)
                     ENDDO
                  ENDIF
               ENDIF
               k9 = k9 + 1
               iyex = NYLDE(IEXP,k9) + ILE(k9) - 1
               ile2 = ILE(k9)
               DO l = ile2 , iyex
                  IF ( JSKIP(IEXP).NE.0 ) THEN
                     idc = IY(l,k9)
                     IF ( idc.GE.1000 ) THEN
                        idc = idc/1000
                        ll1 = IY(l,k9) - idc*1000
                        YGN(idc) = YGN(idc) + YGN(ll1)
                     ENDIF
                     IF ( Itemp.EQ.1 ) THEN
                        CORF(l,k9) = CORF(l,k9)/(YGN(idc)+1.E-24)
                     ELSE
                        CORF(l,k9) = YGN(idc)
                        IF ( IMIN.LE.1 .AND. l.EQ.iyex ) CNOR(k9,IEXP)
     &                       = YEXP(k9,l)/YGN(idc)
                     ENDIF
                  ENDIF
               ENDDO
 40         ENDDO
            RETURN
         ENDIF
      ENDIF
      DO jj = 1 , NEXPT
         IF ( JSKIP(jj).NE.0 ) THEN
            kc = NDST(jj)
            DO jk = 1 , kc
               cnr(jk,jj) = -.5*part(jk,jj,2)/part(jk,jj,1)
               IF ( INNR.NE.0 ) CNOR(jk,jj) = cnr(jk,jj)
            ENDDO
            IF ( INNR.NE.1 ) THEN
               d = 0.
               g = 0.
               DO jj1 = jj , NEXPT
                  IF ( LNORM(jj1).EQ.jj ) THEN
                     k = NDST(jj1)
                     DO jk = 1 , k
                        d = d + YNRM(jk,jj1)*part(jk,jj1,1)*YNRM(jk,jj1)
                        g = g - .5*YNRM(jk,jj1)*part(jk,jj1,2)
                     ENDDO
                  ENDIF
               ENDDO
               IF ( LNORM(jj).EQ.jj ) THEN
                  CNOR(1,jj) = g*YNRM(1,jj)/d
                  k = NDST(jj)
                  IF ( k.NE.1 ) THEN
                     DO jk = 2 , k
                        CNOR(jk,jj) = YNRM(jk,jj)*CNOR(1,jj)/YNRM(1,jj)
                     ENDDO
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF ( INNR.NE.1 ) THEN
         DO jj = 1 , NEXPT
            IF ( LNORM(jj).NE.jj ) THEN
               iw = LNORM(jj)
               k = NDST(jj)
               DO jk = 1 , k
                  CNOR(jk,jj) = CNOR(1,iw)*YNRM(jk,jj)/YNRM(1,iw)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
      IF ( Icall.EQ.7 ) Chisq = 0.
      DO jj = 1 , NEXPT
         k = NDST(jj)
         DO jk = 1 , k
            IF ( MCFIX .EQ. 0) CNOR(jk,jj) = CNOR1(jk,jj)
            Chilo = Chilo + partl(jk,jj,1)*LOG(CNOR(jk,jj))
     &              **2 + partl(jk,jj,2)*2.*LOG(CNOR(jk,jj))
            Chisq = Chisq + CNOR(jk,jj)*CNOR(jk,jj)*part(jk,jj,1)
     &              + CNOR(jk,jj)*part(jk,jj,2)
         ENDDO
      ENDDO
      Chisq = Chisq + sumpr
      Chilo = Chilo + sum3
      RETURN
99012 FORMAT (1X,1I2,2X,32(1E10.4,1X))
99013 FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,1F6.4,6X,1E9.4,6X,
     &        1E9.4,3X,1F6.1,5X,1F4.1,10X,1A4)
      END
