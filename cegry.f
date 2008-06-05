 
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
C      ENDEC  - energy difference for each matrix element
C      FIEX   - phi range of particle detector
C      ICLUST - cluster number for each experiment and detector
C      IDRN   - index of normalising transition for yields
C      IEXP   - number of experiment
C      IFMO   - include correction to angular distance for finite recoil distance.
C      IGRD   -
C      ILE    - yield number for each detector
C      IMIN   -
C      INM    - index of matrix element
C      INNR   - independent normalisation switch (see OP,CONT INR,)
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      IRAWEX -
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      IWF    - warning flag
C      IY     - index for yields
C      JSKIP  - Experiments to skip during minimisation.
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      KVAR   -
C      LASTCL - index of last detector in cluster
C      LFL    -
C      LNORM  - normalization constant control
C      LP2    - maximum number of matrix elements (1500)
C      LP6    - 32
C      LP10   - 600
C      MCFIX  - fixing parameter (added for gosia2)
C      NANG   - number of gamma-ray detectors for each experiment
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NLIFT  - number of lifetimes
C      NMAX   - number of levels
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      SGW    - number of standard deviations to generate warning (see control option WRN,X)
C      SPIN   - spin of level
C      SUBCH1 - partial chisqr
C      SUBCH2 - partial chisqr
C      SUMCL  - sum of yields for clusters
C      TAU    - lifetime in picoseconds
C      TREP   - theta of recoiling nucleus
C      UPL    - upper limits for all gamma detectors
C      VACDP  - G_k for each level
C      YEXP   - experimental yield
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      YGP    - gamma yield calculated with correction to angular distribution from finite recoil distance
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
      REAL*8 ccc , ccd , Chilo , Chisq , cnr , cocos , d , decen
      REAL*8 dl , effi , EMH , EP , fi0 , fi1 , fic , figl , 
     &       fm , g
      REAL*8 gth , part , partl , rik , rl , rx , ry , 
     &       rys , rz , sf , sgm , SGW , SUBCH1 , SUBCH2 , sum3
      REAL*8 sumpr , TACOS , tetrc , tfac , thc , TLBDG , 
     &       VINF , wf , XA , XA1
      INTEGER*4 iabc , Icall , id , idc , Idr ,  ifdu , ifxd , IGRD , 
     &          ii , ile2 , IMIN , inclus , INM
      INTEGER*4 ipd , IPRM , Iredv , Issp , 
     &          Itemp , ITS , iva , iw , IWF , ixl , ixm , 
     &          iyex , IZ , IZ1 , jj
      INTEGER*4 jj1 , jk , jpc , JSKIP , k , k9 , kc , kj , kk , 
     &          KVAR , l , l1 , LFL , LFL1 , LFL2 , lic , 
     &          licz , ll1
      INTEGER*4 LNORM , lth , lu , luu , na , NDIM , NEXPT
      INTEGER*4 nf , nf1 , ni , ni1 , NLIFT , NMAX , NMAX1 , Nwyr
      CHARACTER*4 wupl , war
      REAL*8 CNOR1 , CNOR2 ! For gosia2
      INTEGER*4 MCFIX ! For gosia2
      CHARACTER*4 op2 ! For gosia2
      DIMENSION part(32,50,2) , lic(32) , lth(1500) , cnr(32,50) , 
     &          partl(32,50,2)
      INCLUDE 'clust.inc'
      INCLUDE 'odch.inc'
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      INCLUDE 'tra.inc'
      INCLUDE 'brec.inc'
      INCLUDE 'dimx.inc'
      INCLUDE 'vac.inc'
      INCLUDE 'cinit.inc'
      COMMON /PRT   / IPRM(20)
      COMMON /LIFE  / NLIFT
      INCLUDE 'lev.inc'
      COMMON /IGRAD / IGRD
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /MINNI / IMIN , LNORM(50)
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      INCLUDE 'yteor.inc'
      COMMON /SEL   / KVAR(1500)
      INCLUDE 'mgn.inc'
      INCLUDE 'ccc.inc'
      INCLUDE 'yexpt.inc'
      INCLUDE 'kin.inc'
      COMMON /WARN  / SGW , SUBCH1 , SUBCH2 , IWF
      INCLUDE 'coex.inc'
      COMMON /SKP   / JSKIP(50)
      COMMON /TRB   / ITS
      INCLUDE 'tcm.inc'
      INCLUDE 'cccds.inc'
      COMMON /CINIT2/ CNOR1(32,75) , CNOR2(32,75) , MCFIX ! For gosia2
      SAVE part, partl ! For gosia2
      DATA sum3/0./,sumpr/0./

      op2 = '    ' ! For gosia2
      ifxd = 0
      tetrc = TREP(IEXP) ! Theta of recoiling nucleus

C     If the user set print flag 13 to +1, it is set to -1 by OP,EXIT and then
C     if it is -1, it is set to -2 in MINI, which is called from there, which
C     in turn calls FTBM, which calls this function. In other words, this
C     routine is called with IPRM(13) set to -2 if the user sets IPRM(13) to 1
C     with CONT:PRT, and then does OP,EXIT

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
         ENDDO ! Loop on experiments
      ENDIF ! if Icall.EQ.4 .AND. IPRM(13).EQ.-2

      DO jpc = 1 , LP6 ! LP6 is 32
         lic(jpc) = 0
      ENDDO

      IF ( Icall.NE.7 ) THEN
         IF ( Itemp.EQ.0 ) THEN
            Nwyr = 0
            IF ( IGRD.NE.1 ) THEN
               IF ( IEXP.EQ.1 ) sumpr = 0.
               IF ( IEXP.EQ.1 ) sum3 = 0.
               DO jj = 1 , LP6 ! LP6 is 32
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

            fi0 = FIEX(IEXP,1) ! Lower phi limit
            fi1 = FIEX(IEXP,2) ! Upper phi limit
            na = NANG(IEXP) ! Number of detector angles

            DO k = 1 , LP2 ! LP2 is 1500
               DO k9 = 1 , 20
                  SUMCL(k9,k) = 0.
               ENDDO
            ENDDO

            k9 = 0
            DO k = 1 , na ! For each detector angle
               gth = AGELI(IEXP,k,1) ! theta
               figl = AGELI(IEXP,k,2) ! phi
               ifxd = 0
               fm = (fi0+fi1)/2.
               IF ( Icall.EQ.4 ) ifxd = 1
               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k,op2) ! Changed in gosia2

C              Correct for finite recoil
               IF ( IFMO.NE.0 ) THEN
                  id = ITMA(IEXP,k) ! Get identity for detector
                  d = ODL(id) ! Results of OP,GDET for this detector
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k,op2) ! Changed in gosia2
                  DO ixl = 1 , Idr ! For each decay
                     ixm = KSEQ(ixl,3) ! Initial level of ixl'th decay
                     tfac = TAU(ixm) ! Get lifetime
                     YGN(ixl) = YGN(ixl) + .01199182*tfac*BETAR(IEXP)
     &                          *(sf*YGP(ixl)-YGN(ixl))
                  ENDDO ! Loop on decays ixl
               ENDIF ! If correction for finite recoil

               IF ( IRAWEX(IEXP).NE.0 ) THEN
                  ipd = ITMA(IEXP,k) ! Get identity for detector
                  DO l = 1 , Idr
                     decen = ENDEC(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETAR(IEXP)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLUST(IEXP,k) ! Cluster number for detector k
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr ! For each decay
                        SUMCL(inclus,l) = SUMCL(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LASTCL(IEXP,inclus) ) GOTO 20 ! If it is not the last detector in the cluster
                     DO l = 1 , Idr ! For each decay
                        YGN(l) = SUMCL(inclus,l)
                     ENDDO
                  ENDIF
               ENDIF
               k9 = k9 + 1 ! Increment detector number
               IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                  WRITE (22,99006) IEXP , k9
99006             FORMAT (1X//5X,
     &                 'CALCULATED AND EXPERIMENTAL YIELDS   EXPERIMENT'
     &                 ,1X,1I2,1X,'DETECTOR',1X,1I2//6X,'NI',5X,'NF',7X,
     &                 'II',8X,'IF',9X,'ENERGY(MEV)',6X,'YCAL',8X,
     &                 'YEXP',7X,'PC. DIFF.',2X,'(YE-YC)/SIGMA')
               ENDIF
               lu = ILE(k9) ! Yield number for detector k9
               DO iabc = 1 , LP2 ! LP2 = 1500
                  lth(iabc) = 0
               ENDDO
               DO l = 1 , Idr ! For each decay
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
                           IF ( ABS(sgm).GE.SGW ) war = '*?!*'
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
                           IF ( ABS(sgm).GE.SGW ) war = '*?!*'
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
                              jk = (luu-1)/LP10 + 1 ! LP10 is 600
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
                        IF ( IWF.NE.0 ) THEN ! If warning flag is set
                           WRITE (22,99009) IEXP , ni , nf , 
     &                            ry/UPL(k9,IEXP)
99009                      FORMAT (5X,'WARNING-EXP.',1I2,2X,'TRANS. ',
     &                             1I2,'--',1I2,5X,
     &                             'EXCEEDS UPPER LIMIT (RATIO=',1E14.6,
     &                             ')')
                        ENDIF
                     ENDIF
                  ENDIF
               ENDDO ! Loop on decays l
               IF ( IEXP.EQ.NEXPT ) IWF = 0 ! Turn off warnings now
               IF ( Icall.EQ.4 .AND. IPRM(8).EQ.-2 ) THEN
                  WRITE (22,99010) SUBCH1 - SUBCH2
99010             FORMAT (1X/50X,'CHISQ SUBTOTAL = ',E14.6)
                  SUBCH2 = SUBCH1
               ENDIF
 20         ENDDO ! Loop on detector angles k

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
            DO k = 1 , LP2 ! LP2 is 1500
               DO kj = 1 , 20
                  SUMCL(kj,k) = 0
               ENDDO
            ENDDO
            k9 = 0
            DO k = 1 , na
               gth = AGELI(IEXP,k,1)
               figl = AGELI(IEXP,k,2)
               fm = (fi0+fi1)/2.

               CALL ANGULA(YGN,Idr,ifxd,fi0,fi1,tetrc,gth,figl,k,op2) ! Changed in gosia2

C              Correct for finite recoil
               IF ( IFMO.NE.0 ) THEN
                  id = ITMA(IEXP,k) ! Get identity for detector
                  d = ODL(id) ! Get results of OP,GDET for detector
                  rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(tetrc)*COS(fm)
                  ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(tetrc)*SIN(fm)
                  rz = d*COS(gth) - .25*COS(tetrc)
                  rl = SQRT(rx*rx+ry*ry+rz*rz)
                  sf = d*d/rl/rl
                  thc = TACOS(rz/rl)
                  fic = ATAN2(ry,rx)
                  CALL ANGULA(YGP,Idr,ifxd,fi0,fi1,tetrc,thc,fic,k,op2) ! Changed in gosia2
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
               ENDIF ! if correction for finite recoil

               IF ( IRAWEX(IEXP).NE.0 ) THEN
                  ipd = ITMA(IEXP,k) ! Get identity of detector
                  DO l = 1 , Idr
                     decen = ENDEC(l)
                     cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl)
     &                       + COS(tetrc)*COS(gth)
                     decen = decen*(1.+BETAR(IEXP)*cocos)
                     CALL EFFIX(ipd,decen,effi)
                     YGN(l) = YGN(l)*effi
                  ENDDO
                  inclus = ICLUST(IEXP,k) ! Cluster number for detector k
                  IF ( inclus.NE.0 ) THEN
                     DO l = 1 , Idr ! For each decay
                        SUMCL(inclus,l) = SUMCL(inclus,l) + YGN(l)
                     ENDDO
                     IF ( k.NE.LASTCL(IEXP,inclus) ) GOTO 40 ! If it is not the last detector in the cluster
                     DO l = 1 , Idr ! For each decay
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
               ENDDO ! Loop on l
 40         ENDDO ! Loop on k
            RETURN
         ENDIF ! if Itemp.EQ.0
      ENDIF ! if Icall.NE.7

C     Sort out normalisation coefficients
      DO jj = 1 , NEXPT ! For each experiment
         IF ( JSKIP(jj).NE.0 ) THEN
            kc = NDST(jj) ! Number of datasets
            DO jk = 1 , kc ! For each dataset
               cnr(jk,jj) = -.5*part(jk,jj,2)/part(jk,jj,1)
               IF ( INNR.NE.0 ) CNOR(jk,jj) = cnr(jk,jj)
            ENDDO ! Loop on datasets

C           If we want a common normalisation, sort it out here
            IF ( INNR.NE.1 ) THEN
               d = 0.
               g = 0.
               DO jj1 = jj , NEXPT ! For each experiment
                  IF ( LNORM(jj1).EQ.jj ) THEN
                     k = NDST(jj1) ! Number of datasets
                     DO jk = 1 , k ! For each dataset
                        d = d + YNRM(jk,jj1)*part(jk,jj1,1)*YNRM(jk,jj1)
                        g = g - .5*YNRM(jk,jj1)*part(jk,jj1,2)
                     ENDDO ! Loop on datasets
                  ENDIF ! IF ( LNORM(jj1).EQ.jj )
               ENDDO ! Loop on experiment
               IF ( LNORM(jj).EQ.jj ) THEN ! If this is the normalisation transition
                  CNOR(1,jj) = g*YNRM(1,jj)/d
                  k = NDST(jj) ! Number of datasets
                  IF ( k.NE.1 ) THEN
                     DO jk = 2 , k ! For each dataset
                        CNOR(jk,jj) = YNRM(jk,jj)*CNOR(1,jj)/YNRM(1,jj)
                     ENDDO ! Loop on jk
                  ENDIF ! IF ( k.NE.1 )
               ENDIF ! IF ( LNORM(jj).EQ.jj )
            ENDIF ! IF ( INNR.NE.1 )
             
         ENDIF ! IF ( JSKIP(jj).NE.0 )
      ENDDO ! Loop on experiment

C     If there is a common normalisation, normalise to it
      IF ( INNR.NE.1 ) THEN
         DO jj = 1 , NEXPT ! For each experiment
            IF ( LNORM(jj).NE.jj ) THEN
               iw = LNORM(jj) ! Get index of normalisation transition
               k = NDST(jj) ! Get number of datasets
               DO jk = 1 , k ! For each dataset
                  CNOR(jk,jj) = CNOR(1,iw)*YNRM(jk,jj)/YNRM(1,iw)
               ENDDO ! Loop on datasets
            ENDIF ! IF ( LNORM(jj).NE.jj )
         ENDDO ! Loop on experiment
      ENDIF ! IF ( INNR.NE.1 )

C     Calculate chi squared
      IF ( Icall.EQ.7 ) Chisq = 0.
      DO jj = 1 , NEXPT
         k = NDST(jj)
         DO jk = 1 , k
            IF ( MCFIX .EQ. 0) CNOR(jk,jj) = CNOR1(jk,jj) ! For gosia2
            Chilo = Chilo + partl(jk,jj,1)*LOG(CNOR(jk,jj))
     &              **2 + partl(jk,jj,2)*2.*LOG(CNOR(jk,jj))
            Chisq = Chisq + CNOR(jk,jj)*CNOR(jk,jj)*part(jk,jj,1)
     &              + CNOR(jk,jj)*part(jk,jj,2)
         ENDDO ! Loop on datasets
      ENDDO ! Loop on experiment

      Chisq = Chisq + sumpr
      Chilo = Chilo + sum3
      RETURN

99012 FORMAT (1X,1I2,2X,32(1E10.4,1X))
99013 FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,1F6.4,6X,1E9.4,6X,
     &        1E9.4,3X,1F6.1,5X,1F4.1,10X,1A4)
      END
