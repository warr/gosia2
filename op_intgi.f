
C----------------------------------------------------------------------
C SUBROUTINE OP_INTGI
C
C Called by: GOSIA
C Calls: ALLOC, ANGULA, CMLAB, COORD, DECAY, EFFIX, INTG, INVKIN, LAGRAN,
C        LOAD, PATH, SETIN, SIMIN, SNAKE, SPLNER, STING, TACOS, TAPMA,
C        TENB, TENS
C
C Purpose: perform the integration over the meshpoints for both OP,INTG
C          and OP,INTI
C
C Uses global variables:
C      ACCUR  - accuracy required
C      AGELI  - angles of the Ge detectors
C      ELMH   - matrix element
C      BETAR  - recoil beta
C      DS     - integrated rutherford cross-section
C      DSE    - rutherford cross section at given energy integrated over angles
C      DSG    - differential gamma-ray yield at meshpoints
C      EN     - energy of level
C      EP     - bombarding energy
C      ERR    - error flag
C      FIEX   - phi range of particle detector
C      GRAD   - partial derivative of chi squared wrt. each matrix element
C      HLM    - matrix elements before minimisation
C      HLMLM  - old value of matrix element or chi squared
C      IAXS   - axial symmetry flag
C      ICLUST - cluster number for each experiment and detector
C      IDRN   - index of normalising transition for yields
C      IEXP   - experiment number
C      IFMO   - include correction to angular distance for finite recoil distance.
C      ILE    - yield number for each detector
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      IRAWEX - flag to indicate raw uncorrected yield
C      ISKIN  - kinematic flag (0,1)
C      ISO    - isotropic flag
C      ISPL   - spline flag
C      ITMA   - identify detectors according to OP,GDET
C      IZ1    - Z of not-investated nucleus
C      JZB    - unit to read from
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LASTCL - index of last detector in cluster
C      LMAX   - ground-state spin + 1
C      LP2    - maximum number of matrix elements (1500)
C      LP6    - maximum number of Ge detectors 32
C      NANG   - number of gamma-ray detectors for each experiment
C      NCM    - calculate kinematics assuming this state for final state (default = 2)
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      SPIN   - spin of level
C      SUMCL  - sum of yields for clusters
C      TAU    - lifetime in picoseconds
C      TLBDG  - theta of particle detector in lab frame (in degrees)
C      XA     - A of investigated nucleus
C      XA1    - A of not-investated nucleus
C      XI     - xi coupling coefficients
C      XV     - energy meshpoints (sometimes theta meshpoints) where we calculate exact Coulex
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      YGP    - gamma yield calculated with correction to angular distribution from finite recoil distance
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C      ZETA   - various coefficients
C      ZPOL   - dipole term (GDR excitation)
C
C The Inti flag is 0 for OP,INTG and 1 for OP,INTI. The only difference
C between the two is that OP,INTI calls INVKIN and OP,INTG doesn't. In both
C we perform the calculation for each meshpoint and then integrate by
C interpolating over the meshpoints

      SUBROUTINE OP_INTGI(Iretval, Op2, Ipinf, Jpin, Iecd, Izcap,
     &  Lfagg, Idr, Bten, Fiex1, Inti)
      IMPLICIT NONE
      INCLUDE 'brec.inc'
      INCLUDE 'caux0.inc'
      INCLUDE 'ccc.inc'
      INCLUDE 'cccds.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'clcom9.inc'
      INCLUDE 'clm.inc'
      INCLUDE 'clust.inc'
      INCLUDE 'coex.inc'
      INCLUDE 'cx.inc'
      INCLUDE 'cxi.inc'
      INCLUDE 'dimx.inc'
      INCLUDE 'dumm.inc'
      INCLUDE 'hhh.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'lev.inc'
      INCLUDE 'mgn.inc'
      INCLUDE 'prt.inc'
      INCLUDE 'seck.inc'
      INCLUDE 'switch.inc'
      INCLUDE 'tra.inc'
      INCLUDE 'vlin.inc'
      INCLUDE 'yexpt.inc'
      INCLUDE 'yteor.inc'
      INCLUDE 'fconst.inc'

      CHARACTER*4 Op2
      REAL*8 ax , Bten(1600), ccc , ccd , cocos , d , decen , dedx(20) ,
     &  dsig , dst , dsx , dsxm(100,100,100) , effi , emn , emx , enb ,
     &  enh , esp(20) , fi0 , fi1 , fic , Fiex1(100,100,2), figl, fm ,
     &  gth , hen , het , pfi(101), polm , rl , rx , ry , rz , SIMIN ,
     &  sf , TACOS , tetrc , tfac , thc ,  tmn , tmx , todfi , tta ,
     &  tth , tting , txx , wph , wpi(100,2) , wth , wthh , xx , yy , zz
      INTEGER*4 i , icll , id , Idr , Iecd(50) , ija0 , ijaja , ijan ,
     &  ilx , iocc , inclus , Inti, ipd , Ipinf , Iretval , iske ,
     &  iskf , iskin_protect(50) , isko , ixl , ixm , Izcap , j, ja ,
     &  jan1 , jan , jd , jdy , je , jfi , jgd , jj, jjlx , jkloo ,
     &  jktt , jmpin , Jpin(50) , jt , jtp , jyi , jyv , kloop , ktt ,
     &  Lfagg , locat , lpin , lu , lx , mfla , mpin , na , na1 , naa ,
     &  ndum , ne , nf , nfi , nflr , nft , nged , ngpr , ni , npce ,
     &  npce1 , npct , npct1 , npt , nptx , ntt

      Iretval = 100
      DO lx = 1 , NEXPT ! For each experiment store original ISKIN
        iskin_protect(lx) = ISKIN(lx)
      ENDDO
      REWIND 14
      Lfagg = 1
      IF ( SPIN(1).LT..25 ) ISO = 0
      DO lx = 1 , NEXPT ! For each experiment
        lpin = 1
        IF ( Ipinf.NE.0 ) THEN
          IF ( Jpin(lx).NE.0 ) lpin = Jpin(lx)
        ENDIF
        IEXP = lx
        tth = TLBDG(lx)
        enh = EP(lx)
        DO mpin = 1 , lpin ! For each pin diode
          IF ( Iecd(lx).EQ.1 ) THEN ! Circular detector
            READ (JZB,*) ne , ntt , emn , emx , wth , wph , wthh
            mfla = 1
            CALL COORD(wth,wph,wthh,ntt,0,pfi,wpi,tth,lx,tmn,tmx)
          ELSE
            READ (JZB,*) ne , ntt , emn , emx , tmn , tmx
            mfla = 0
            IF ( ntt.LT.0 ) mfla = 1
          ENDIF
          ntt = ABS(ntt)
          jan = NANG(lx)
          jan1 = NDST(lx)
          IF ( IRAWEX(lx).EQ.0 ) jan1 = jan
          IF ( Iecd(lx).EQ.1 ) THEN ! Circular detector
            WRITE (14,*) ne , ntt , emn , emx , tmn , tmx , jan1 ,
     &        wth , wph , wthh
          ELSE
            WRITE (14,*) ne , ntt , emn , emx , tmn , tmx , jan1 ,
     &        tmx , tmx , tmx
          ENDIF
          READ (JZB,*) (XV(i),i=1,ne)
          IF ( Iecd(lx).NE.1 ) READ (JZB,*) (YV(i),i=1,ntt)
          IF ( tth.LT.0. ) ELMH(2*lx-1) = YV(1)
          IF ( tth.LT.0. ) ELMH(2*lx) = YV(ntt)
          DO kloop = 1 , ne ! For each energy meshpoint
            enb = XV(kloop)
            EP(lx) = enb
            DO ktt = 1 , ntt
              tta = YV(ktt)
              IF ( tth.LT.0 .AND. Inti.EQ.1 )
     &          CALL INVKIN(EP(lx),EN(NCM),IZ1(lx),XA,XA1(lx),YV(ktt),
     &          tta,1,ISKIN(lx))
              tta = SIGN(tta, tth)
              IF ( IAXS(lx).NE.0 ) THEN ! If not axial symmetry
                IF ( Iecd(lx).NE.1 ) THEN
                  IF ( kloop.EQ.1 ) THEN
                    READ (JZB,*) nfi ! Number of phi ranges
                    READ (JZB,*) (Fiex1(ktt,jfi,1),Fiex1(ktt,jfi,2),
     &                jfi=1,nfi)
                    IF ( tth.LT.0. ) THEN
                      DO jfi = 1 , nfi ! For each phi angle
                        Fiex1(ktt,jfi,1) = Fiex1(ktt,jfi,1) + 180.D0
                        Fiex1(ktt,jfi,2) = Fiex1(ktt,jfi,2) + 180.D0
                      ENDDO
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF ! If not axial symmetry
              TLBDG(lx) = tta
              IF ( kloop.EQ.1 ) THEN
                IF ( Iecd(lx).NE.0 ) THEN
                  nfi = 1
                  Fiex1(ktt,1,1) = wpi(ktt,1) ! Lower phi limit
                  Fiex1(ktt,1,2) = wpi(ktt,2) ! Upper phi limit
                ENDIF
              ENDIF
              CALL CMLAB(lx,dsig,tetrc)
              IF ( ERR ) THEN
                Iretval=2000 ! Error
                RETURN
              ENDIF
              tting = TLBDG(lx)
              IF ( ERR ) THEN
                Iretval=1900 ! Troubleshoot
                RETURN
              ENDIF
              CALL LOAD(lx,1,1,0.D0,jj)
              CALL ALLOC(ACCUR)
              CALL SNAKE(lx,ZPOL)
              CALL SETIN
              DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                polm = DBLE(j-1) - SPIN(1)
                CALL LOAD(lx,2,1,polm,jj)
                CALL STING(jj)
                CALL PATH(jj)
                CALL INTG(IEXP)
                CALL TENB(j,Bten,LMAX)
              ENDDO
              CALL TENS(Bten)
              CALL DECAY(ccd,0,ccc)
              DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                DO ijan = 1 , 20
                  SUMCL(ijan,j) = 0.
                ENDDO
              ENDDO
              ija0 = 0
              DO ijan = 1 , jan ! For each detector angle
                IF ( IAXS(lx).EQ.0 ) nfi = 1
                DO jyi = 1 , Idr
                  GRAD(jyi) = 0.
                ENDDO
                todfi = 0.
                DO jfi = 1 , nfi ! For each phi angle
                  fi0 = Fiex1(ktt,jfi,1)*pi/180.D0
                  fi1 = Fiex1(ktt,jfi,2)*pi/180.D0
                  gth = AGELI(IEXP,ijan,1)
                  fm = (fi0+fi1)/2.
                  figl = AGELI(IEXP,ijan,2)
                  CALL ANGULA(YGN,Idr,1,fi0,fi1,tetrc,gth,figl,ijan,Op2)
                  IF ( IFMO.NE.0 ) THEN ! If correction due to recoil
                    id = ITMA(IEXP,ijan) ! Get detector identity
                    d = ODL(id) ! Get result of OP,GDET calculation
                    rx = d*SIN(gth)*COS(figl-fm)
     &                - .25D0*SIN(tetrc)*COS(fm)
                    ry = d*SIN(gth)*SIN(figl-fm)
     &                - .25D0*SIN(tetrc)*SIN(fm)
                    rz = d*COS(gth) - .25D0*COS(tetrc)
                    rl = SQRT(rx*rx+ry*ry+rz*rz)
                    sf = d*d/rl/rl
                    thc = TACOS(rz/rl)
                    fic = ATAN2(ry,rx)
                    CALL ANGULA(YGP,Idr,1,fi0,fi1,tetrc,thc,fic,ijan,
     &                Op2)
                    DO ixl = 1 , Idr ! For each decay
                      ixm = KSEQ(ixl,3)
                      tfac = TAU(ixm)
C c/s = 0.11991698 /ps, where s = 0.25 cm, c = speed of light in cm/ps
                      IF ( tfac*BETAR(IEXP).GT. 25.D0) THEN
                        WRITE(22,99011) IEXP,KSEQ(ixl,3),tfac
                        IFMO = 0
                      ELSE
                        YGN(ixl) = YGN(ixl) + .11991698D0*tfac*
     &                    BETAR(IEXP)*(sf*YGP(ixl)-YGN(ixl))
                      ENDIF
                    ENDDO ! Loop on decays
                  ENDIF ! If correction due to recoil
                  IF ( IRAWEX(lx).NE.0 ) THEN
                    ipd = ITMA(lx,ijan) ! Get identity of detector
                    DO jyi = 1 , Idr ! For each decay
                      ni = KSEQ(jyi,3)
                      nf = KSEQ(jyi,4)
                      decen = EN(ni) - EN(nf)
                      cocos = SIN(tetrc)*SIN(gth)*COS(fm-figl) +
     &                  COS(tetrc)*COS(gth)
                      decen = decen*(1.+BETAR(lx)*cocos)
                      CALL EFFIX(IEXP,ipd,decen,effi)
                      YGN(jyi) = YGN(jyi)*effi
                    ENDDO
                    inclus = ICLUST(lx,ijan) ! Cluster number for detector ijan
                    IF ( inclus.NE.0 ) THEN
                      DO jyi = 1 , Idr ! For each decay
                        SUMCL(inclus,jyi) = SUMCL(inclus,jyi) + YGN(jyi)
                      ENDDO
                      IF ( ijan.NE.LASTCL(lx,inclus)) GOTO 132 ! If it is not the last detector in the cluster
                      DO jyi = 1 , Idr ! For each decay
                        YGN(jyi) = SUMCL(inclus,jyi)
                      ENDDO
                    ENDIF
                  ENDIF
                  IF ( jfi.EQ.1 ) ija0 = ija0 + 1
                  DO jyi = 1 , Idr ! For each decay
                    GRAD(jyi) = GRAD(jyi) + YGN(jyi)
                  ENDDO ! Loop on decays jyi
                  todfi = todfi + ABS(fi1-fi0)
                ENDDO ! For each phi angle jfi
                IF ( IAXS(lx).EQ.0 ) todfi = pi * 2.D0
                ax = 1.
                IF ( mfla.EQ.1 ) ax = 1./todfi
                dsx = dsig
                IF ( mfla.NE.1 ) dsx = dsig*todfi
                dsxm(mpin,kloop,ktt) = dsx
                WRITE (17,*) lx , mpin , kloop , ktt , dsx
                WRITE (14,*) lx , enb , tting , ija0 , dsx , 
     &            (GRAD(jyi)*dsig*ax,jyi=1,Idr)
                IF ( IPRM(11).EQ.1 ) THEN
                  WRITE (22,99048) lx , ija0 , enb , tta
                  IF ( tta.LT.0. ) WRITE (22,99017) tting
99017             FORMAT (5X,
     &              'RESPECTIVE TARGET SCATTERING ANGLE='
     &              ,1F7.3,1X,'DEG'/)
                  DO jyi = 1 , Idr
                    ni = KSEQ(jyi,3)
                    nf = KSEQ(jyi,4)
                    WRITE (22,99049) ni , nf , SPIN(ni) , SPIN(nf) ,
     &                GRAD(jyi)*dsig*ax , GRAD(jyi)/GRAD(IDRN)
                  ENDDO ! Loop on decays jyi
                ENDIF ! If printout of yields at meshpoints
 132          CONTINUE
            ENDDO ! Loop on detector angles ijan
          ENDDO ! Loop on theta angles ktt
        ENDDO ! Loop on energy meshpoints kloop
      ENDDO ! Loop on pin diodes mpin
      
      EP(lx) = enh
      TLBDG(lx) = tth
      ENDDO ! Loop on experiments lx
      REWIND 14
      REWIND 15
      iske = 0
      DO na = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
        ILE(na) = 1
      ENDDO
      ilx = 0
C     We have now performed the full coulex calculation at each of the
C     meshpoints, so now we start the integration
      DO lx = 1 , NEXPT ! Loop over experiments
C       Read tape 17
        REWIND 17
        DO ijaja = 1 , 300000
          READ (17,*,END=134) jjlx , jmpin , jkloo , jktt , dsx
          IF ( jjlx.EQ.lx ) dsxm(jmpin,jkloo,jktt) = dsx
        ENDDO
 134    na = NANG(lx)
        IF ( lx.NE.1 ) THEN
          DO na1 = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
            ILE(na1) = ILE(na1) + NYLDE(lx-1,na1)
          ENDDO
        ENDIF
        READ (JZB,*) nptx ! Number of meshpoints for stopping powers
        IF ( nptx.NE.0 ) THEN
          READ (JZB,*) (esp(i),i=1,nptx) ! Energy
          READ (JZB,*) (dedx(i),i=1,nptx) ! Stopping power
          npt = nptx
        ENDIF
        READ (JZB,*) npce , npct
        mfla = 0
        IF ( npct.LT.0 ) mfla = 1
        IF ( Iecd(lx).EQ.1 ) mfla = 1
        npct = ABS(npct)
        IF ( npct.GT.100 )
     &    STOP 'ABS(NI2) is limited to 100!'
        npce = npce + MOD(npce,2)
        npct = npct + MOD(npct,2)
        mpin = 1
        IF ( Ipinf.NE.0 ) THEN
          IF ( Jpin(lx).NE.0 ) mpin = Jpin(lx)
        ENDIF
        dst = 0.
        DO lpin = 1 , mpin ! Loop over pin diodes
          ilx = ilx + 1
          IF ( ilx.NE.1 )
     &      CALL TAPMA(lx,iske,isko,iskf,nflr,Idr,0,nft,enb)
          READ (14,*) ne , ntt , emn , emx , tmn , tmx , jan , wth ,
     &      wph , wthh
          iocc = (ne+ntt)*Idr
          IF ( iocc.GT.Izcap ) THEN
            Iretval=1800
            RETURN
          ENDIF
          hen = (emx-emn)/npce
          npce1 = npce + 1
          het = (tmx-tmn)/npct ! Step in theta in degrees
          npct1 = npct + 1
          IF ( Iecd(lx).EQ.1 ) ! Circular detector
     &      CALL COORD(wth,wph,wthh,npct1,1,pfi,wpi,TLBDG(lx),lx,tmn,
     &      tmx)
          IF ( Iecd(lx).NE.1 ) THEN
            IF ( mfla.EQ.1 ) READ (JZB,*)
     &        (pfi(j),j=1,npct1)
          ENDIF
          het = het*pi/180.D0 ! Step in theta in radians
          
C         Interpolate stopping power for each of the energies
C         that we need. esp is an array of energies and dedx is
C         an array containing the stopping powers at those
C         energies. Function is unweighted sqrt. The energies
C         are not the energies we gave for the meshpoints, but
C         the range over which we integrate the bombarding energy
C         with the number of steps specified.
          DO j = 1 , npce1
            xx = (j-1)*hen + emn
            IF ( ISPL.EQ.0 )
     &        CALL LAGRAN(esp,dedx,npt,1,xx,yy,3,1)
            IF ( ISPL.EQ.1 )
     &        CALL SPLNER(esp,dedx,npt,xx,yy,3)
            HLMLM(j) = 1.D0/yy
          ENDDO
          
C         Now we calculate for all the mesh points. 
          naa = NDST(lx)
          IF ( IRAWEX(lx).EQ.0 ) naa = NANG(lx)
          iskf = naa - 1
          DO ja = 1 , naa ! Loop over detector angles
            icll = 3 ! Weighting mode
            DO je = 1 , ne ! ne = number of energy mesh points
              lu = ILE(ja)
              isko = (je-1)*naa*ntt + ja - 1
              CALL TAPMA(lx,iske,isko,iskf,ntt,Idr,1,nft,enb)
              IF ( nft.EQ.1 ) THEN
                Iretval=1900 ! Troubleshoot
                RETURN
              ENDIF
              DO jd = 1 , Idr ! For each decay
                DO jtp = 1 , ntt ! ntt = number of theta meshpoints
                  IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &              DSG(jtp) = dsxm(lpin,je,jtp)
                  jyv = (jtp-1)*Idr + jd
                  YV(jtp) = ZETA(jyv) ! Point yield
                ENDDO ! Loop on theta meshpoints jtp
                DO jt = 1 , npct1 ! number of equal divisions in theta for interpolation
                  xx = (jt-1)*het + tmn*pi/180.D0
                  IF ( ISPL.EQ.0 )
     &              CALL LAGRAN(XV,YV,ntt,jt,xx,yy,2,icll) ! interpolate point yield at theta = xx
                  IF ( ISPL.EQ.1 )
     &              CALL SPLNER(XV,YV,ntt,xx,yy,2) ! interpolate point yield at theta = xx
                  IF ( ISPL.EQ.0 )
     &              CALL LAGRAN(XV,DSG,ntt,jt,xx,zz,2,icll) ! interpolate gamma yield at theta = xx
                  IF ( ISPL.EQ.1 )
     &              CALL SPLNER(XV,DSG,ntt,xx,zz,2) ! interpolate gamma yield at theta = xx
                  IF ( mfla.EQ.1 ) yy = yy*pfi(jt)*pi/180.D0
                  IF ( yy.LE.0. ) yy = 1.D-15
                  IF ( mfla.EQ.1 ) zz = zz*pfi(jt)*pi/180.D0
                  XI(jt) = yy*SIN(xx) ! yy = integral of point yields over phi
                  IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &              = zz*SIN(xx) ! zz = integral over phi of Rutherford cross section
                ENDDO ! Loop on equal theta divisions jt
                icll = 4
                locat = ntt*Idr + (je-1)*Idr + jd
C               Integrate point yields over theta using Simpson's rule
                ZETA(locat) = SIMIN(npct1,het,XI)
C               If it is first decay and angle, integrate Rutherford cross section over theta
                IF ( jd.EQ.1 .AND. ja.EQ.1 ) DSE(je)
     &            = SIMIN(npct1,het,HLM)
                ZV(je) = enb
              ENDDO ! Loop on decays jd
            ENDDO ! Loop on energy meshpoints je
            
C    Interpolation over energy:
C    The array ZV contains the energies of the meshpoints and the elements of the YV
C    array are set to the angle-integrated yield for each decay at the corresponding
C    energy, while DSE contains the Rutherford cross section for those energies. Since
C    the energies of the meshpoints are not necessarily equally spaced, we need to
C    interpolate to a set of equally spaced energies separated by "hen" starting from
C    "emn". To get the contribution from each energy, dE = 1 / (stopping power). Note
C    that we only evaluate the Rutherford cross section for the first decay and first
C    angle, since it is the same for all.

            icll = 3
            DO jd = 1 , Idr ! For each decay
              DO jtp = 1 , ne ! For each energy meshpoint
                jyv = (jtp-1)*Idr + jd + ntt*Idr
                YV(jtp) = ZETA(jyv)
              ENDDO ! Loop on energy meshpoints jtp
              DO jt = 1 , npce1 ! npce1 is number of equal energy steps
                xx = (jt-1)*hen + emn
                
C               Interpolate the angle-integrated yield for this energy
                IF ( ISPL.EQ.0 )
     &            CALL LAGRAN(ZV,YV,ne,jt,xx,yy,2,icll)
                IF ( ISPL.EQ.1 )
     &            CALL SPLNER(ZV,YV,ne,xx,yy,2)
                
C               Interpolate Rutherford cross-section for this energy
                IF ( jd.EQ.1 .AND. ja.EQ.1 .AND. ! Only for first decay and angle
     &            ISPL.EQ.0 )
     &            CALL LAGRAN(ZV,DSE,ne,jt,xx,zz,2,icll) ! Interpolate for this energy
                IF ( jd.EQ.1 .AND. ja.EQ.1 .AND. ISPL.EQ.1 )
     &            CALL SPLNER(ZV,DSE,ne,xx,zz,2) ! Interpolate for this energy
                IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &            = zz*HLMLM(jt) ! HLMLM = 1 / stopping power
                XI(jt) = yy*HLMLM(jt)
              ENDDO ! Loop on equal energy steps
              
C   So now after this loop, we have XI containing the angle-integrated yield times dE for 
C   a set of equally spaced energies, so we use Simpson's rule to integrate them and store
C   in GRAD(jd). The first time, we also have in HLM a set of Rutherford cross-sections for
C   equally spaced energies, which we integrate in the same way.
              icll = 4
              IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &          DS = SIMIN(npce1,hen,HLM) ! integrate
              GRAD(jd) = SIMIN(npce1,hen,XI)
            ENDDO ! Loop over decays jd
            
            IF ( ja.EQ.1 ) dst = dst + DS
            IF ( ja.EQ.1 ) WRITE (22,99018) DS , lx
99018       FORMAT (1X/////5X,
     &        'INTEGRATED RUTHERFORD CROSS SECTION='
     &        ,1E9.4,2X,'FOR EXP.',1I2///)
            
            WRITE (22,99019) lx , ja , emn , emx , tmn , tmx
99019       FORMAT (1X,//50X,'INTEGRATED YIELDS'//5X,
     &        'EXPERIMENT ',1I2,2X,'DETECTOR ',
     &        1I2/5X,'ENERGY RANGE ',1F8.3,'---',
     &        1F8.3,1X,'MEV',3X,
     &        'SCATTERING ANGLE RANGE ',1F7.3,
     &        '---',1F7.3,1X,'DEG'//5X,'NI',5X,
     &        'NF',5X,'II',5X,'IF',5X,'YIELD',5X,
     &        'NORMALIZED YIELD'/)
            DO jd = 1 , Idr
              WRITE (15,*) GRAD(jd)
            ENDDO
            DO jd = 1 , Idr
              ni = KSEQ(jd,3)
              nf = KSEQ(jd,4)
              WRITE (22,99049) ni , nf , SPIN(ni) , SPIN(nf) ,
     &          GRAD(jd) , GRAD(jd)/GRAD(IDRN) ! IDRN is the normalising transition
            ENDDO
          ENDDO ! Loop over detector angles ja
          
          IF ( Iecd(lx).EQ.1 ) THEN ! Circular detector
            IF ( Jpin(lx).EQ.0 ) THEN
              CALL COORD(wth,wph,wthh,1,2,pfi,wpi,TLBDG(lx),lx,txx,txx)
              WRITE (22,99020) FIEX(lx,1)*180.D0/pi ,
     &          FIEX(lx,2)*180.D0/pi , lx
99020         FORMAT (//5X,'WARNING: THE PHI ANGLE WAS REPLACED BY'
     &          ,1X,F8.3,1X,'TO',F8.3,3X,'FOR EXPERIMENT',2X,I3)
              IF ( TLBDG(lx).LT.0 ) THEN
                FIEX(lx,1) = FIEX(lx,1) + pi
                FIEX(lx,2) = FIEX(lx,2) + pi
              ENDIF ! If theta_lab < 0
            ENDIF ! If no pin diodes
          ENDIF ! If circular detector
          iske = iske + ne*ntt*naa
        ENDDO ! Loop over pin diodes
        IF ( mpin.GT.1 ) WRITE (22,99021) dst , lx
99021   FORMAT (1x//2x,
     &    'Total integrated Rutherford cross section='
     &    ,1E8.3,' for exp. ',1I2/)
      ENDDO
      REWIND 17 ! Added PJN (17Jul2009)
      IF ( Ipinf.NE.0 ) THEN
        ngpr = 0
        DO lx = 1 , NEXPT ! For each experiment
          nged = NDST(lx)
          IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
          IF ( lx.NE.1 ) ngpr = ngpr + Idr*Jpin(lx-1)*NDST(lx-1)
          lpin = Jpin(lx)
          IF ( lpin.EQ.0 ) lpin = 1
          DO jgd = 1 , nged ! For each angle or dataset
            DO jd = 1 , Idr
              GRAD(jd) = 0.
            ENDDO
            DO mpin = 1 , lpin ! For each pin diode
              REWIND 15
              ndum = ngpr + (jgd-1)*Idr + (mpin-1)*nged*Idr ! Was jgd instead of nged (PJN 17Jul2009)
              IF ( ndum.NE.0 ) THEN
                DO jd = 1 , ndum
                  READ (15,*) xx
                ENDDO
              ENDIF
              DO jd = 1 , Idr ! For each decay
                READ (15,*) xx
                GRAD(jd) = GRAD(jd) + xx
              ENDDO ! Loop on decays jd
            ENDDO ! Loop on pin diodes mpin
            WRITE (17,*) (GRAD(jd),jd=1,Idr)
          ENDDO ! Loop on angle or dataset jgd
        ENDDO ! Loop on experiment lx
        REWIND 15
        REWIND 17
        DO lx = 1 , NEXPT ! For each experiment
          nged = NDST(lx)
          IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
          DO ija0 = 1 , nged ! For each angle or dataset
            READ (17,*) (GRAD(jdy),jdy=1,Idr)
            DO jd = 1 , Idr ! For each decay
              WRITE (15,*) GRAD(jd)
            ENDDO ! Loop on decays jd
          ENDDO ! Loop on angle or dataset ija0
        ENDDO ! Loop on experiments lx
      ENDIF
      DO lx = 1 , NEXPT ! For each experiment restore original ISKIN
        ISKIN(lx) = iskin_protect(lx)
      ENDDO
99011 FORMAT (1X,/,2X,'DURING THE INTEGRATION',1X,
     &  'IT WAS NECESSARY TO SWITCH OFF THE TIME-OF-FLIGHT CORRECTION '
     &  'FOR EXPT ',I3,' LEVEL ', I3, ' TAU=',1E9.4)
99048 FORMAT (1X//50X,'CALCULATED YIELDS'//5X,'EXPERIMENT ',1I2,2X,
     &        'DETECTOR ',1I2/5X,'ENERGY ',1F10.3,1X,'MEV',2X,'THETA ',
     &        1F7.3,1X,'DEG'//5X,'NI',5X,'NF',5X,'II',5X,'IF',5X,
     &        'YIELD',5X,'NORMALIZED YIELD'/)
99049 FORMAT (4X,1I3,4X,1I3,3X,1F4.1,3X,1F4.1,3X,1E11.5,3X,1E11.5)
      END
