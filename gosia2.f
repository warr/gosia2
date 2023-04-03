      INCLUDE 'header.txt'
C
C PROGRAM GOSIA2
C
C Calls: ADHOC, ALLOC, ANGULA, ARCCOS, ARCTG, CMLAB, DECAY, DJMM,
C        EFFIX, ELMT, FAKP, FHIP, FTBM, INTG, INVKIN, KLOPOT, KONTUR,
C        LOAD, MINI, MIXR, MIXUP, OPENF, PATH, PRELM, PTICC, QFIT, READY,
C        SETIN, SNAKE, STING, TACOS, TEMB, TENS, WSIXJ,
C        WTHREJ
C
C Uses global variables:
C      ABC    - absorption coefficients
C      ACCA   - accuracy
C      ACCUR  - accuracy required
C      AGELI  - angles of Ge detectors
C      AKAVKA - efficiency curve parameters
C      ARM    - excitation amplitudes of substates.
C      AVJI   - average J (N.B. here it is G(1))
C      B      - table of factorials
C      BEQ    - identifier for angle for rotations
C      BETAR  - recoil beta
C      CAT    - substates of levels (n_level, J, m)
C      CC     - conversion coefficients
C      CNOR   - normalization factors
C      CNOR1  - normalization factors for target (for gosia2)
C      CNOR2  - normalization factors for projectile (for gosia2)
C      CORF   - internal correction factors
C      DEVD   -
C      DEVU   -
C      DIPOL  - E1 polarization parameter
C      DIX    - Ge parameters (inner & outer radius, length, distance)
C      DLOCK  - limit derivative below which matrix element is fixed if LOCKS=1
C      DS     - integrated rutherford cross-section
C      DSE    - rutherford cross section at given energy integrated over angles
C      DSG    - differential gamma-ray yield at meshpoints
C      DSIGS  - dsigma for each experiment
C      DYEX   - error on experimental yield
C      EAMX   - known matrix elements and their errors
C      ELM    - matrix elements
C      ELM25  - matrix elements storage for particle1 (for gosia2)
C      ELM26  - matrix elements storage for particle2 (for gosia2)
C      ELMH   -
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      EMMA   - Controls number of magnetic substates in full coulex calc.
C      EN     - energy of level
C      EP     - bombarding energy
C      ERR    - error flag
C      EXPO   - adiabatic exponential
C      FIEL   - K (N.B. here it is G(6))
C      FIEX   - phi range of particle detector
C      GAMMA  - Gamma (N.B. here it is G(2))
C      GFAC   - g (N.B. here it is G(5))
C      GRAD   - partial derivative of chi squared wrt. each matrix element
C      HLM    - matrix elements before minimisation
C      HLMLM  - old value of matrix element or chi squared
C      IAMX   - index of matrix element for known matrix element
C      IAMY   - level indices of pair of levels for which matrix element is known
C      IAX    - axial symmetry flag
C      IBPS   - target/projectile switch (gosia2)
C      IBYP   - flag to indicate whether we calculate <\alpha_k>
C      ICLUST - cluster number for each experiment and detector
C      ICS    - read internal correction factors flag (OP,CONT switch CRF,)
C      IDIVE  - number of subdivisions
C      IDRN   - index of normalising transition for yields
C      IEXP   - experiment number
C      IFAC   - spin/parity phase factor
C      IFBFL  - calculate derivatives with forward-backward method
C      IFMO   - include correction to angular distance for finite recoil distance.
C      ILE    - yield number for each detector
C      IMIN   -
C      INHB   - inhibit error flag (LERF) setting in POMNOZ
C      INNR   - independent normalisation switch (see OP,CONT INR,)
C      INTERV - default accuracy check parameter for Adams-Moulton (see OP,CONT:INT)
C      INTR   - flag to swap chisqr and log(chisqr)
C      IP     - table of prime numbers
C      IPRM   - various flags to control output
C      IPS1   - terminate after calculating and storing internal correction factors
C      IRAWEX - flag to indicate raw uncorrected yield
C      ISEX   -
C      ISKIN  - kinematic flag (0,1)
C      ISMAX  - number of substates used
C      ISO    - isotropic flag
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      ITTE   - thick target experiment flag
C      IUNIT3 - unit for TAPE3
C      IVAR   - indicates a limit or correlation is set
C      IWF    - warning flag
C      IY     - index for yields
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of non-investigated nucleus
C      JENTR  - flag set to 0 normally, 1 in OP,ERRO
C      JSKIP  - Experiments to skip during minimisation.
C      JZB    - unit to read from
C      KFERR  - error flag for minimization
C      KSEQ   - index of level
C      KVAR   -
C      LAMAX  - number of multipolarities to calculate
C      LAMBDA - list of multipolarities to calculate
C      LASTCL - index of last detector in cluster
C      LDNUM  - number of matrix elements with each multipolarity populating each level
C      LEAD   - pair of levels involved in each matrix element
C      LIFCT  - index for lifetimes
C      LMAX   - ground-state spin + 1
C      LMAXE  - maximum multipolarity needed for calculation
C      LNORM  - normalisation constant control
C      LNY    - use logs to calculate chi squared
C      LOCKF  - flag to fix matrix elements with most significant derivative
C      LOCKS  - lock flag. If LOCKS=1, fix at first stage of minimization
C      LP1    - maximum number of experiments (50)
C      LP10   - maximum number of substates (1200)
C      LP11   - LP8 - 1 (2800)
C      LP12   - number of steps of omega (365)
C      LP13   - LP9 + 1 (47901)
C      LP14   - maximum space for collision functions (4900)
C      LP2    - maximum number of matrix elements (1500)
C      LP3    - maximum number of levels (100)
C      LP4    - maximum number of yields (1500)
C      LP6    - maximum number of gamma detectors (32)
C      LP7    - start of collision functions (45100)
C      LP8    - (104)
C      LP9    - length of ZETA - 2100 (47900)
C      MAGA   - number of magnetic substates in approximate calculation
C      MAGEXC - flag: 0 means no magnetic excitations, 1 means with mag. exc.
C      MCFIX  - flag to say if normalization coefficient is fixed (for gosia2)
C      MEMAX  - number of matrix elements
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      MULTI  - number of matrix elements having given multipolarity
C      NAMX   - number of known matrix elements
C      NANG   - number of gamma-ray detectors for each experiment
C      NBRA   - number of branching ratios
C      NCM    - calculate kinematics assuming this state for final state (default = 2)
C      NDIM   - maximum number of levels
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NLIFT  - number of lifetimes
C      NLOCK  - number of elemnts to fix if LOCKF=1
C      NMAX   - number of levels
C      NMAX1  - number of levels with decays
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      PARX   - [for maps]
C      PARXM  - [for maps]
C      POWER  - x (N.B. here it is G(7))
C      QAPR   - approximate Coulomb amplitudes
C      SA     - ratio of matrix elements for correlated elements
C      SE     - seed for random number generator of OP,RAND
C      SGW    - number of standard deviations to generate warning (see control option WRN,X)
C      SPIN   - spin of level
C      SUBCH1 - partial chisqr
C      SUBCH2 - partial chisqr
C      SUMCL  - sum of yields for clusters
C      TAU    - lifetime in picoseconds
C      THICK  - thickness of each absorber type
C      TIMEC  - Tau_C (N.B. here it is G(4))
C      TIMEL  - lifetimes and their errors
C      TLBDG  - theta of particle detector in lab frame (in degrees)
C      TREP   - theta of recoiling nucleus (in radians)
C      UPL    - upper limits for all gamma detectors
C      VINF   - speed of projectile at infinty
C      WNOR   - weighting factor for the target normalisation (0.5 = equal weight, 1.0 = all weight to target)
C      XA     - A of investigated nucleus
C      XA1    - A of non-investigated nucleus
C      XI     - xi coupling coefficients
C      XIR    - [for maps]
C      XLAMB  - Lambda* (N.B. here it is G(3))
C      XV     - energy meshpoints (sometimes theta meshpoints) where we calculate exact Coulex
C      YEXP   - experimental yields
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      YGP    - gamma yield calculated with correction to angular distribution from finite recoil distance
C      YNRM   - relative normalisation for gamma detectors
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C      ZETA   - various coefficients
C      ZPOL   - dipole term (GDR excitation)
C      ZV     - energy meshpoints

      PROGRAM GOSIA2
      IMPLICIT NONE
      REAL*8 acof , ap , ARCCOS , ARCTG , arg , bcof , be2 , be2a ,
     &       be2b , be2c , bl , bmx , bten , bu , ccc , ccd ,
     &       cf , chilo , chiok , chis0 , chisl , chisq , chiss , cnst ,
     &       cocos , conu , d , decen , dsd , dsig , effi , eh1 , elmi ,
     &       emhl1 , esd , ess , fi0 , fi1 , fic , fiex1 , figl ,
     &       fipo1 , fm , gth , p , ph1 , ph2 , po1 , po2 , polm ,
     &       pop1 , pr , pv , q1 , q2 , qc , qfac , qr , qui , r , r1 ,
     &       r2 , r3 , r4 , rem , remax , rl , rlr , rm , rx , ry , rz ,
     &       s , s11 , s12 , s21 , s22 , sbe , sf , sh , sh1 , sh2 ,
     &       slim , summm , sz1 , sz2 , TACOS , test , tfac , thc ,
     &       title , ttttt , txx , u , val , waga , WSIXJ , WTHREJ ,
     &       xep , xl1 , xlk , xtest , xw , xxi , ycorr , yyd1 , yydd ,
     &       yyy , zmir , zp , zz
      INTEGER*4 i , i122 , iapx , icg , ict , ictl , id , ideff , idf ,
     &          idr , iecd , ient , ifbp , ifc , ifm , ifwd , ig1 ,
     &          ig2 , ih1 , ih2 , ihlm , ihuj , ii , ij , ijk , ijx ,
     &          ile1 , ilx , im , imode , in1 , in2 , inclus , ind ,
     &          indx , inko , inm1 , inm2 , inn , inpo , intend ,
     &          intvh , inx1 , iobl , iopri , iosr , ipd , iph , ipine ,
     &          ipinf , ipo1 , ipo2 , ipo3 , ipp , iprc , ipri , irea ,
     &          irep , iretval , irfix , irix , isip , iskok , isoh ,
     &          ispa , ispb , itno , itp , iuy , iva , iva1 , ivarh ,
     &          ivari , ivrh , ixj , ixl , ixm , iyr , izcap , j , jd ,
     &          jde , jex , jexp , jfre , jgl , jgl1 , jgr , jgs , jj ,
     &          jj1 , jjjj , jjx , jk , jl , jmm , jp , jphd , jpin ,
     &          jrls , js , jyi , jyi1 , jyi2 , jz , k , kclust , kerf ,
     &          kex , kh , kh1 , kh2 , kk , kk1 , kk2 , kkk , kmat ,
     &          kq , kuku , l , la , la1 , lam , lamh , lck1 , lck2 ,
     &          levmax , lex , lexp , lfagg , lfini , lh1 , lh2 ,
     &          liscl , lkj , lkj1 , ll , lli , lll , lmax1 , lmaxh ,
     &          loct , lp0 , ltrn , ltrn1 , ltrn2 , lu , lxd , magh ,
     &          MEM , memax1 , memh , memx4 , mend , mexl , mm , ms ,
     &          n , na , nallow , naxfl , nch , ndima , ne , nf , nfd ,
     &          nfdd , ni , nksi , nl , nmaxh , nmemx , nogeli , nptl ,
     &          ns1 , ns2 , ntap , numcl , nval , nz
      CHARACTER*4 oph , op1 , opcja , op2
      CHARACTER*80 line
      CHARACTER*1 prp
      REAL*8 chir , chp , ccch1 ,ccch2 , cht ! For gosia2
      INTEGER*4 nawr , mmmm , kkkk , mres1 , mrepf  ! For gosia2
      INTEGER*4 mret , mawr ! For gosia2
      DIMENSION ihlm(32) , bten(1600) , ! bten dimension = 16 * maxlevels
     &          fiex1(100,100,2) , title(20) , zmir(6,2,50) ,
     &          iecd(50) , xl1(7) , qui(8,10) , cf(8,2) , ivarh(1500) ,
     &          liscl(200) , ivari(1500) , jpin(50) , ideff(50)
      INCLUDE 'clust.inc'
      INCLUDE 'cccds.inc'
      INCLUDE 'inhi.inc'
      INCLUDE 'ident.inc'
      INCLUDE 'efcal.inc'
      INCLUDE 'tcm.inc'
      INCLUDE 'brec.inc'
      INCLUDE 'adbxi.inc'
      INCLUDE 'dimx.inc'
      INCLUDE 'tra.inc'
      INCLUDE 'cinit.inc'
      INCLUDE 'xra.inc'
      INCLUDE 'hhh.inc'
      INCLUDE 'vac.inc'
      INCLUDE 'me2d.inc'
      INCLUDE 'life1.inc'
      INCLUDE 'dftb.inc'
      INCLUDE 'erran.inc'
      INCLUDE 'mgn.inc'
      INCLUDE 'seck.inc'
      INCLUDE 'vlin.inc'
      INCLUDE 'dumm.inc'
      INCLUDE 'brnch.inc'
      INCLUDE 'yexpt.inc'
      INCLUDE 'yteor.inc'
      INCLUDE 'lev.inc'
      INCLUDE 'map.inc'
      INCLUDE 'ccc.inc'
      INCLUDE 'ggg.inc'
      INCLUDE 'az.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'cxi.inc'
      INCLUDE 'clcom.inc'
      INCLUDE 'coex.inc'
      INCLUDE 'minni.inc'
      INCLUDE 'cx.inc'
      INCLUDE 'cexc.inc'
      INCLUDE 'prt.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'cb.inc'
      INCLUDE 'clm.inc'
      INCLUDE 'clcom0.inc'
      INCLUDE 'clcom8.inc'
      INCLUDE 'clcom9.inc'
      INCLUDE 'comme.inc'
      INCLUDE 'coex2.inc'
      INCLUDE 'cexc9.inc'
      INCLUDE 'caux0.inc'
      INCLUDE 'pth.inc'
      INCLUDE 'aprcat.inc'
      INCLUDE 'warn.inc'
      INCLUDE 'thtar.inc'
      INCLUDE 'fit.inc'
      INCLUDE 'aprx.inc'
      INCLUDE 'skp.inc'
      INCLUDE 'trb.inc'
      INCLUDE 'sel.inc'
      INCLUDE 'ercal.inc'
      INCLUDE 'logy.inc'
      INCLUDE 'fakul.inc'
      INCLUDE 'life.inc'
      INCLUDE 'switch.inc'
      INCLUDE 'fconst.inc'
      INCLUDE 'cinit2.inc' ! For gosia2
      INCLUDE 'resc.inc' ! For gosia2
      INCLUDE 'nist.inc'
      DATA q1/0./,q2/0./,iph/0/

C     Initialize prime numbers
      IP(1) = 2
      IP(2) = 3
      IP(3) = 5
      IP(4) = 7
      IP(5) = 11
      IP(6) = 13
      IP(7) = 17
      IP(8) = 19
      IP(9) = 23
      IP(10) = 29
      IP(11) = 31
      IP(12) = 37
      IP(13) = 41
      IP(14) = 43
      IP(15) = 47
      IP(16) = 53
      IP(17) = 59
      IP(18) = 61
      IP(19) = 67
      IP(20) = 71
      IP(21) = 73
      IP(22) = 79
      IP(23) = 83
      IP(24) = 89
      IP(25) = 97
      IP(26) = 101

C     Initialize pointers
      lp0 = 155600 ! Size of ZETA array
      LP1 = 50 ! Maximum number of experiments
      LP2 = 1500 ! Maximum number of matrix elements
      LP3 = 100 ! Maximum number of levels
      LP4 = 1500
      LP6 = 32 ! Maximum number of gamma detectors
      LP7 = lp0 - 4900 ! Start of collision coefficients in ZETA
      LP8 = LP3*28 + 1
      LP9 = lp0 - LP3*28
      LP10 = 1200 ! Maximum number of substates
      LP11 = LP8 - 1
      LP12 = 365 ! Maximum number of steps of omega (dimension of ADB, SH, CH)
      LP13 = LP9 + 1
      LP14 = 4900 ! Maximum number of collision coefficients

      JZB = 5
      idr = 0

C---- gosia2 changes start
C     Read target/projectile switch and first option from standard input
      READ (*,*) IBPS
      READ (*,99001) op1 , op2

C     Open files
      CALL OPENF
      nawr = 0

C     Use input unit 25 for particle 1 and 26 for particle 2
      IBPS = IBPS - 1
      JZB = 25
      IF ( IBPS.EQ.1 ) JZB = 26
C---- gosia2 changes end

C     Initialize normalization to 1. and weighting to 0.5
      WNOR = 0.5D0
      DO i = 1 , LP3 ! LP3 = 100 (maximum number of levels)
         DO j = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
            CNOR(j,i) = 1.
            CNOR1(j,i) = 1. ! Added for gosia2
            CNOR2(j,i) = 1. ! Added for gosia2
         ENDDO
      ENDDO

C---- gosia2 changes start
      IUNIT3 = 33 ! Is 3 in gosia
      MCFIX = 1
      mres1 = 0
      mrepf = 0
      mret = -1
      chir = 0.
      mawr = 0

 2200 DO i = 1 , 50 ! This is where we come back when we iterate over target/beam
        LIFCT(i) = 0
      ENDDO

      chp = chir
      IF ( mret.EQ.1 .AND. JZB.EQ.26 ) chir = 0.
      REWIND 25 ! Rewind input files
      REWIND 26
      READ (JZB,*) IBPS
      IBPS = IBPS - 1
      JZB = 25
      IF ( IBPS.EQ.1 ) JZB = 26
C---- gosia2 changes end

      IBYP = 0
      INHB = 0
      BEQ = -983872.D0 ! This is just a magic value to indicate it was not
                       ! yet initialised
      ipinf = 0
      iyr = 0
      INNR = 0
      itno = 0
      chisq = 0.
      chilo = 0.
      IWF = 1 ! Turn on warnings
      ifm = 0 ! Fast minimisation switch off by default
      IPS1 = 11
      ifwd = -1
      INTR = 0
      LNY = 0
      JENTR = 0 ! Flag to indicate we are not in OP,ERRO
      ICS = 0
      ISPL = 0 ! Flag to indicate we should use LAGRAN not SPLNER

      DO i = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         ideff(i) = 0
         jpin(i) = 0
         iecd(i) = 0
      ENDDO
      txx = 0.
      SGW = 3.D0
      SUBCH1 = 0.
      SUBCH2 = 0.
      ITS = 0 ! Create tape 18 flag
      iosr = 0
      LOCKS = 0
      DLOCK = 1.1D0
      kerf = 0
      IFBFL = 0
      NLOCK = 0
      LOCKF = 0
      DO i = 1 , LP4 ! LP4 = 1500
         DO j = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
            CORF(i,j) = 1.D0
         ENDDO
      ENDDO
      DO i = 1 , 20
         IPRM(i) = 1
      ENDDO
      DO i = 1 , 50
         DO j = 1 , 5
            CC(i,j) = 0.
         ENDDO
      ENDDO
      IPRM(4) = -2
      IPRM(5) = 11111
      IPRM(6) = 11111
      IPRM(7) = 0
      IPRM(16) = 0
      IPRM(17) = 0
      IPRM(18) = 0
      IPRM(19) = 0
      IPRM(20) = 0
      IPRM(21) = 0
      DO i = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         DO j = 1 , 5
            IF ( j.NE.5 ) THEN
               DO k = 1 , 10
                  DO kuku = 1 , 6
                     PARXM(i,j,k,kuku) = 0.
                  ENDDO
               ENDDO
            ENDIF
            DO k = 1 , 12
               PARX(i,k,j) = 0.
            ENDDO
         ENDDO
      ENDDO
      DO k = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         IDIVE(k,1) = 1
         IDIVE(k,2) = 1
         DO iuy = 1 , 6
            XIR(iuy,k) = 0.
         ENDDO
      ENDDO
      iobl = 0
      lfagg = 0
      izcap = 12800
      KFERR = 0
      NDIM = LP3 ! LP3 = 100 (maximum number of levels)
      ISO = 1
      B(1) = 1.D0
      DO i = 2 , 20
         B(i) = B(i-1)*(i-1)
      ENDDO
      LMAXE = 0
      CALL FAKP
      CALL FHIP
      NCM = 2 ! Default final state for kinematics calculation (OP,CONT NCM,)
      DO ijx = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         INTERV(ijx) = 1000
      ENDDO
      la = 0
      ipo3 = 1
      indx = 0
      ACCUR = .00001D0
      icg = 1
      ient = 1
      jphd = 1 ! Print header flag
      DIPOL = 0 ! WAS 0.005 (see function LOAD, where this is now set)
      MAGEXC = 0 ! Initially flag that we don't need magnetic excitations
      LAMMAX = 0
      DO lam = 1 , 8
         DO lexp = 1 , LP3 ! LP3 = 100 (maximum number of levels)
            LDNUM(lam,lexp) = 0
         ENDDO
         MULTI(lam) = 0
         LAMDA(lam) = 0
      ENDDO
      DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
         EXPO(j) = (1.,0.)
         KVAR(j) = 1
         ELM(j) = 0.
      ENDDO
      DO j = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
         JSKIP(j) = 1
         ISKIN(j) = 0
      ENDDO
      DO j = 1 , LP3 ! LP3 = 100 (maximum number of levels)
         ISEX(j) = 1111
      ENDDO
      ISEX(1) = 0
      ACCA = .00001D0
      oph = '    '
      nmemx = LP2 + 9 ! LP2 = 1500 (maximum number of matrix elements)
      IEXP = 1
      IMIN = 0
      i122 = 0
      DO j = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
         DO k = 1 , 2
            DO l = 1 , 7
               QAPR(j,k,l) = 0.
            ENDDO
         ENDDO
         ivari(j) = 0
      ENDDO
      ERR = .FALSE.
      opcja = '    '
      levmax = 0
      intend = 0 ! End of initialization

C.............................................................................
C     Start reading input file.
 100  READ (JZB,99001) op1 , op2
99001 FORMAT (1A3,1A4)

C     If line doesn't start with OP,XXXX, we have a syntax error, so stop
      IF ( op1.NE.'OP, ' ) THEN
         WRITE(*,*) 'Line does not have the OP,XXXX syntax: ', op1, op2
         STOP 'SYNTAX ERROR'
      ENDIF

      IF ( op2.EQ.'GOSI' ) THEN
        oph = op2
        opcja = op2
      ENDIF

C     Treat OP,FILE (attach files to fortran units)
      IF ( op2.EQ.'FILE' ) THEN
        CALL OPENF
        GOTO 100 ! End of OP,FILE - back to input loop
      ENDIF

C     Print header
      IF ( jphd.EQ.1 ) WRITE (22,99002)
99002 FORMAT ('1'/1X,125('*')/1X,125('*')/1X,50('*'),25X,50('*')/1X,
     &  50('*'),10X,'GOSIA2',9X,50('*')/1X,50('*'),25X,50('*')
     &  /1X,125('*')/1X,125('*')////)
      IF ( jphd.EQ.1 ) WRITE (22,99003)
99003 FORMAT (1X/20X,'ROCHESTER COULOMB EXCITATION DATA ANALYSIS ',
     &  'CODE BY T.CZOSNYKA,D.CLINE AND C.Y.WU'/50X,
     &  'LATEST REVISION- JUNE  2006'//////)
      jphd = 0 ! Set print header flag to zero, so we don't repeat header

C     Treat OP,BRIC
      IF ( op2.EQ.'BRIC' ) THEN
        GOTO 70000

C     Treat OP,CORR
      ELSE IF ( op2.EQ.'CORR' ) THEN
        GOTO 70001

C     Treat OP,COUL
      ELSE IF ( op2.EQ.'COUL' ) THEN
        GOTO 70002

C     Treat OP,ERRO (calculate errors)
      ELSE IF ( op2.EQ.'ERRO' ) THEN
        GOTO 70003

C     Treat OP,EXIT
      ELSE IF ( op2.EQ.'EXIT' ) THEN
        GOTO 70004

C     Treat OP,GDET (germanium detectors)
      ELSE IF ( op2.EQ.'GDET' ) THEN
        GOTO 70005

C     Treat OP,GOSI
      ELSE IF ( op2.EQ.'GOSI' ) THEN
        GOTO 70002

C     Treat OP,INTG
      ELSE IF ( op2.EQ.'INTG' ) THEN
        GOTO 70006

C     Treat OP,INTI
      ELSE IF ( op2.EQ.'INTI' ) THEN
        GOTO 70007

C     Treat OP,MAP
      ELSE IF ( op2.EQ.'MAP ' ) THEN
        iobl = 1
        GOTO 70008 ! End of OP,MAP

C     Treat OP,MINI
      ELSE IF ( op2.EQ.'MINI' ) THEN
        GOTO 70009

C     Treat OP,POIN
      ELSE IF ( op2.EQ.'POIN' ) THEN
        GOTO 70008

C     Treat OP,RAND (randomise matrix elements)
      ELSE IF ( op2.EQ.'RAND' ) THEN
        GOTO 70010

C     Treat OP,RAW (raw uncorrected gamma yields)
      ELSE IF ( op2.EQ.'RAW ' ) THEN
        GOTO 70011

C     Treat OP,RE,A (release A)
      ELSE IF ( op2.EQ.'RE,A' ) THEN
        jfre = 0
        irfix = 0
        GOTO 70012

C     Treat OP,RE,C (release C)
      ELSE IF ( op2.EQ.'RE,C' ) THEN
        jfre = 1
        irfix = 0
        GOTO 70012 ! End of OP,RE,C

C     Treat OP,RE,F (release F)
      ELSE IF ( op2.EQ.'RE,F' ) THEN
        jfre = 0
        irfix = 1
        GOTO 70012

C     Treat OP,REST (restart)
      ELSE IF ( op2.EQ.'REST' ) THEN
        GOTO 70013

C     Treat OP,SELE
      ELSE IF ( op2.EQ.'SELE' ) THEN
        GOTO 70014

C     Treat OP,SIXJ
      ELSE IF ( op2.EQ.'SIXJ' ) THEN
        GOTO 70015

C     Treat OP,STAR
      ELSE IF ( op2.EQ.'STAR' ) THEN
        GOTO 70008

C     Treat OP,THEO
      ELSE IF ( op2.EQ.'THEO' ) THEN
        GOTO 70016

C     Treat OP,TITL (title)
      ELSE IF ( op2.EQ.'TITL' ) THEN
        GOTO 70017

C     Treat OP,TROU (troubleshooting)
      ELSE IF ( op2.EQ.'TROU' ) THEN
        GOTO 70018

C     Treat OP,YIEL
      ELSE IF ( op2.EQ.'YIEL' ) THEN
        GOTO 70019

      ENDIF ! End cases for different options

      WRITE (22,99022) op1 , op2
99022 FORMAT (5X,'UNRECOGNIZED OPTION',1X,1A3,1A4)
      GOTO 2000 ! Normal end of execution
C.............................................................................

C     Handle OP,ERRO
 400  IF ( ICS.EQ.1 ) THEN
         REWIND 11
         DO kh1 = 1 , LP4
            READ (11) (CORF(kh1,kh2),kh2=1,LP6) ! LP6 = 32 (maximum number of gamma detectors)
         ENDDO
      ELSE
         CALL FTBM(0,chiss,idr,0,chilo,bten)
         REWIND 11
         DO kh1 = 1 , LP4
            WRITE (11) (CORF(kh1,kh2),kh2=1,LP6) ! LP6 = 32 (maximum number of gamma detectors)
         ENDDO
      ENDIF

      CALL FTBM(3,chiss,idr,1,chilo,bten)
      chis0 = chiss
      WRITE (22,99028) chis0
99028 FORMAT (1X///10X,'***** CENTRAL CHISQ=',1E12.4,1X,'*****'//)
      INHB = 1
      chisl = chiss
      DO kh = 1 , MEMAX
         HLM(kh) = ELM(kh)
      ENDDO
      IF ( idf.EQ.1 ) THEN
         IFBFL = 1
         IF ( irep.NE.2 ) GOTO 700
         IF ( iosr.EQ.0 ) GOTO 700
         REWIND IUNIT3
         READ (IUNIT3,*) ll , mm , kk , inn
         DO inn = 1 , ll
            READ (IUNIT3,*) mm , yyy , zz
         ENDDO
         DO inn = 1 , MEMAX
            READ (IUNIT3,*) mm , ll , kk
         ENDDO
         DO inn = 1 , MEMAX
            READ (IUNIT3,*) mm , yyy
         ENDDO
 450     READ (IUNIT3,*) mm , ll
         IF ( mm.EQ.0 ) THEN
            BACKSPACE IUNIT3
            GOTO 700
         ELSE
            READ (IUNIT3,*) kk , ll , yyy
            READ (IUNIT3,*) (SA(mm),mm=1,MEMAX)
            GOTO 450
         ENDIF
      ELSE
         naxfl = 0
         IF ( ms.EQ.0 ) mend = MEMAX
         IF ( ms.EQ.0 ) ms = 1
         DO kh = ms , mend ! Loop over matrix elements
            DO ij = 1 , 2
               pv = (ELMU(kh)-ELML(kh))/100.D0
               IF ( ij.NE.1 .OR. (ELM(kh)-ELML(kh)).GE.pv ) THEN
                  IF ( ij.NE.2 .OR. (ELMU(kh)-ELM(kh)).GE.pv ) THEN
                     DO kh1 = 1 , MEMAX
                        SA(kh1) = 0.
                     ENDDO
                     IF ( IVAR(kh).EQ.0 ) GOTO 500
                     SA(kh) = 1.D0*(-1)**ij
                     kh1 = kh
                     CALL KONTUR(idr,chis0,chisl,ifbp,-1,kh1,sh,bten,
     &                           rem)
                     ELM(kh) = HLM(kh)
                  ENDIF
               ENDIF
            ENDDO
            REWIND 15
            WRITE (15,*) (DEVD(ij),DEVU(ij),ij=1,MEMAX)
 500        CONTINUE
         ENDDO
      ENDIF
 600  IF ( ifbp.EQ.1 ) THEN
         REWIND 17
         DO lkj = 1 , MEMAX
            READ (17,*) ELM(lkj)
         ENDDO
         WRITE (22,99029)
99029    FORMAT (1X///20X,'*** BEST POINT FOUND (TAPE17) ***'///)
         CALL PRELM(3)
      ENDIF
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      WRITE (22,99030)
99030 FORMAT (40X,'ESTIMATED ERRORS'//5X,'INDEX',5X,'NI',5X,'NF',5X,
     &        'ME AND ERRORS'//)
      DO kh1 = 1 , MEMAX
         IF ( IVAR(kh1).NE.0 .AND. IVAR(kh1).LE.999 ) THEN
            WRITE (22,99031) kh1 , LEAD(1,kh1) , LEAD(2,kh1) , HLM(kh1)
     &                       , DEVD(kh1) , DEVU(kh1) , DEVD(kh1)
     &                       *100./ABS(HLM(kh1)) , DEVU(kh1)
     &                       *100./ABS(HLM(kh1))
99031       FORMAT (6X,1I3,5X,1I3,4X,1I3,5X,1F9.5,2X,'(',1F9.5,' ,',
     &              1F9.5,')','......',1F7.1,' ,',1F7.1,1X,'PC')
         ENDIF
      ENDDO
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      WRITE (22,99032)
99032 FORMAT (40X,'ESTIMATED ERRORS',//5X,'INDEX',5X,'NI',5X,'NF',5X,
     &        'B(E,ML)(OR QUADRUPOLE MOMENT)',' AND ERRORS'//)

      DO kh2 = 1 , MEMAX
         IF ( IVAR(kh2).NE.0 .AND. IVAR(kh2).LE.999 ) THEN
            ispa = LEAD(2,kh2)
            IF ( LEAD(1,kh2).NE.LEAD(2,kh2) ) THEN
               sbe = 2.*SPIN(ispa) + 1.
               be2 = HLM(kh2)*HLM(kh2)/sbe
               be2a = HLM(kh2) + DEVD(kh2)
               be2b = HLM(kh2) + DEVU(kh2)
               be2c = be2b
               IF ( ABS(be2a).GT.ABS(be2b) ) be2b = be2a
               IF ( ABS(be2a-be2c).GT.1.D-6 ) be2a = be2c
               IF ( be2a/HLM(kh2).LE.0. .OR. be2b/HLM(kh2).LE.0. )
     &              be2a = 0.
               be2a = be2a**2/sbe
               be2b = be2b**2/sbe
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , be2 ,
     &                          be2a - be2 , be2b - be2
            ELSE
               ispb = INT(SPIN(ispa))*2
C 3.170662D0 = sqrt(16 pi / 5)
               qfac = SQRT(16.D0*pi/5.D0)*
     &                WTHREJ(ispb,4,ispb,-ispb,0,ispb)
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) ,
     &                          HLM(kh2)*qfac , DEVD(kh2)*qfac ,
     &                          DEVU(kh2)*qfac
            ENDIF
         ENDIF
      ENDDO
      GOTO 2000 ! Normal end of execution

 700  irea = 0
      IF ( ms.LT.0 ) irea = 1
      IF ( ms.EQ.0 ) mend = MEMAX
      IF ( ms.EQ.0 ) ms = 1
 800  naxfl = 1
      IF ( irea.EQ.1 ) READ (JZB,*) ms , mend
      IF ( ms.NE.0 ) THEN
         DO kh = ms , mend ! For matrix elements
            IF ( ifc.NE.1 ) THEN ! Use correlation matrix if IFC = 0
               REWIND 18
               DO kh1 = 1 , kh
                  READ (18,*) (KVAR(jyi),jyi=1,MEMAX)
               ENDDO
               DO kh1 = 1 , MEMAX ! For each matrix element
                  ivrh = IVAR(kh1)
                  IF ( KVAR(kh1).EQ.0 ) IVAR(kh1) = 0
                  KVAR(kh1) = ivrh
               ENDDO
            ENDIF
            DO ij = 1 , 2 ! Lower and upper
               sh = DEVU(kh)
               IF ( ij.EQ.1 ) sh = DEVD(kh)
               IF ( ABS(sh).LT.1.D-6 ) sh = (-1)**ij*ABS(HLM(kh))/10.D0
               ELM(kh) = HLM(kh) + 1.5D0*sh
               mm = 0
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IF ( ifc.EQ.1 ) KVAR(kh1) = IVAR(kh1)
                  mm = mm + IVAR(kh1)
               ENDDO
               IF ( mm.EQ.0 ) WRITE (22,99033) kh
99033          FORMAT (10X,'ME=',1I3,5X,'NO FREE MATRIX ELEMENTS')
               IF ( mm.NE.0 ) THEN
                  KFERR = 1
                  IF ( iosr.EQ.1 ) WRITE (IUNIT3,*) kh , kh ! For sigma program
                  IF ( iosr.EQ.1 ) WRITE (IUNIT3,*) kh , ij , ELM(kh)
                  LOCKS = 1
                  DLOCK = .05D0
                  CALL MINI(chiss,-1.D0,2,.0001D0,1000,idr,100000.D0,0,
     &                      iosr,kh,bten)
                  DO kh1 = 1 , MEMAX ! For each matrix element
                     SA(kh1) = (ELM(kh1)-HLM(kh1))/ABS(sh)
                  ENDDO
                  CALL KONTUR(idr,chis0,chisl,ifbp,inpo,kh,sh,bten,rem)
               ENDIF
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IF ( ifc.EQ.1 ) IVAR(kh1) = KVAR(kh1)
                  ELM(kh1) = HLM(kh1)
               ENDDO
            ENDDO ! Loop on lower and upper
            IF ( ifc.NE.1 ) THEN
               DO kh1 = 1 , MEMAX ! For each matrix element
                  IVAR(kh1) = KVAR(kh1)
               ENDDO
            ENDIF
            REWIND 15
            WRITE (15,*) (DEVD(kh1),DEVU(kh1),kh1=1,MEMAX)
         ENDDO ! Loop on matrix elements kh
         IF ( irea.EQ.1 ) GOTO 800
      ENDIF
      IF ( iosr.NE.0 ) THEN
         im = 0
         WRITE (IUNIT3,*) im , im
      ENDIF
      GOTO 600

C.............................................................................
 1500 WRITE (22,99043)
99043 FORMAT (5X,'ERROR-M.E. DOES NOT BELONG TO THE UPPER TRIANGLE')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1600 WRITE (22,99044)
99044 FORMAT (5X,'ERROR-WRONG SEQUENCE OF MULTIPOLARITIES')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1700 WRITE (22,99045)
99045 FORMAT (5X,'ERROR-REPEATED APPEARANCE OF THE STATE')
      GOTO 1900 ! Troubleshoot

C.............................................................................
 1800 WRITE (22,99046)
99046 FORMAT (1X///10X,'ERROR-INSUFFICIENT SPACE FOR E-THETA INTEGR ',
     &        'ATION')
      GOTO 1900 ! Troubleshoot

C.............................................................................
C     Troubleshooting
 1900 IF ( ITS.NE.0 ) THEN
         iva = 0
         WRITE (18,*) iva , iva , iva , chisq
         IF ( ITS.NE.2 ) THEN
            WRITE (15,*) iva , chisq , chisq , chisq , chisq
            CALL KLOPOT(kmat,rlr) ! Troubleshooting
         ENDIF
      ENDIF

C     End of execution
 2000 WRITE (22,99047)
99047 FORMAT (15X,'********* END OF EXECUTION **********')

C---- gosia2 addition start
C     Decide if we have to loop again for beam/projectile
 2100 IF ( mrepf.EQ.1 ) THEN
         mrepf = 2
         JZB = 26
         GOTO 2200
      ENDIF

      IF ( mret.EQ.1 ) THEN
         JZB = 25
         GOTO 2200
      ENDIF

C     Write normalization coefficients
      DO mmmm = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
         DO kkkk = 1 , LP1 ! LP1 is 50 (maximum number of experiments)
            WRITE (13,*) CNOR1(mmmm,kkkk)
         ENDDO
      ENDDO
C---- gosia2 addition end
      STOP

C--------------------------------
C OP,BRIC
70000 CALL BRICC
      GOTO 100 ! End of OP,BRIC - back to input loop

C--------------------------------
C OP,CORR
70001 CALL READY(idr,ntap,0)
      REWIND 3
      REWIND 15
      REWIND 4
      GOTO 70008 ! End of OP,CORR

C--------------------------------
C OP,COUL and OP,GOSI
70002 READ (JZB,99023) op1 ! Read the suboption
99023 FORMAT (1A4)
      IF ( op1.EQ.'    ' ) THEN
        IF (levmax .ne. NMAX) THEN
          WRITE(*,*) 'ERROR: Max. level index (',levmax,
     &      ') .NE. number of levels (',NMAX,')'
          STOP 'INPUT ERROR'
        ENDIF
        GOTO 100 ! Back to input loop
      ENDIF

C     Treat suboption LEVE (levels)
      IF ( op1.EQ.'LEVE' ) THEN
         NMAX = 0
         IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99024)
99024    FORMAT (1X/40X,'LEVELS',//5X,'INDEX',5X,'PARITY',9X,'SPIN',11X,
     &           'ENERGY(MEV)')
         ndima = NDIM + 1
         DO k = 1 , ndima
            READ (JZB,*) ipo1 , ipo2 , po2 , po1 ! level number, parity, spin, energy
            IF ( ipo1.EQ.0 ) GOTO 70002
            IF ( ipo1.EQ.1 .AND. ABS(po2).LT.1.D-6 ) ISO = 0
            NMAX = NMAX + 1
            SPIN(ipo1) = po2
            IF ( k.EQ.1 ) iph = ipo2
            iprc = ipo2 - iph
            IF ( iprc.NE.0 ) iprc = 1
            IFAC(ipo1) = (-1)**(iprc-INT(po2-SPIN(1)))
            EN(ipo1) = po1
            prp = '+'
            IF ( ipo2.EQ.-1 ) prp = '-'
            if (ipo1 .GT. levmax) levmax = ipo1
            IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99025) ipo1 , prp ,
     &           SPIN(ipo1) , EN(ipo1)
99025       FORMAT (5X,1I3,11X,1A1,10X,1F4.1,8X,1F10.4)
         ENDDO
         WRITE(*,*) 'Too many levels! - limit is ',NDIM
         STOP 'ERROR'

C     Treat suboption ME (matrix elements)
      ELSEIF ( op1.EQ.'ME  ' ) THEN
         DO k = 1 , nmemx
            IF ( op2.EQ.'GOSI' ) THEN
               READ (JZB,*) ipo1 , ipo2 , po1 , bl , bu ! lamda, 0, 0, 0, 0 OR ind1, ind2, me, lo, hi
               iopri = 2
               icg = 2
            ELSE
               iopri = 1
               READ (JZB,*) ipo1 , ipo2 , po1 ! lambda, 0, 0 OR ind1, ind2, me
            ENDIF
            IF ( ipo1.NE.0 ) THEN
               IF ( ipo2.EQ.0 ) THEN
                  IF ( ipo1.LE.la ) GOTO 1600 ! Error - wrong sequence of multipolarities
                  LAMMAX = LAMMAX + 1
                  LAMDA(LAMMAX) = ipo1
                  ipo3 = 0
                  IF ( indx.EQ.0 ) GOTO 220
               ELSE
                  MULTI(la) = MULTI(la) + 1
                  indx = indx + 1
                  IF ( ipo1.GT.NMAX ) THEN
                    WRITE(*,*) 'ERROR: Invalid initial state ',ipo1
                    STOP 'INPUT ERROR'
                  ENDIF
                  IF ( ABS(ipo2).GT.NMAX ) THEN
                    WRITE(*,*) 'ERROR: Invalid final state ',ipo2
                    STOP 'INPUT ERROR'
                  ENDIF
                  IF ( ipo1.GT.ABS(ipo2) ) GOTO 1500 ! Error - M.E. does not belong to the upper triangle
                  IF ( ipo1.NE.ipo3 ) THEN
                     IF ( ipo1.LT.ipo3 ) GOTO 1700 ! Error - repeated appearance of the state
                     ipo3 = ipo1
                  ENDIF
                  ELM(indx) = po1
                  IF ( JZB.EQ.25 ) ELM25(indx) = po1 ! Added for gosia2
                  IF ( JZB.EQ.26 ) ELM26(indx) = po1 ! Added for gosia2
                  MLT(indx) = la
                  LEAD(1,indx) = ipo1
                  LEAD(2,indx) = ABS(ipo2)
                  LDNUM(la,ipo1) = LDNUM(la,ipo1) + 1
                  IF ( op2.EQ.'GOSI' ) THEN
                     IF ( ipo2.LT.0 ) THEN ! If negative, bl and bu are indices
                                           ! to which we fix this element
                        IVAR(indx) = 10000*INT(bl) + INT(bu)
                     ELSE                 ! Otherwise they are limits
                        ELMU(indx) = bu
                        ELML(indx) = bl
                        IF ( ABS(bl-bu).LT.1.D-6 ) THEN
                           IVAR(indx) = 0 ! Fixed
                        ELSE
                           IVAR(indx) = 2
                           IF ( la.GT.4 ) IVAR(indx) = 1
                        ENDIF
                     ENDIF
                     isip = ISEX(ipo1) + 1
                     ISEX(ABS(ipo2)) = MIN(isip,ISEX(ABS(ipo2)))
                  ENDIF
                  GOTO 250
               ENDIF
            ENDIF
            DO kk = 1 , indx
               IF ( ABS(ELM(kk)).LE.1.D-6 ) ELM(kk) = 1.D-6
               IF ( IVAR(kk).GE.10000 ) THEN ! Correlated
                  kk1 = IVAR(kk)/10000
                  kk2 = IVAR(kk) - 10000*kk1
                  la1 = la
                  IF ( kk2.GE.100 ) THEN
                     la1 = kk2/100
                     kk2 = kk2 - 100*la1
                  ENDIF
                  inx1 = MEM(kk1,kk2,la1)
                  ELML(kk) = ELML(inx1)*ELM(kk)/ELM(inx1)
                  ELMU(kk) = ELMU(inx1)*ELM(kk)/ELM(inx1)
                  SA(kk) = ELM(kk)/ELM(inx1)
                  ivari(kk) = IVAR(kk)
                  IVAR(kk) = 1000 + inx1
                  IF ( ELMU(kk).LE.ELML(kk) ) THEN
                     elmi = ELMU(kk)
                     ELMU(kk) = ELML(kk)
                     ELML(kk) = elmi
                  ENDIF
               ENDIF
            ENDDO
            IF ( ipo1.EQ.0 ) GOTO 300
 220        la = ipo1
            IF ( la.GT.LMAXE .AND. la.LE.6 ) LMAXE = la
 250        CONTINUE
         ENDDO
 300     MEMAX = indx
         IF ( la.GT.6 ) MAGEXC = 1 ! Flag that we need magnetic excitations
         memx4 = MULTI(1) + MULTI(2) + MULTI(3) + MULTI(4)
         MEMX6 = memx4 + MULTI(5) + MULTI(6)
         IF ( ABS(IPRM(1)).EQ.1 ) CALL PRELM(iopri)
         DO kh = 1 , NMAX
            IF ( ISEX(kh).EQ.1111 ) ISEX(kh) = 1
         ENDDO
         DO kh = 1 , MEMAX
            ivarh(kh) = IVAR(kh)
         ENDDO

C     Treat suboption CONT (control)
      ELSEIF ( op1.EQ.'CONT' ) THEN
 350     READ (JZB,'(A80)') line
         READ (line,99026,ERR=353) op1 , fipo1
 353     READ (line,99056,ERR=354) op1 , ipo1
         fipo1 = ipo1
99026    FORMAT (1A4,1F7.1)
99056    FORMAT (1A4,I5)
 354     ipo1 = INT(fipo1)
         IF ( op1.EQ.'ACC,' ) THEN
           ACCUR = 10.D0**(-fipo1)
         ELSE IF ( op1.EQ.'ACP,' ) THEN
           ACCA = 10.D0**(-fipo1)
         ELSE IF ( op1.EQ.'CCF,' ) THEN
           IPS1 = ipo1
         ELSE IF ( op1.EQ.'CRF,' ) THEN
           ICS = 1
         ELSE IF ( op1.EQ.'CRD,' ) THEN
           DO jjx = 1 , ipo1
             READ (JZB,*) ipo2
             iecd(ipo2) = 1
           ENDDO
         ELSE IF ( op1.EQ.'DIP,' ) THEN
           DIPOL = 0.001D0*fipo1
         ELSE IF ( op1.EQ.'EFF,' ) THEN
           DO jjx = 1 , ipo1
             READ (JZB,*) ipo2 , ijx
             ideff(ipo2) = ijx
           ENDDO
         ELSE IF ( op1.EQ.'END,' ) THEN
           GOTO 70002
         ELSE IF ( op1.EQ.'FIX,' ) THEN
           READ (JZB,*) nallow
           DO jjx = 1 , nallow
             READ (JZB,*) ijk
             IVAR(ijk) = -IVAR(ijk)
           ENDDO
           DO jjx = 1 , MEMAX
             IF ( IVAR(jjx).GE.0 ) THEN
               IF ( IVAR(jjx).LE.999 ) IVAR(jjx) = 0
             ENDIF
           ENDDO
           DO jjx = 1 , MEMAX
             IF ( IVAR(jjx).LT.0 ) IVAR(jjx) = -IVAR(jjx)
             ivarh(jjx) = IVAR(jjx)
           ENDDO
         ELSE IF ( op1.EQ.'FMI,' ) THEN
           ifm = 1
         ELSE IF ( op1.EQ.'INR,' ) THEN
           INNR = 1
         ELSE IF ( op1.EQ.'INT,' ) THEN
           DO jjx = 1 , ipo1
             READ (JZB,*) ipo2 , ijx
             INTERV(ipo2) = ijx
           ENDDO
         ELSE IF ( op1.EQ.'LCK,' ) THEN
 352       READ (JZB,*) lck1 , lck2
           IF ( lck1.EQ.0 ) GOTO 350
           DO jjx = lck1 , lck2
             ivarh(jjx) = 0
             IVAR(jjx) = 0
           ENDDO
           GOTO 352
         ELSE IF ( op1.EQ.'NCM,' ) THEN
           NCM = ipo1
         ELSE IF ( op1.EQ.'PIN,' ) THEN
           ipine = ipo1
           ipinf = 1
           DO ipp = 1 , ipine
             READ (JZB,*) ig1 , ig2
             jpin(ig1) = ig2
           ENDDO
         ELSE IF ( op1.EQ.'PRT,' ) THEN
           DO jjx = 1 , 20
             READ (JZB,*) inm1 , inm2
             IF ( inm1.EQ.0 ) GOTO 350
             IPRM(inm1) = inm2
           ENDDO
         ELSE IF ( op1.EQ.'SEL,' ) THEN
           ITS = 2
         ELSE IF ( op1.EQ.'SKP,' ) THEN
           DO jjx = 1 , ipo1
             READ (JZB,*) ijx
             JSKIP(ijx) = 0
           ENDDO
         ELSE IF ( op1.EQ.'SMR,' ) THEN
           iosr = 1
         ELSE IF ( op1.EQ.'SPL,' ) THEN
           ISPL = ipo1
         ELSE IF ( op1.EQ.'TEN,' ) THEN
           itno = 1
         ELSE IF ( op1.EQ.'VAC,' ) THEN
           DO jjx = 1 , 7
             READ (JZB,*) ijx , val
             IF ( ijx.EQ.0 ) GOTO 350
             G(ijx) = val
           ENDDO
         ELSE IF ( op1.EQ.'WRN,' ) THEN
           SGW = fipo1
         ELSE
           WRITE(*,*) 'Did not understand CONT option ',op1
         ENDIF
         GOTO 350 ! Back to beginning of CONT loop

C     Treat suboption EXPT
      ELSEIF ( op1.EQ.'EXPT' ) THEN
         READ (JZB,*) NEXPT , IZ , XA
         G(1) = 3.D0           ! AVJI
         G(2) = .02D0          ! GAMMA
         G(3) = .0345D0        ! XLAMB
         G(4) = 3.5D0          ! TIMEC
         G(5) = DBLE(IZ)/XA    ! GFAC
         G(6) = 6.D-06         ! FIEL
         G(7) = .6D0           ! POWER
         DO k = 1 , NEXPT ! Zn, An, E_p, THETA_lab, M_c, M_A, IAX, phi1, phi2, ikin, ln
            READ (JZB,*) IZ1(k) , XA1(k) , EP(k) , TLBDG(k) , EMMA(k) ,
     &           MAGA(k) , IAXS(k) , fi0 , fi1 , ISKIN(k) , LNORM(k)
            ITTE(k) = 0
            IF ( XA1(k).LT.0. ) ITTE(k) = 1
            XA1(k) = ABS(XA1(k))
            FIEX(k,1) = fi0*pi/180.D0 ! Convert to radians
            FIEX(k,2) = fi1*pi/180.D0
            IF ( TLBDG(k).LT.0. ) THEN
               FIEX(k,1) = FIEX(k,1) + pi
               FIEX(k,2) = FIEX(k,2) + pi
            ENDIF
         ENDDO

C     Else we don't recognize the suboption
      ELSE
         WRITE (22,99027) op1
99027    FORMAT (5X,'UNRECOGNIZED SUBOPTION',1X,1A4)
         GOTO 2000 ! Normal end of execution
      ENDIF
      GOTO 70002 ! Get next suboption

C--------------------------------
C OP,ERRO
70003 READ (JZB,*) idf , ms , mend , irep , ifc , remax
C---- gosia2 changes start
      MCFIX = 0
      DO mmmm = 1 , 32
        DO kkkk = 1 , 50
          READ (13,*) CNOR1(mmmm,kkkk)
        ENDDO
      ENDDO
C---- gosia2 changes end
      rem = LOG(remax)
      LOCKS = 0
      LOCKF = 0
      JENTR = 1 ! Flag to indicate we are in OP,ERRO
      sh = 1.
      ifbp = 0
      inpo = 1
      inko = 1
      IF ( iosr.NE.0 .AND. idf.NE.0 ) THEN
        inn = 0
        ij = MULTI(1)
        IF ( ij.NE.0 ) THEN
          DO ij = 1 , NMAX
            lxd = LDNUM(1,ij)
            IF ( lxd.NE.0 ) THEN
              DO ijk = 1 , lxd
                inn = inn + 1
              ENDDO
            ENDIF
          ENDDO
          inpo = inn + 1
        ENDIF
        DO ij = 1 , NMAX
          lxd = LDNUM(2,ij)
          IF ( lxd.NE.0 ) THEN
            DO ijk = 1 , lxd
              inn = inn + 1
            ENDDO
          ENDIF
        ENDDO
        inko = inn
        IF ( irep.NE.2 ) THEN
          WRITE (IUNIT3,*) NMAX , MEMAX , inpo , inko
          DO inn = 1 , NMAX
            WRITE (IUNIT3,*) inn , SPIN(inn) , EN(inn)
          ENDDO
          DO inn = 1 , MEMAX
            WRITE (IUNIT3,*) inn , LEAD(1,inn) , LEAD(2,inn)
          ENDDO
          DO inn = 1 , MEMAX
            WRITE (IUNIT3,*) inn , ELM(inn)
          ENDDO
        ENDIF ! IF ( irep.NE.2 )
      ENDIF ! IF ( iosr.NE.0 .AND. idf.NE.0 )
      IF ( irep.NE.0 ) THEN
        REWIND 15
        READ (15,*) (DEVD(kh1),DEVU(kh1),kh1=1,MEMAX)
      ELSE
        DO kh1 = 1 , MEMAX
          DEVD(kh1) = ELML(kh1) - ELM(kh1)
          DEVU(kh1) = ELMU(kh1) - ELM(kh1)
        ENDDO
      ENDIF
      IF ( IMIN.EQ.0 ) CALL CMLAB(0,dsig,ttttt)
      IF ( ERR ) GOTO 2000 ! Error
      IF ( IMIN.NE.0 ) GOTO 400
      GOTO 1300 ! End of OP,ERRO

C--------------------------------
C OP,EXIT
C---- gosia2 changes start
70004 IF ( mret.EQ.1 .AND. JZB.EQ.26 ) nawr = 0
      IF ( mret.EQ.1 ) THEN
        GOTO 2100
      ENDIF
      IF ( IPRM(18).NE.0 ) CALL PTICC(idr)
C---- gosia2 changes end
      IF ( oph.EQ.'GOSI' ) THEN
         IF ( lfagg.NE.1 ) THEN
            IF ( IMIN.NE.0 ) THEN
               IF ( IPRM(4).EQ.-1 ) IPRM(4) = 111111
               iskok = IPRM(7) + IPRM(8) + IPRM(13) + IPRM(14)
               IF ( iskok.NE.0 .OR. IPRM(4).NE.111111 ) THEN
                  IF ( iskok.NE.0 ) THEN
                     IF ( IPRM(7).EQ.1 ) IPRM(7) = -1
                     IF ( IPRM(8).EQ.1 ) IPRM(8) = -1
                     IF ( IPRM(3).EQ.1 .AND. NBRA.NE.0 ) IPRM(3) = -1
                     IF ( IPRM(13).EQ.1 ) IPRM(13) = -1
                     IF ( IPRM(14).EQ.1 ) IPRM(14) = -1
                  ENDIF
                  CALL MINI(chisq,chiok,+1,conu,2000,idr,xtest,2,0,0,
     &                      bten)
               ENDIF
            ENDIF
            CALL MIXR(iva,1,chisq,chilo)
            IF ( IPRM(15).NE.0 .AND. KFERR.NE.1 .AND. iyr.NE.0 ) THEN
               WRITE (22,99011)
99011          FORMAT (1X//20X,'CALCULATED LIFETIMES'//5X,'LEVEL',5X,
     &                 'LIFETIME(PSEC)',5X,'EXP',8X,'ERROR'/)
               DO iva = 2 , NMAX
                  DO iva1 = 1 , NLIFT
                     IF ( LIFCT(iva1).EQ.iva ) GOTO 122
                  ENDDO
                  WRITE (22,99012) iva , TAU(iva)
99012             FORMAT (6X,1I3,7X,1E10.4)
                  GOTO 124
 122              WRITE (22,99013) iva , TAU(iva) , TIMEL(1,iva1) ,
     &                   -TIMEL(2,iva1) , TIMEL(3,iva1)
99013             FORMAT (6X,1I3,7X,1E10.4,5X,1E10.4,4X,1E10.4,4X,
     &              1E10.4)
 124              IF ( iva.EQ.NMAX ) THEN
                    IF ( NAMX.GE.1 ) THEN
                      WRITE (22,99014)
99014                 FORMAT (5x,//,
     &         'CALCULATED AND EXPERIMENTAL MATRIX ELEMENTS'
     &                    ,//)
                        WRITE (22,99015)
99015                   FORMAT (5x,'NI ','NF ',' EXP. ME   ',
     &                     'CURRENT ME','   SIGMA')
                        DO kq = 1 , NAMX
                           ni = IAMY(kq,1)
                           nf = IAMY(kq,2)
                           ind = IAMX(kq)
                           ess = ELM(ind)
                           esd = EAMX(kq,1)
                           IF ( ess.LT.esd) THEN
                             dsd = EAMX(kq,2)
                           ELSE
                             dsd = EAMX(kq,3)
                           ENDIF
                           WRITE (22,99016) ni , nf , esd , ess ,
     &                        (ess-esd)/dsd
99016                      FORMAT (4x,1I3,1x,1I3,1x,1F9.4,1x,1F9.4,1x,
     &                              1F9.4)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF ( IMIN.NE.0 ) CALL PRELM(3)
         ENDIF
      ENDIF
      GOTO 1900 ! End of OP,EXIT - troubleshoot

C--------------------------------
C OP,GDET
70005 nl = 7
      READ (JZB,*) nfdd ! number of physical detectors

      nfd = ABS(nfdd) ! Negative value means graded absorber
      IF ( nfdd.LE.0 ) THEN
        REWIND 8
        DO i = 1 , nl
          WRITE (8,*) (tau2(l,i),l=1,10)
        ENDDO
        WRITE (8,*) (eng(l),l=1,10)
      ENDIF

C     Write file for gamma-ray energy dependence of Ge solid-angle
C     attenuation coefficients
      REWIND 9
      WRITE (9,*) nfd
      DO i = 1 , nfd ! For each detector
        READ (JZB,*) (DIX(k),k=1,4) ! radius of core, outer radius, length, distance
        READ (JZB,*) (xl1(k),k=1,nl) ! thicknesses of 7 kinds of absorber
        IF ( DIX(1).LE.0. ) DIX(1) = .01D0
        WRITE (9,*) DIX(4) ! length
        IF ( nfdd.LE.0 ) WRITE (8,*) (xl1(k),k=1,nl)
        ind = 1
        IF ( xl1(5).GT.0. ) ind = 3
        IF ( xl1(6).GT.0. ) ind = 4
        IF ( xl1(7).GT.0. ) ind = 5
        WRITE (9,*) eng(ind) ! First energy
        CALL QFIT(qui,tau1,tau2,eng,xl1,cf,nl,ind)
        WRITE (22,99004) i
99004   FORMAT (10X,'DETECTOR',1X,1I2)
        DO k = 1 , 8
          WRITE (22,99005) k , cf(k,1) , cf(k,2)
99005     FORMAT (1X,//5X,'K=',1I1,2X,'C1=',1E14.6,2X,'C2=',
     &      1E14.6/5X,'ENERGY(MEV)',5X,'FITTED QK',5X,
     &      'CALC.QK',5X,'PC.DIFF.'/)
          WRITE (9,*) cf(k,1) , cf(k,2) , qui(k,ind)
          DO l = 1 , 10
            arg = (eng(l)-eng(ind))**2
            qc = (qui(k,ind)*cf(k,2)+cf(k,1)*arg)/(cf(k,2)+arg)
            WRITE (22,99006) eng(l) , qc , qui(k,l) ,
     &        100.D0*(qc-qui(k,l))/qui(k,l)
99006       FORMAT (8X,1F4.2,6X,1F9.4,5X,1F9.4,3X,1E10.2)
          ENDDO
        ENDDO
      ENDDO
      GOTO 100 ! End of OP,GDET - back to input loop

C--------------------------------
C OP,INTG
70006 CALL OP_INTGI(iretval, op2, ipinf, jpin, iecd, izcap, lfagg,
     &  idr, bten, fiex1, 0)
      IF ( iretval.eq.1800 ) GOTO 1800
      IF ( iretval.eq.1900 ) GOTO 1900
      IF ( iretval.eq.2000 ) GOTO 2000
      GOTO 100 ! End of OP,INTG - back to input loop

C--------------------------------
C OP,INTI
70007 CALL OP_INTGI(iretval, op2, ipinf, jpin, iecd, izcap, lfagg,
     &  idr, bten, fiex1, 1)
      IF ( iretval.eq.1800 ) GOTO 1800
      IF ( iretval.eq.1900 ) GOTO 1900
      IF ( iretval.eq.2000 ) GOTO 2000
      GOTO 100 ! End of OP,INTI - back to input loop

C--------------------------------
C OP,POIN
70008 CALL CMLAB(0,dsig,ttttt) ! Options MAP, STAR, POINT, MINI etc.
      IF ( ERR ) GOTO 2000 ! Error
      IF ( op2.EQ.'POIN' ) READ (JZB,*) ifwd , slim
      ient = 1
      icg = 1
      IF ( SPIN(1).LT.1.D-6 ) ISO = 0
      IF ( iobl.LT.1 ) THEN
         IF ( op2.NE.'GOSI' ) THEN
            iapx = 0
            DO ii = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
               ILE(ii) = 1
            ENDDO
            nch = 0
            DO jexp = 1 , NEXPT ! For each experiment
               IEXP = jexp
               ttttt = TREP(IEXP)
               dsig = DSIGS(IEXP)
               IF ( op2.NE.'STAR' ) THEN
                  jmm = IEXP
                  IF ( IEXP.NE.1 ) THEN
                     DO lli = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
                        ILE(lli) = ILE(lli) + NYLDE(IEXP-1,lli)
                     ENDDO
                  ENDIF
               ENDIF
               fi0 = FIEX(IEXP,1) ! Lower phi limit
               fi1 = FIEX(IEXP,2) ! Upper phi limit
               CALL LOAD(IEXP,1,icg,0.D0,jj)
               CALL ALLOC(ACCUR)
               CALL SNAKE(IEXP,ZPOL)
               CALL SETIN
               DO j = 1 , LMAX ! For each spin up to ground-state spin + 1
                  polm = DBLE(j-1) - SPIN(1)
                  CALL LOAD(IEXP,2,icg,polm,jj)
                  CALL STING(jj)
                  CALL PATH(jj)
                  CALL INTG(IEXP)
                  CALL TENB(j,bten,LMAX)
                  pr = 0.
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99034) (DBLE(j)-1.-SPIN(1)) , IEXP
99034             FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',
     &                    1F5.1,5X,'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,
     &                    'SPIN',2X,'M',5X,'REAL AMPLITUDE',2X,
     &                    'IMAGINARY AMPLITUDE'//)
                  DO k = 1 , ISMAX ! For substates
                     pr = pr + DBLE(ARM(k,5))**2 + DIMAG(ARM(k,5))**2
                     IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                    WRITE (22,99035) INT(CAT(k,1)) , CAT(k,2) ,
     &                    CAT(k,3) , DBLE(ARM(k,5)) , DIMAG(ARM(k,5))
99035                FORMAT (7X,1I2,3X,1F4.1,2X,1F5.1,2X,1E14.6,2X,
     &                       1E14.6)
                  ENDDO ! Loop on substates k
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99036) pr
99036             FORMAT (1X/5X,'SUM OF PROBABILITIES=',1E14.6)
               ENDDO ! Loop over spins j
               CALL TENS(bten)
               IF ( itno.NE.0 ) THEN ! write statistical tensors on tape 17
                  DO k = 2 , NMAX
                     WRITE (17,*) k
                     DO kk = 1 , 4
                        in1 = (k-1)*28 + 1 + (kk-1)*7
                        in2 = in1 + 2*kk - 2
                        WRITE (17,*) (ZETA(kkk),kkk=in1,in2)
                     ENDDO
                  ENDDO
               ENDIF
               summm = 0.
               DO jgl = 2 , NMAX
                  loct = (jgl-1)*28 + 1
                  summm = summm + ZETA(loct)
               ENDDO
               pop1 = 1. - summm
               jgl = 1
               IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 ) WRITE (22,99053)
     &              jgl , pop1
               DO jgl = 2 , NMAX
                  loct = (jgl-1)*28 + 1
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99053) jgl , ZETA(loct)
               ENDDO
               IF ( op2.NE.'STAR' ) THEN
                  CALL DECAY(ccd,0,ccc)
                  nogeli = NANG(IEXP) ! Number of detector angles for expt
                  jgl1 = 0
                  DO js = 1 , LP2 ! LP2 = 1500 (maximum number of matrix elements)
                     DO jgl = 1 , 20
                        SUMCL(jgl,js) = 0.
                     ENDDO
                  ENDDO
                  DO jgl = 1 , nogeli ! For each detector angle
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                       WRITE (23,99037) IEXP , jgl , EP(IEXP) ,
     &                       ABS(TLBDG(IEXP))
99037                   FORMAT (1x//50x,'CALCULATED YIELDS'//5x,
     &                          'EXPERIMENT ',1I2,2x,'DETECTOR ',1I2/5x,
     &                          'ENERGY ',1F10.3,1x,'MEV',2x,'THETA ',
     &                          1F7.3,1x,'DEG'//5x,'NI',5x,'NF',5x,'II',
     &                          5x,'IF',5x,'E(MeV)',5x,'EFFICIENCY'/)
                     ENDIF
                     gth = AGELI(IEXP,jgl,1)
                     figl = AGELI(IEXP,jgl,2)
                     fm = (fi0+fi1)/2.
                     CALL ANGULA(YGN,idr,1,fi0,fi1,ttttt,gth,figl,jgl,
     &                 op2)
                     IF ( IFMO.NE.0 ) THEN
                        id = ITMA(IEXP,jgl) ! Get identity of detector
                        d = ODL(id) ! Get results of OP,GDET for that detector
                        rx = d*SIN(gth)*COS(figl-fm) -
     &                       .25D0*SIN(ttttt)*COS(fm)
                        ry = d*SIN(gth)*SIN(figl-fm) -
     &                       .25D0*SIN(ttttt)*SIN(fm)
                        rz = d*COS(gth) - .25D0*COS(ttttt)
                        rl = SQRT(rx*rx+ry*ry+rz*rz)
                        thc = TACOS(rz/rl)
                        sf = d*d/rl/rl
                        fic = ATAN2(ry,rx)
                        CALL ANGULA(YGP,idr,1,fi0,fi1,ttttt,thc,fic,jgl,
     &                    op2)
                        DO ixl = 1 , idr
                           ixm = KSEQ(ixl,3)
                           tfac = TAU(ixm)
C c/s = 0.11991698 /ps, where s = 0.25 cm, c = speed of light in cm/ps
                           IF ( tfac*BETAR(IEXP).GT. 25.D0 ) THEN
                             WRITE(22,99057) IEXP,KSEQ(ixl,3),tfac
                             IFMO = 0
                           ELSE
                             YGN(ixl) = YGN(ixl)
     &                         + .11991698D0*tfac*BETAR(IEXP)
     &                         *(sf*YGP(ixl)-YGN(ixl))
                           ENDIF
                        ENDDO
                     ENDIF
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        ipd = ITMA(IEXP,jgl) ! Get identity of detector
                        DO jyi = 1 , idr ! For each decay
                           ni = KSEQ(jyi,3)
                           nf = KSEQ(jyi,4)
                           decen = EN(ni) - EN(nf)
                           cocos = SIN(ttttt)*SIN(gth)*COS(fm-figl)
     &                             + COS(ttttt)*COS(gth)
                           decen = decen*(1.D0+BETAR(IEXP)*cocos)
                           CALL EFFIX(IEXP,ipd,decen,effi)
                           IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                          WRITE (23,99049) ni , nf , SPIN(ni) ,
     &                                 SPIN(nf) , decen , effi
                           YGN(jyi) = YGN(jyi)*effi
                        ENDDO
                        inclus = ICLUST(IEXP,jgl) ! Cluster number for detector jgl
                        IF ( inclus.NE.0 ) THEN
                           DO jyi = 1 , idr ! For each decay
                              SUMCL(inclus,jyi) = SUMCL(inclus,jyi)
     &                           + YGN(jyi)
                           ENDDO
                           IF ( jgl.NE.LASTCL(IEXP,inclus) ) GOTO 1205 ! If it is not the last detector in the cluster
                           DO jyi = 1 , idr ! For each decay
                              YGN(jyi) = SUMCL(inclus,jyi)
                           ENDDO
                        ENDIF
                     ENDIF
                     jgl1 = jgl1 + 1
                     lu = ILE(jgl1)
                     IF ( op2.EQ.'POIN' .OR. IPRM(11).EQ.1 )
     &                    WRITE (22,99048) IEXP , jgl1 , EP(IEXP) ,
     &                    TLBDG(IEXP)
                     jmm = 0
                     DO jyi = 1 , idr
                        ni = KSEQ(jyi,3)
                        nf = KSEQ(jyi,4)
                        IF ( op2.EQ.'POIN' .OR. IPRM(11).EQ.1 )
     &                       WRITE (22,99049) ni , nf , SPIN(ni) ,
     &                       SPIN(nf) , YGN(jyi) , YGN(jyi)/YGN(IDRN)
                        IF ( ifwd.EQ.1 ) THEN
                           IF ( (YGN(jyi)/YGN(IDRN)).GE.slim ) THEN
                              IF ( jgl1.EQ.1 ) sh1 = YGN(IDRN)
                              jmm = jmm + 1
                              CORF(jmm,1) = DBLE(ni)
                              CORF(jmm,2) = DBLE(nf)
                              CORF(jmm,3) = YGN(jyi) ! Was divided by sh1 in gosia
                              IF ( YGN(jyi).GE.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)/20.D0
                              IF ( YGN(jyi).LT.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)
     &                               *(.05D0+.2D0*(1.D0-YGN(jyi)/
     &                               YGN(IDRN)))
                           ENDIF
                        ENDIF
                        IF ( op2.EQ.'CORR' ) THEN
                           READ (15,*) yydd
                           nch = nch + 1
                           jjjj = IY(lu,jgl1)/1000
                           jyi1 = IY(lu,jgl1) - jjjj*1000
                           IF ( IY(lu,jgl1).EQ.jyi .OR. jjjj.EQ.jyi .OR.
     &                          jyi1.EQ.jyi ) THEN
                              IF ( IY(lu,jgl1).GE.1000 ) THEN
                                 jyi2 = jyi1 - jjjj
                                 IF ( jyi2.LE.0 ) GOTO 1202
                                 DO ihuj = 1 , jyi2
                                    READ (15,*) yyd1
                                 ENDDO
                                 yydd = yydd + yyd1
                                 YGN(jyi) = YGN(jyi) + YGN(jyi1)
                                 REWIND 15
                                 DO ihuj = 1 , nch
                                    READ (15,*) yyd1
                                 ENDDO
                              ENDIF
                              IF ( IEXP.EQ.1 .AND. lu.EQ.NYLDE(1,1)
     &                             .AND. jgl1.EQ.1 )
     &                             cnst = yydd/YGN(jyi)
                              CORF(lu,jgl1) = YEXP(jgl1,lu)
                              YEXP(jgl1,lu) = YEXP(jgl1,lu)
     &                           /yydd*YGN(jyi)
                              DYEX(jgl1,lu) = DYEX(jgl1,lu)
     &                           /yydd*YGN(jyi)
                              lu = lu + 1
                           ENDIF
                        ENDIF
 1202                   CONTINUE
                     ENDDO
                     IF ( ifwd.EQ.1 ) THEN
                        xw = 1.
                        WRITE (4,*) IEXP , jgl1 , ABS(IZ1(IEXP)) ,
     &                              ABS(XA1(IEXP)) , ABS(EP(IEXP)) ,
     &                              jmm , xw
                        DO jyi = 1 , jmm
                           WRITE (4,*) INT(CORF(jyi,1)) ,
     &                                 INT(CORF(jyi,2)) , CORF(jyi,3) ,
     &                                 CORF(jyi,4)
                        ENDDO
                     ENDIF
 1205                CONTINUE
                  ENDDO ! Loop on detector angles jgl
                  IF ( op2.EQ.'CORR' ) THEN
                     jgl1 = 0
                     DO jgl = 1 , nogeli ! For each detector
                        IF ( IRAWEX(jexp).NE.0 ) THEN
                           inclus = ICLUST(jexp,jgl) ! Cluster number for detector jgl
                           IF ( inclus.NE.0 ) THEN
                              IF ( jgl.NE.LASTCL(jexp,inclus) ) ! If detector is not the last in the cluster
     &                             GOTO 1206
                           ENDIF
                        ENDIF
                        jgl1 = jgl1 + 1
                        READ (3,*) ne , na , zp , ap , xep , nval , waga
                        WRITE (4,*) ne , na , zp , ap , EP(IEXP) ,
     &                              nval , waga
                        WRITE (22,99038) IEXP , jgl1
99038                   FORMAT (///10X,'EXPERIMENT',1X,I2,8X,'DETECTOR',
     &                          1X,I2,//9X,'NI',5X,'NF',5X,'YEXP',8X,
     &                          'YCOR',8X,'COR.F'/)
                        ile1 = ILE(jgl1)
                        DO itp = 1 , nval
                           READ (3,*) ns1 , ns2 , fiex1(1,1,1) ,
     &                                fiex1(1,1,2)
                           ltrn = IY(ile1+itp-1,jgl1)
                           IF ( ltrn.LT.1000 ) THEN
                              ns1 = KSEQ(ltrn,3)
                              ns2 = KSEQ(ltrn,4)
                           ELSE
                              ltrn1 = ltrn/1000
                              ns1 = KSEQ(ltrn1,3)*100
                              ns2 = KSEQ(ltrn1,4)*100
                              ltrn2 = ltrn - ltrn1*1000
                              ns1 = ns1 + KSEQ(ltrn2,3)
                              ns2 = ns2 + KSEQ(ltrn2,4)
                           ENDIF
                           ycorr = YEXP(jgl1,ile1+itp-1) ! Changed in gosia2
                           WRITE (4,*) ns1 , ns2 , ycorr ,
     &                                 DYEX(jgl1,ile1+itp-1) ! Changed in gosia2
                           WRITE (22,99039) ns1 , ns2 ,
     &                            CORF(ile1+itp-1,jgl1) , ycorr ,
     &                            ycorr/CORF(ile1+itp-1,jgl1)
99039                      FORMAT (5X,I4,5X,I4,3X,E8.3,4X,E8.3,4X,E8.3)
                        ENDDO ! Loop over itp
 1206                   CONTINUE
                     ENDDO ! Loop over jgl
                  ENDIF ! if ( op2.EQ. 'CORR')
               ENDIF
            ENDDO ! Loop over jexp
            IF ( op2.EQ.'STAR' ) oph = op2
            IF ( op2.NE.'STAR' ) THEN
               IF ( op2.EQ.'CORR' ) THEN
                  ntap = 4
                  CALL READY(idr,ntap,ipri)
                  REWIND ntap
               ENDIF
            ENDIF
            GOTO 100 ! Back to input loop
         ENDIF ! if (op2 .NE. 'GOSI') if statement
      ENDIF ! if ( iobl.LT.1 ) if statement

 1300 IF ( iobl.GE.1 ) THEN ! OP,ERRO
         ient = 1
         icg = 2
         nmaxh = NMAX
         lmax1 = LMAX
         sh1 = SPIN(1) ! Save ground-state spin
         sh2 = SPIN(2) ! Save spin of first excited state
         ih1 = IFAC(1)
         ih2 = IFAC(2)
         magh = MAGEXC
         lmaxh = LMAXE
         isoh = ISO
         ISO = 0
         eh1 = ELM(1)
         lh1 = LEAD(1,1)
         lh2 = LEAD(2,1)
         lamh = LAMMAX
         memh = MEMAX
         DO kh = 1 , 8 ! For each multipolarity
            ihlm(kh) = MULTI(kh)
            ihlm(kh+24) = LDNUM(kh,2)
            ihlm(kh+8) = LAMDA(kh)
            ihlm(kh+16) = LDNUM(kh,1)
         ENDDO
         DO jexp = 1 , NEXPT ! For each experiment
            IEXP = jexp
            intvh = INTERV(IEXP)
            DO jgs = 1 , MEMAX
               DO jgr = 1 , 7
                  QAPR(jgs,1,jgr) = 0.
               ENDDO
            ENDDO
            DO iuy = 1 , 6
               XIR(iuy,IEXP) = 0.
            ENDDO
            emhl1 = EMMA(IEXP)
            EMMA(IEXP) = DBLE(MAGA(IEXP))
            jde = 2
            IF ( MAGA(IEXP).EQ.0 ) jde = 1
            DO iuy = 1 , 6
               zmir(iuy,1,IEXP) = 0.
               zmir(iuy,2,IEXP) = 0.
            ENDDO
            CALL LOAD(IEXP,1,2,0.D0,jj)
            DO jgs = 1 , LMAX ! For each spin up to ground-state spin + 1
               polm = DBLE(jgs-1) - SPIN(1)
               CALL LOAD(IEXP,3,2,polm,jj)
               CALL PATH(jj)
               CALL LOAD(IEXP,2,2,polm,jj)
               ictl = 1
               DO kk = 1 , 6
                  ll = ihlm(kk)
                  IF ( ll.NE.0 ) THEN
                     lfini = ll + ictl - 1
                     ict = ictl
                     DO lll = ict , lfini
                        ictl = ictl + 1
                        IF ( jgs.EQ.1 ) XIR(kk,IEXP)
     &                       = MAX(XIR(kk,IEXP),ABS(XI(lll)))
                        r1 = ABS(QAPR(lll,1,1))
                        r2 = ABS(QAPR(lll,1,4))
                        r3 = ABS(QAPR(lll,1,7))
                        rm = MAX(r1,r2,r3)
                        bmx = MAX(ABS(ELMU(lll)),ABS(ELML(lll)))
                        zmir(kk,2,IEXP)
     &                     = MAX(zmir(kk,2,IEXP),rm*bmx/ABS(ELM(lll)),
     &                     rm)
                        r1 = ABS(QAPR(lll,1,2))
                        r2 = ABS(QAPR(lll,1,3))
                        r3 = ABS(QAPR(lll,1,5))
                        r4 = ABS(QAPR(lll,1,6))
                        rm = MAX(r1,r2,r3,r4)
                        zmir(kk,1,IEXP)
     &                     = MAX(zmir(kk,1,IEXP),rm*bmx/ABS(ELM(lll)),
     &                     rm)
                     ENDDO
                     IF ( zmir(kk,1,IEXP).LT..5D0 )
     &                 zmir(kk,1,IEXP) = .5D0
                     IF ( zmir(kk,2,IEXP).LT..5D0 )
     &                 zmir(kk,2,IEXP) = .5D0
                  ENDIF
               ENDDO
            ENDDO
            DO kk = 1 , 6
               XIR(kk,IEXP) = XIR(kk,IEXP)*1.01D0
               DO kh = 1 , 8
                  MULTI(kh) = 0
                  LAMDA(kh) = 0
                  LDNUM(kh,2) = 0
                  LDNUM(kh,1) = 0
               ENDDO
               NMAX = 2
               ELM(1) = 1.
               LEAD(1,1) = 1
               LEAD(2,1) = 2
               SPIN(1) = 0.
               IFAC(1) = 1
               LAMMAX = 1
               MEMAX = 1
               MAGEXC = 0
               kkk = 0
               icg = 1
               IF ( ihlm(kk).NE.0 ) THEN
                  MULTI(kk) = 1
                  LAMDA(1) = kk
                  SPIN(2) = DBLE(kk)
                  IFAC(2) = 1
                  LDNUM(kk,1) = 1
                  icg = 1
                  CALL LOAD(IEXP,1,icg,0.D0,jj)
                  CALL LOAD(IEXP,2,icg,0.D0,jj)
                  CALL PATH(1)
                  sz1 = MIN(zmir(kk,1,IEXP),10.D0)
                  sz2 = zmir(kk,2,IEXP)/50.D0
C 2.4009604D-3 = 2 / 833
                  acof = 2.D0/833.D0/zmir(kk,2,IEXP)
C 8.163265D-4 = 1 / 1225
                  bcof = 1.D0/1225.D0
                  DO jd = 1 , jde
                     nksi = 5
                     IF ( jd.EQ.2 ) nksi = 10
                     IF ( MAGA(IEXP).EQ.0 ) nksi = 10
                     DO jk = 1 , 3
                        ZETA(jk) = 0.
                     ENDDO
                     nz = 50
                     IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) nz = 1
                     DO jk = 1 , nksi
                        XI(1) = XIR(kk,IEXP)*(jk-1)/(nksi-1)
                        IF ( jk.EQ.1 ) XI(1) = .02D0
                        s11 = 0.
                        s21 = 0.
                        s12 = 0.
                        s22 = 0.
                        ph1 = 0.
                        ph2 = 0.
                        DO jz = 1 , nz
                           ZETA(jd) = sz2*jz
                           IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) ZETA(jd)
     &                          = sz1
                           IF ( ZETA(jd).LT..1 ) INTERV(IEXP) = 1000
                           IF ( ZETA(jd).GE..1 ) INTERV(IEXP) = intvh
                           CALL ALLOC(ACCUR)
                           CALL SNAKE(IEXP,ZPOL)
                           CALL SETIN
                           CALL STING(1)
                           IF ( kk.GT.2 ) THEN
                              ARM(1,5) = (.9999999D0,0.)
                              ARM(2,5) = (1.2D-6,0.)
                              ARM(1,6) = (.9999998D0,0.)
                              ARM(2,6) = (.9D-6,0.)
                              DO kh = 1 , 4
                                 ARM(1,kh) = (-1.D-6,0.)
                                 ARM(2,kh) = (1.D-6,0.)
                              ENDDO
                           ENDIF
                           CALL INTG(IEXP)
                           jp = 2
                           IF ( MAGA(IEXP).NE.0 .AND. jd.EQ.2 ) jp = 3
                           p = DBLE(ARM(1,5))
                           r = DIMAG(ARM(1,5))
                           qr = DBLE(ARM(jp,5))
                           s = DIMAG(ARM(jp,5))
                           test = p*p + r*r + qr*qr + s*s
                           p = p/SQRT(test)
                           s = ABS(r/s)
                           IF ( jk.EQ.1 ) THEN
                              IF ( MAGA(IEXP).EQ.0 ) THEN
                                 q1 = 0.
                                 GOTO 1302
                              ELSEIF ( jd.EQ.2 .OR. MAGA(IEXP).EQ.0 )
     &                                 THEN
                                 q1 = 0.
                                 GOTO 1302
                              ENDIF
                           ENDIF
                           q1 = ARCTG(s,ph1)
                           ph1 = q1
 1302                      IF ( jk.EQ.1 ) THEN
                              IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) THEN
                                 q2 = 0.
                                 GOTO 1304
                              ENDIF
                           ENDIF
                           q2 = ARCCOS(p,ph2)
                           ph2 = q2
 1304                      q1 = q1/ZETA(jd)/2.
                           q2 = q2/ZETA(jd)
                           IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) q2 = -q2
                           IF ( jd.NE.1 .OR. MAGA(IEXP).EQ.0 ) THEN
                              s11 = s11 + q1
                              s12 = s12 + q1*jz
                              s21 = s21 + q2
                              s22 = s22 + jz*q2
                           ENDIF
                        ENDDO
                        IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) THEN
                           PARX(IEXP,2*kk-1,jk) = q1
                           PARX(IEXP,2*kk,jk) = q2
                        ELSE
                           PARXM(IEXP,1,jk,kk) =
     &                       acof*(2.D0*s12-51.D0*s11)
                           PARXM(IEXP,2,jk,kk) =
     &                       bcof*(101.D0*s11-3.D0*s12)
                           PARXM(IEXP,3,jk,kk) =
     &                       acof*(2.D0*s22-51.D0*s21)
                           PARXM(IEXP,4,jk,kk) =
     &                       bcof*(101.D0*s21-3.D0*s22)
                        ENDIF
                     ENDDO ! Loop over jk
                  ENDDO ! Loop over jd
               ENDIF
            ENDDO ! Loop over kk
            EMMA(IEXP) = emhl1
            NMAX = nmaxh
            SPIN(1) = sh1 ! Restore ground-state spin
            SPIN(2) = sh2 ! Restore spin of first excited state
            IFAC(1) = ih1
            IFAC(2) = ih2
            MAGEXC = magh
            ISO = isoh
            ELM(1) = eh1
            LEAD(1,1) = lh1
            LEAD(2,1) = lh2
            LAMMAX = lamh
            MEMAX = memh
            DO kh = 1 , 8 ! For each multipolarity
               LDNUM(kh,2) = ihlm(kh+24)
               MULTI(kh) = ihlm(kh)
               LAMDA(kh) = ihlm(kh+8)
               LDNUM(kh,1) = ihlm(kh+16)
            ENDDO
            INTERV(IEXP) = intvh
         ENDDO ! Loop over experiments jexp

         irix = 7
         IF ( IBPS.EQ.1 ) irix = 27 ! unit 7 for particle 1, 27 for particle 2
         REWIND irix
         DO iuy = 1 , 6
            WRITE (irix,*) (XIR(iuy,jj),jj=1,NEXPT)
            WRITE (irix,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT ! For each experiment
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  WRITE (irix,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               WRITE (irix,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
C---- gosia2 changes start
C        Added for gosia2
         irix = 12
         IF ( IBPS.EQ.1 ) irix = 32
         REWIND irix
         DO kuku = 1 , MEMAX
            WRITE (irix,*) ELM(kuku)
         ENDDO
C---- gosia2 changes end
         DO jj = 1 , 2
            DO jj1 = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
               IDIVE(jj1,jj) = 1
            ENDDO
         ENDDO
      ELSE ! iobl .lt. 1
         irix = 7
         IF ( IBPS.EQ.1 ) irix = 27 ! unit 7 for particle 1, 27 for particle 2 (gosia2 only)
         REWIND irix
         DO iuy = 1 , 6
            READ (irix,*) (XIR(iuy,jj),jj=1,NEXPT)
            READ (irix,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT ! For each experiment
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  READ (irix,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               READ (irix,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
         DO jgs = 1 , MEMAX ! For each matrix element
            DO jgr = 1 , 7
               QAPR(jgs,1,jgr) = 0.
            ENDDO
         ENDDO
      ENDIF

C     Handle map
      IF ( IPRM(12).NE.0 .OR. op2.EQ.'MAP ' ) THEN ! changed for gosia2
         IPRM(12) = 0
         DO jex = 1 , NEXPT
            DO lex = 1 , 6
               IF ( MULTI(lex).NE.0 ) THEN
                  WRITE (22,99040) jex , XIR(lex,jex)
99040             FORMAT (1X//30X,'EXPERIMENT',1X,1I2,10X,'MAX.XI=',
     &                    1F6.4)
                  WRITE (22,99041) lex , zmir(lex,2,jex)
99041             FORMAT (1X/30X,'E',1I1,8X,'MI=0',5X,'MAX.ZETA=',
     &                    1F6.3//)
                  WRITE (22,99054)
                  DO kex = 1 , 10
                     xxi = XIR(lex,jex)*(kex-1)/9.D0
                     WRITE (22,99055) xxi ,
     &                                (PARXM(jex,ilx,kex,lex),ilx=1,4)
                  ENDDO
                  IF ( MAGA(jex).NE.0 ) THEN
                     WRITE (22,99042) lex , zmir(lex,1,jex)
99042                FORMAT (1X//30X,'E',1I1,8X,'MI=+/-1',5X,
     &                       'MAX.ZETA=',1F6.3//)
                     WRITE (22,99054)
                     DO kex = 1 , 5
                        xxi = XIR(lex,jex)*(kex-1)/4.D0
                        u = 0.
                        WRITE (22,99055) xxi , u , PARX(jex,2*lex-1,kex)
     &                         , u , PARX(jex,2*lex,kex)
                     ENDDO ! Loop on kex
                  ENDIF ! if maga(jex).ne.0
               ENDIF ! if multi(lex).ne.0
            ENDDO ! Loop on lex
         ENDDO ! Loop on jex
      ENDIF ! IPRM(12).ne.0 .or. op2.eq.'MAP '
      IF ( op2.NE.'GOSI' .AND. op2.NE.'ERRO' ) GOTO 100 ! Back to input loop
      IF ( op2.EQ.'ERRO' ) GOTO 400

 1400 mret = ABS(mret) ! Added for gosia2
      DO kh1 = 1 , MEMAX
         HLM(kh1) = ELM(kh1)
      ENDDO
      lfagg = 0
      DO kh1 = 1 , MEMAX
         IVAR(kh1) = ivarh(kh1)
      ENDDO
C---- gosia2 changes start
      IF ( nawr.EQ.1 ) THEN
         MCFIX = 0 ! Calculate chisq using CNOR1 not CNOR
         CALL MINI(chisq,chiok,nptl,conu,imode,idr,xtest,0,0,0,bten)

         IF ( IBPS.EQ.0 ) THEN
            JZB = 25
            mrepf = 1
         ELSE ! i.e. IBPS is 1
            IBPS = IBPS + 1
            IF ( IBPS.EQ.1 ) JZB = 26 ! But IBPS must be 2 here! Ooops!
            IF ( IBPS.EQ.1 ) GOTO 2200
         ENDIF
         IF ( IPS1.EQ.0 ) GOTO 2000
         IMIN = IMIN + 1
         DO iva = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
            JSKIP(iva) = 1
         ENDDO
         irix = 12
         IF ( IBPS.EQ.2 ) irix = 32 ! unit 12 for particle 1, 32 for particle 2
         REWIND irix
         DO lkj = 1 , MEMAX
            WRITE (irix,*) ELM(lkj)
         ENDDO
         IF ( ifm.EQ.1 ) CALL PRELM(3)
         IF ( ifm.NE.1 ) GOTO 70004 ! changed for gosia2
         GOTO 2000
      ENDIF

C     Calculate chi squared and normalization without minimizing
      MCFIX = 1
C---- gosia2 changes end
      CALL MINI(chisq,1.D+38,nptl,conu,imode,idr,xtest,0,0,0,bten)
C---- gosia2 changes start

C     Set CNOR1 to the weighted average of CNOR1 and CNOR2
      DO kh1 = 1 , LP6 ! LP6 = 32 (maximum number of gamma detectors)
         DO kh2 = 1 , LP3 ! LP3 = 75 (maximum number of levels)
            IF ( JZB.EQ.25 ) THEN ! If it is the second nucleus
               CNOR1(kh1,kh2) = CNOR(kh1,kh2)
            ELSE
               CNOR2(kh1,kh2) = CNOR(kh1,kh2)
               CNOR1(kh1,kh2) = CNOR1(kh1,kh2)*(1.0-WNOR)+
     &           CNOR2(kh1,kh2)*WNOR
            ENDIF
         ENDDO
      ENDDO

C     Increment iteration counter
      mawr = mawr + 1

C     Store chisq
      IF ( JZB.EQ.25 ) ccch1 = chisq
      IF ( JZB.EQ.26 ) ccch2 = chisq
      WRITE (*,*) 'ITER = ' , mawr , ' CHISQ1 = ' , ccch1 ,
     &  ' CHISQ2 = ' , ccch2
      chir = ccch1 + ccch2

C     Decide whether to terminate here
      IF ( JZB.EQ.26 ) THEN
         cht = ABS(chir-chp) ! Difference in old and new chisq
         WRITE (*,*) 'CHISQ DIFF = ' , chir , ' - ' , chp ,
     &     ' = ' , cht
         IF ( (cht.LT.0.00001 .AND. mawr.GE.2) .OR. mawr.GE.100 )
     &     mret = 0 ! terminate
         nawr = 1
      ENDIF

      IF ( IBPS.EQ.0 .AND. JZB.EQ.25 ) JZB = 26
      IF ( IBPS.NE.0 ) JZB = 25
      GOTO 2200 ! Back to beginning
C---- gosia2 changes end

C--------------------------------
C OP,MINI
70009 READ (JZB,*) imode , nptl , chiok , conu , xtest , LOCKF ,
     &  NLOCK , IFBFL , LOCKS , DLOCK
      op2 = opcja
      IMIN = IMIN + 1
      IF ( IMIN.NE.1 ) GOTO 1400
      GOTO 70008 ! End of OP,MINI

C--------------------------------
C OP,RAND
70010 READ (JZB,*) SE ! Seed for random number generator
      CALL MIXUP
      WRITE (22,99007)
99007 FORMAT (1X///5X,'MATRIX ELEMENTS RANDOMIZED...'///)
      CALL PRELM(2)
      GOTO 100 ! End of OP,RAND - back to input loop

C--------------------------------
C OP,RAW
C     Read absorber coefficients from unit 8
70011 REWIND 8
      DO l = 1 , 8
        READ (8,*) (ABC(l,j),j=1,10) ! Absorption coefficients
        DO j = 1 , 10
          ABC(l,j) = LOG(ABC(l,j))
        ENDDO
      ENDDO
      DO l = 1 , nfd
        READ (8,*) (THICK(l,j),j=1,7) ! thickness of absorbers
      ENDDO
      DO l = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
        DO j = 1 , 200
          ICLUST(l,j) = 0
        ENDDO
        DO j = 1 , 20
          LASTCL(l,j) = 0
        ENDDO
        IRAWEX(l) = 0
      ENDDO

C     Read input from standard input
      DO l = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
        READ (JZB,*) mexl ! experiment number
        IF ( mexl.EQ.0 ) GOTO 100 ! Back to input loop
        IRAWEX(mexl) = 1
        n = NANG(mexl)
        DO j = 1 , n
          jj = ITMA(mexl,j) ! Get identity of detector
          READ (JZB,*) (AKAVKA(mexl,k,jj),k=1,8) ! efficiency curve parameters
          AKAVKA(mexl,9,jj) = ideff(mexl)
        ENDDO
        READ (JZB,*) kclust ! number of clusters
        IF ( kclust.NE.0 ) THEN
          DO j = 1 , kclust
            READ (JZB,*) numcl ! Number of detectors for this cluster
            READ (JZB,*) (liscl(k),k=1,numcl) ! Indices of logical detectors
            LASTCL(l,j) = liscl(numcl) ! Index of last detector in cluster
            DO k = 1 , numcl
              kk = liscl(k)
              ICLUST(l,kk) = j ! Set cluster number
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      GOTO 100 ! End of OP,RAW - back to input loop

C--------------------------------
C OP,RE
70012 DO jrls = 1 , MEMAX ! For each matrix element
         IF ( IVAR(jrls).NE.0 .OR. irfix.NE.1 ) THEN
            IF ( IVAR(jrls).GT.999 ) THEN
               IF ( jfre.EQ.1 ) GOTO 1100
            ENDIF
            IVAR(jrls) = 2
            ELML(jrls) = -ABS(ELML(jrls))
            ELMU(jrls) = ABS(ELMU(jrls))
            IF ( jrls.GT.MEMX6 ) IVAR(jrls) = 1
         ENDIF
 1100    CONTINUE
      ENDDO ! For each matrix element jrls
      DO jrls = 1 , MEMAX
         ivarh(jrls) = IVAR(jrls)
      ENDDO
      GOTO 100 ! Back to input loop

C--------------------------------
C OP,REST
70013 irix = 12
      IF ( IBPS.EQ.1 ) irix = 32 ! unit 12 for particle 1, 32 for particle 2 (gosia2 only)
      REWIND irix
      memax1 = MEMAX + 1
      DO lkj = 1 , MEMAX
        READ (irix,*) ELM(lkj)
      ENDDO
      DO lkj = 1 , memax1
        READ (JZB,*) lkj1 , xlk
        IF ( lkj1.EQ.0 ) GOTO 120
        IF ( mres1.EQ.0 ) ELM(lkj1) = xlk ! condition added for gosia2
      ENDDO
 120  WRITE (22,99008)
99008 FORMAT (1X///5X,'*****',2X,
     &  'RESTART-MATRIX ELEMENTS OVERWRITTEN',2X,'*****'///)
C---- gosia2 changes start
      DO lkj = 1 , MEMAX
        IF ( mres1.EQ.0 .AND. JZB.EQ.25 ) ELM25(lkj) = ELM(lkj)
        IF ( mres1.EQ.0 .AND. JZB.EQ.26 ) ELM26(lkj) = ELM(lkj)
      ENDDO
      IF ( JZB.EQ.26 ) mres1 = 1
C---- gosia2 changes end
      DO kk = 1 , MEMAX
        la = MLT(kk)
        IF ( ivari(kk).GE.10000 ) THEN
          kk1 = ivari(kk)/10000
          kk2 = ivari(kk) - 10000*kk1
          la1 = la
          IF ( kk2.GE.100 ) THEN
            la1 = kk2/100
            kk2 = kk2 - 100*la1
          ENDIF
          inx1 = MEM(kk1,kk2,la1)
C      ELML(KK)=ELML(INX1)*ELM(KK)/ELM(INX1)
C      ELMU(KK)=ELMU(INX1)*ELM(KK)/ELM(INX1)
          SA(kk) = ELM(kk)/ELM(inx1)
          IVAR(kk) = 1000 + inx1
          IF ( ELMU(kk).LE.ELML(kk) ) THEN
            elmi = ELMU(kk)
            ELMU(kk) = ELML(kk)
            ELML(kk) = elmi
          ENDIF
        ENDIF
      ENDDO
      CALL PRELM(4) ! changed from 2 to 4 for gosia2
      GOTO 100 ! End of OP,REST - back to input loop

C--------------------------------
C OP,SELE
70014 CALL SELECT
      GOTO 2000 ! End of execution

C--------------------------------
C OP,SIXJ
70015 DO k = 1 , 2
        l = 4*k
        DO j = 1 , 80
          ixj = j - 1
          DO ms = 1 , 5
            mend = 2*(ms-3) + ixj
            WRITE (14,*) WSIXJ(l,4,4,ixj,mend,ixj-4) ,
     &        WSIXJ(l,4,4,ixj,mend,ixj-2) ,
     &        WSIXJ(l,4,4,ixj,mend,ixj) ,
     &        WSIXJ(l,4,4,ixj,mend,ixj+2) ,
     &        WSIXJ(l,4,4,ixj,mend,ixj+4)
          ENDDO
        ENDDO
      ENDDO
      GOTO 2000 ! End of OP,SIXJ - normal end of execution

C--------------------------------
C OP,THEO
70016 CALL OP_THEO
      GOTO 100 ! End of OP,THEO - back to input loop

C--------------------------------
C OP,TITL
70017 READ (JZB,99009) (title(k),k=1,20)
99009 FORMAT (20A4)
      WRITE (22,99010) (title(k),k=1,20)
99010 FORMAT (10X,20A4/10X,100('-'))
      GOTO 100 ! End of OP,TITL - back to input loop

C--------------------------------
C OP,TROU
70018 ITS = 1 ! Create tape 18 flag
      READ (JZB,*) kmat , rlr
      GOTO 100 ! End of OP,TROU - back to input loop

C--------------------------------
C OP,YIEL
70019 CALL ADHOC(oph,idr,nfd,ntap,iyr)
      GOTO 100 ! End of OP,YIEL - back to input loop

C--------------------------------
99048 FORMAT (1X//50X,'CALCULATED YIELDS'//5X,'EXPERIMENT ',1I2,2X,
     &        'DETECTOR ',1I2/5X,'ENERGY ',1F10.3,1X,'MEV',2X,'THETA ',
     &        1F7.3,1X,'DEG'//5X,'NI',5X,'NF',5X,'II',5X,'IF',5X,
     &        'YIELD',5X,'NORMALIZED YIELD'/)
99049 FORMAT (4X,1I3,4X,1I3,3X,1F4.1,3X,1F4.1,3X,1E11.5,3X,1E11.5)
99050 FORMAT (1X///44X,'OVERALL')
99051 FORMAT (1X///43X,'DIAGONAL')
99052 FORMAT (6X,1I3,5X,1I3,4X,1I3,5X,1F10.5,2X,'(',1F10.5,' ,',1F10.5,
     &        ')')
99053 FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
99054 FORMAT (5X,'XI',13X,'Q1',22X,'Q2'///13X,'SLOPE',2X,'INTERCEPT',7X,
     &        'SLOPE',5X,'INTERCEPT'//)
99055 FORMAT (2X,1F6.4,3X,1E8.2,2X,1E8.2,6X,1E8.2,2X,1E8.2)
99057 FORMAT (1X,/,2X,'DURING THE INTEGRATION',1X,
     &  'IT WAS NECESSARY TO SWITCH OFF THE TIME-OF-FLIGHT CORRECTION '
     &  'FOR EXPT ',I3,' LEVEL ', I3, ' TAU=',1E9.4)
      END
