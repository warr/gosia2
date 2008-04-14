C      PROGRAM GOSIA2(INPUT,OUTPUT,TAPE6=OUTPUT,TAPE3,TAPE14,TAPE1,
C        *TAPE15,TAPE17,TAPE18,TAPE4,TAPE2,TAPE7,TAPE9)
C************************
C HANDLES SIMULTANEOUS PROJECTILE/TARGET EXCITATION
C   THIS IS A 32-GE DETECTOR 386+WEITEK SCALAR VERSION FROM MAR., 1989
C FULL Q MAPS CALCULATED FOR E1 THROUGH E4 - UPDATE FROM AUG. 1989
C NOTE: THERE IS NO LONGER THE DOMINANT MULTIPOLARITY SELECTED
C FULL MAPS AND FITTING OF E5,E6 ADDED - APRIL 1990
C XI AND ZETA RANGES CALCULATED FOR EACH MULTIPOLARITY INDEPENDENTLY
C
C May 1995 additions - matrix elements generator following general
C structure of matrix elements as given in Bohr-Mottelson book.
C Generator is activated by OP,THEO
C Using OP,THEO the starting values are set according to the model while
C all couplings are kept as given using ME suboption. OP,THEO writes
C matrix elements to restart file which should be used subsequently.
C
C Eh. Vogt's addition - OP,FILE - allows to specify OPEN statements,
C if convenient.
C
C November 1990 update - number of levels = 75, yields=32 x 1500
C magnetic substates = 600 matrix elements=500
C April 1991 update - OP,RAW added - handles non-efficiency-corrected
C spectra, allows to form Ge detector clusters. Up to 20 clusters
C can be defined. Number of physical Ge's increased to 200, while
C number of datasets (i.e. single detector+cluster gamma yields)
C is still limited to 32
C Output is written on unit 22 to avoid mixing it with system
C messages on some computers
C PIN diode option added -Sept. 96
C July 1997 - known matrix elements input in OP,YIEL extended to
C all multipolarities. Note change in the input, now:
C LAMDA, NS1, NS2, ME, DME
C************************
C            GOSIA HAS BEEN DEVELOPED BY T.CZOSNYKA,D.CLINE AND C.Y.WU
C            AT NUCLEAR STRUCTURE RESEARCH LABORATORY,UNIV. OF ROCHESTER
C
C            SOME CONCEPTS USED COME FROM 1978 WINTHER/DE BOER CODE
C            C O U L E X AND FROM NSRL DEEXCITATION CODE C E G R Y.
C            HOWEVER,THE PARTS TAKEN FROM BOTH CODES ARE IN MOST
C            CASES COMPLETELY REWRITTEN,SO THE SIMILARITY OF
C            VARIABLE AND ROUTINE NAMES MAY BE MISLEADING.
C
C            VALUABLE CONTRIBUTIONS WERE ALSO MADE BY:
C            L.HASSELGREN AND A.BACKLIN  (UPPSALA)
C            J.SREBRNY  (WARSAW)
C            B.KOTLINSKI  (WARSAW AND ROCHESTER)
C
C
C            FOR INFORMATION,PLEASE CONTACT:
C            TOMASZ CZOSNYKA,SLCJ, WARSAW UNIVERSITY, WARSZAWA, POLAND
C            02-297 WARSZAWA, BANACHA 4-----PHONE (22)-222-123
C            DOUGLAS CLINE,DEPARTMENT OF PHYSICS,THE UNIVERSITY OF ROCHESTER
C            ROCHESTER,NY14627,USA------PHONE(585)275-4943
C
C            REF.----   UR/NSRL REPORT 308/1986
C
C************************  VERSION FROM SEPT, 2005  *******************
C
C**********************************************************************
C PROGRAM GOSIA2
C
C Calls: ADHOC, ALLOC, ANGULA, ARCCOS, ARCTG, CMLAB, COORD, DECAY, DJMM,
C        EFFIX, ELMT, FAKP, FHIP, FTBM, INTG, KLOPOT, KONTUR, LAGRAN, LOAD,
C        MINI, MIXR, MIXUP, OPENF, PATH, PRELM, PTICC, QFIT, READY, SETIN,
C        SIMIN, SNAKE, STING, TACOS, TAPMA, TEMB, TENS, WSIXJ, WTHREJ
C
C Uses global variables:
C      ABC    - absorption coefficients
C      ACCA   - accuracy
C      ACCUR  - accuracy required
C      AGELI  - angles of Ge detectors
C      AKAVKA - efficiency curve parameters
C      ARM    - reduced matrix elements
C      AVJI   - average J (N.B. here it is G(1))
C      B      - table of factorials
C      BEQ    - beta
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
C      DSIGS  -
C      DYEX   - error on experimental yield
C      EAMX   - known matrix elements and their errors
C      ELM    - matrix elements
C      ELM25  - matrix elements storage for target (for gosia2)
C      ELM26  - matrix elements storage for projectile (for gosia2)
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
C      HLM    -
C      HLMLM  -
C      IAMX   -
C      IAMY   -
C      IAX    - axial symmetry flag
C      IBPS   - target/projectile switch
C      ICLUST -
C      ICS    - read internal correction factors flag (OP,CONT switch CRF,)
C      IDIVE  - number of subdivisions
C      IDRN   - normalising transition for yields
C      IEXP   - experiment number
C      IFAC   -
C      IFBFL  - calculate derivatives with forward-backward method
C      IFMO   - include correction to angular distance for finite recoil distance.
C      ILE    -
C      IMIN   -
C      INHB   - inhibit error flag (LERF) setting in POMNOZ
C      INNR   - independent normalisation switch (see OP,CONT INR,)
C      INTERV - default accuracy check parameter for Adams-Moulton (see OP,CONT:INT)
C      INTR   -
C      IP     - table of prime numbers
C      IPRM   - various flags to control output
C      IPS1   - terminate after calculating and storing internal correction factors
C      IRAWEX - flag to indicate raw uncorrected yield
C      ISEX   -
C      ISKIN  - kinematic flag (0,1)
C      ISMAX  -
C      ISO    -
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      ITTE   - thick target experiment flag
C      IUNIT3 - unit for TAPE3
C      IVAR   - indicates a limit or correlation is set
C      IWF    -
C      IY     - index for yields
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of non-investigated nucleus
C      JENTR  -
C      JSKIP  - Experiments to skip during minimisation.
C      JZB    - unit to read from
C      KFERR  - error flag for minimization
C      KSEQ   - index of level
C      KVAR   -
C      LAMAX  - number of multipolarities to calculate
C      LAMBDA - list of multipolarities to calculate
C      LASTCL -
C      LDNUM  - number of matrix elements with each multipolarity populating each level
C      LEAD   - pair of levels involved in each matrix element
C      LIFCT  - index for lifetimes
C      LMAX   -
C      LMAXE  - maxmium multipolarity needed for calculation
C      LNORM  - normalisation constant control
C      LNY    - use logs to calculate chi squared
C      LOCKF  - flag to fix matrix elements with most significant derivative
C      LOCKS  - lock flag. If LOCKS=1, fix at first stage of minimization
C      LP1    - maximum number of experiments (50)
C      LP10   - maximum number of reduced matrix elements (600)
C      LP11   - LP8 - 1 (103)
C      LP12   - number of steps of omega (365)
C      LP13   - LP9 + 1 (47901)
C      LP14   - maximum space for collision functions (4900)
C      LP2    - maximum number of matrix elements (500)
C      LP3    - maximum number of levels (75)
C      LP4    - maximum number of yields (1500)
C      LP6    - maximum number of gamma detectors (32)
C      LP7    - start of collision functions (45100)
C      LP8    - (104)
C      LP9    - length of ZETA - 2100 (47900)
C      MAGA   - number of magnetic substates in approximate calculation
C      MAGEXC -
C      MCFIX  - flag to say if normalization coefficient is fixed (for gosia2)
C      MEMAX  - number of matrix elements
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      MULTI  - number of matrix elements having given multipolarity
C      NAMX   - number of known matrix elements
C      NANG   - number of gamma-ray detectors for each experiment
C      NBRA   - number of branching ratios
C      NCM    - calculate kinematics assuming this spin for final state (default = 2.0)
C      NDIM   - maximum number of levels
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NLOCK  - number of elemnts to fix if LOCKF=1
C      NMAX   - number of levels
C      NMAX1  -
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      PARX   - [for maps]
C      PARXM  - [for maps]
C      POWER  - x (N.B. here it is G(7))
C      QAPR   - approximate Coulomb amplitudes
C      SA     - ratio of matrix elements for correlated elements
C      SE     -
C      SPIN   - spin of level
C      SUBCH1 -
C      SUBCH2 -
C      SUMCL  -
C      SWG    -
C      TAU    - 
C      THICK  - thickness of each absorber type
C      TIMEC  - Tau_C (N.B. here it is G(4))
C      TIMEL  - lifetimes and their errors
C      TLBGD  - theta of particle detector
C      TREP   - theta of recoiling nucleus
C      UPL    - upper limits for all gamma detectors
C      VINF   - speed of projectile at infinty
C      XA     - A of investigated nucleus
C      XA1    - A of non-investigated nucleus
C      XI     - xi coupling coefficients
C      XIR    - [for maps]
C      XLAMB  - Lambda* (N.B. here it is G(3))
C      XV     - energy meshpoints where we calculate exact Coulex
C      YEXP   - experimental yields
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      YGP    - gamma yield calculated with correction to angular distribution from finite recoil distance
C      YNRM   - relative normalisation for gamma detectors
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C      ZETA   - various coefficients
C      ZPOL   -
C      ZV     - energy meshpoints

      PROGRAM GOSIA2
      IMPLICIT NONE
      REAL*8 ABC , ACCA , ACCUR , acof , AGELI , AKAVKA , AKS , ap , 
     &       ARCCOS , ARCTG , arg , ax , B , bcof , be2 , be2a , be2b , 
     &       be2c , BEQ , BETAR
      REAL*8 bk , bl , bm , bmx , BRAT , bten , bu , CAT , CC , ccc , 
     &       ccch1 , ccch2 , ccd , cf , chilo , chiok , chir , chis0 , 
     &       chisl , chisq
      REAL*8 chiss , chp , cht , CNOR , CNOR1 , CNOR2 , cnst , cocos , 
     &       conu , CORF , d , decen , dedx , DELTA , DEVD , DEVU , 
     &       DIPOL , DIX , DLOCK , DQ
      REAL*8 DS , dsd , DSE , DSG , dsig , DSIGS , dst , dsx , dsxm , 
     &       DYEX , EAMX , effi , EG , eh1 , ELM , ELM25 , ELM26 , 
     &       ELMH , elmi , ELML
      REAL*8 ELMT , ELMU , emhl1 , EMMA , emn , emx , EN , enb , ENDEC , 
     &       eng , enh , ENZ , EP , EPS , EROOT , esd , esp , ess , 
     &       fi0 , fi1
      REAL*8 fic , FIEX , fiex1 , figl , fipo1 , fm , G , GRAD , gth , 
     &       hen , het , HLM , HLMLM , ODL , p , PARX , PARXM , pfi , 
     &       ph1 , ph2
      REAL*8 pi , PILOG , po1 , po2 , polm , pop1 , pr , pv , Q , q1 , 
     &       q2 , QAPR , qc , QCEN , qfac , qr , qui , r , r1 , r2
      REAL*8 r3 , r4 , rem , remax , rl , rlr , rm , rx , ry , rz , s , 
     &       s11 , s12 , s21 , s22 , SA , sbe , SE , sf , SGW
      REAL*8 sh , sh1 , sh2 , SIMIN , slim , SPIN , SUBCH1 , SUBCH2 , 
     &       SUMCL , summm , sz1 , sz2 , TACOS , TAU , tau1 , tau2 , 
     &       test , TETACM , tetrc , tfac
      REAL*8 thc , THICK , TIMEL , title , TLBDG , tmn , tmx , todfi , 
     &       TREP , tta , tth , tting , ttttt , txx , u , UPL , VACDP , 
     &       val , VINF , waga
      REAL*8 wph , wpi , WSIXJ , wth , wthh , WTHREJ , XA , XA1 , xep , 
     &       XI , xi1 , xi2 , XIR , xk1 , xk2 , xl1 , xlevb , xlk , 
     &       xm1 , xm2
      REAL*8 xm3 , XNOR , xtest , XV , xw , xx , xxi , ycorr , YEXP , 
     &       YGN , YGP , YNRM , YV , yy , yyd1 , yydd , yyy , ZETA , 
     &       zmir , zp
      REAL*8 ZPOL , ZV , zz
      INTEGER*4 i , i122 , IAMX , IAMY , IAPR , iapx , IAXS , ib , 
     &          ibaf , IBPS , IBRC , IBYP , icg , icll , ICLUST , ICS , 
     &          ict , ictl , id , idf
      INTEGER*4 IDIVE , idr , IDRN , iecd , ient , IEXP , IFAC , IFBFL , 
     &          ifbp , ifc , ifm , IFMO , ifwd , ig1 , ig2 , ih1 , ih2 , 
     &          ihlm , ihuj , ii
      INTEGER*4 ij , ija0 , ijaja , ijan , ijk , ijx , ILE , ile1 , 
     &          ilevls , ilx , im , IMIN , imode , in1 , in2 , inclus , 
     &          ind , ind1 , ind2 , indx
      INTEGER*4 INHB , inko , inm1 , inm2 , inn , INNR , inpo , intend , 
     &          INTERV , INTR , intvh , inva , inx1 , iobl , iocc , 
     &          iopri , iosr , IP , IPATH , ipd
      INTEGER*4 iph , IPI , ipine , ipinf , ipo1 , ipo2 , ipo3 , ipp , 
     &          iprc , ipri , IPRM , IPS1 , IRAWEX , irea , irep , 
     &          irfix , irix , ISEX , isip , iske
      INTEGER*4 iskf , ISKIN , isko , iskok , ISMAX , ISO , isoh , 
     &          ispa , ispb , ITMA , itno , itp , ITS , ITTE , IUNIT3 ,
     &          iuy , iva , iva1 , IVAR , ivarh , ivari
      INTEGER*4 ivrh , IWF , ixj , ixl , ixm , IY , iyr , IZ , IZ1 , 
     &          izcap , j , ja , jan , jan1 , jb , jb1 , jb2 , jd , 
     &          jde , jdy
      INTEGER*4 je , JENTR , jex , jexp , jfi , jfre , jgd , jgl , 
     &          jgl1 , jgr , jgs , jj , jj1 , jjjj , jjlx , jjx , jk , 
     &          jkloo , jkmx , jkmy
      INTEGER*4 jktt , jl , jmm , jmpin , jp , jphd , jpin , jrls , js , 
     &          JSKIP , jt , jtp , jyi , jyi1 , jyi2 , jyv , jz , JZB , 
     &          k , kb
      INTEGER*4 kclust , kerf , kex , KF , KFERR , kh , kh1 , kh2 , kk , 
     &          kk1 , kk2 , kkk , kkkk , kl , kloop , kmat , kq , KSEQ , 
     &          ktt , kuku
      INTEGER*4 KVAR , l , la , la1 , lam , lamd , LAMDA , lamh , 
     &          LAMMAX , LASTCL , lb , lck1 , lck2 , LDNUM , LEAD , 
     &          LERF , levl , lex , lexp , lfagg
      INTEGER*4 lfini , lh1 , lh2 , LIFCT , liscl , lkj , lkj1 , ll , 
     &          lli , lll , LMAX , lmax1 , LMAXE , lmaxh , LNORM , LNY , 
     &          locat , LOCKF , LOCKS , loct
      INTEGER*4 lp0 , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , lpin , ltrn , 
     &          ltrn1 , ltrn2 , lu , lx
      INTEGER*4 lxd , LZETA , m25 , m26 , MAGA , MAGEXC , magh , mawr , 
     &          MCFIX , mdupa , MEM , MEMAX , memax1 , memh , memx4 , 
     &          MEMX6 , mend , mexl , mfla , mlt
      INTEGER*4 mm , mmmm , mpin , mrepf , mres1 , mret , ms , MULTI , 
     &          n , na , na1 , naa , nallow , NAMX , NANG , nawr , 
     &          naxfl , nb1 , nb2 , nbands
      INTEGER*4 NBRA , nch , NCM , NDIM , ndima , NDST , ndum , ne , 
     &          NEXPT , nf , nfd , nfdd , nfi , nflr , nft , nged , 
     &          ngpr , ni , NICC , nksi
      INTEGER*4 nl , NLOCK , NMAX , NMAX1 , nmaxh , nmemx , nnl , 
     &          nogeli , npce , npce1 , npct , npct1 , npt , nptl , 
     &          nptx , ns1 , ns2 , ntap , ntt , numcl
      INTEGER*4 nval , NYLDE , nz
      LOGICAL ERR
      COMPLEX*16 ARM , EXPO
      CHARACTER*4 oph , op1 , opcja , op2
      CHARACTER*1 prp
      DIMENSION ihlm(32) , esp(20) , dedx(20) , bten(1200) , 
     &          fiex1(11,20,2) , title(20) , pfi(101) , zmir(6,2,50) , 
     &          iecd(50) , wpi(11,2) , tau1(10) , eng(10) , tau2(10,7) , 
     &          xl1(7) , qui(8,10) , cf(8,2) , ivarh(500) , liscl(200) , 
     &          dsxm(100,20,20) , levl(50) , xlevb(50,2) , bm(8,20,20,3)
     &          , mlt(500) , ivari(500) , jpin(50)
      COMMON /CLUST / ICLUST(50,200) , LASTCL(50,20) , SUMCL(20,500) , 
     &                IRAWEX(50)
      COMMON /CCCDS / NDST(50)
      COMMON /INHI  / INHB
      COMMON /IDENT / BEQ
      COMMON /EFCAL / ABC(8,10) , AKAVKA(8,200) , THICK(200,7)
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      COMMON /BREC  / BETAR(50)
      COMMON /ADBXI / EXPO(500)
      COMMON /DIMX  / DIX(4) , ODL(200)
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /CINIT / CNOR(32,75) , INNR , MCFIX , CNOR1(32,75) , 
     &                CNOR2(32,75)
      COMMON /XRA   / SE
      COMMON /HHH   / HLM(500)
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      COMMON /LIFE1 / LIFCT(50) , TIMEL(2,50)
      COMMON /DFTB  / DEVD(500) , DEVU(500)
      COMMON /ERRAN / KFERR
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /SECK  / ISKIN(50)
      COMMON /VLIN  / XV(51) , YV(51) , ZV(20) , DSG(20) , DSE(20) , DS
      COMMON /DUMM  / GRAD(500) , HLMLM(500) , ELMH(500)
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /YTEOR / YGN(500) , YGP(500) , IFMO
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /MAP   / PARX(50,12,5) , PARXM(50,4,10,6) , XIR(6,50)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      COMMON /GGG   / G(7)
      COMMON /AZ    / ARM(600,7)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CXI   / XI(500)
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /MINNI / IMIN , LNORM(50)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /PRT   / IPRM(20)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CB    / B(20)
      COMMON /CLM   / LMAX
      COMMON /CLCOM0/ IFAC(75)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /CLCOM9/ ERR
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CEXC9 / INTERV(50)
      COMMON /CAUX0 / EMMA(75) , NCM
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /WARN  / SGW , SUBCH1 , SUBCH2 , IWF
      COMMON /THTAR / ITTE(50)
      COMMON /FIT   / LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
      COMMON /APRX  / LERF , IDIVE(50,2)
      COMMON /SKP   / JSKIP(50)
      COMMON /TRB   / ITS
      COMMON /SEL   / KVAR(500)
      COMMON /ERCAL / JENTR , ICS
      COMMON /LOGY  / LNY , INTR , IPS1
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      COMMON /SWITCH/ JZB , IBPS , IUNIT3
      COMMON /RESC  / ELM25(500) , ELM26(500) ! For gosia2
      DATA (eng(k),k=1,10)/.05 , .06 , .08 , .1 , .15 , .2 , .3 , .5 , 
     &      1. , 1.5/
C     Absorption coefficients in units of 1/cm for Ge
      DATA (tau1(k),k=1,10)/17.656 , 10.726 , 5.076 , 2.931 , 1.3065 , 
     &      .8828 , .5959 , .4357 , .3041 , .2472/
C     Absorption coefficients in units of 1/cm for Al, C, Fe, Cu, Ag/Cd/Sn, Ta
C     and Pb at the energies 0.05, 0.06, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5, 1, 1.5
C     MeV
      DATA (tau2(k,1),k=1,10)/.9883 , .7473 , .5442 , .4592 , .3718 , 
     &      .3302 , .2814 , .2278 , .1657 , .1350/
      DATA (tau2(k,2),k=1,10)/1.014 , .7443 , .5195 , .4261 , .3362 , 
     &      .2967 , .2518 , .2038 , .1479 , .1204/
      DATA (tau2(k,3),k=1,10)/15.167 , 9.405 , 4.652 , 2.889 , 1.525 , 
     &      1.135 , .8643 , .6592 , .4703 , .3830/
      DATA (tau2(k,4),k=1,10)/23.184 , 14.182 , 6.777 , 4.059 , 1.970 , 
     &      1.384 , .9936 , .7473 , .5274 , .4297/
      DATA (tau2(k,5),k=1,10)/84.351 , 51.445 , 23.822 , 13.070 , 
     &      4.774 , 2.605 , 1.339 , .7925 , .5005 , .4032/
      DATA (tau2(k,6),k=1,10)/93.364 , 58.559 , 125.96 , 70.713 , 
     &      25.302 , 12.541 , 5.193 , 2.215 , 1.077 , .8176/
      DATA (tau2(k,7),k=1,10)/89.809 , 56.338 , 27.009 , 62.966 , 
     &      22.933 , 11.334 , 4.540 , 1.813 , .8020 , .5900/

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
      lp0 = 50000 ! Size of ZETA array
      LP1 = 50 ! Maximum number of experiments
      LP2 = 500 ! Maximum number of matrix elements
      LP3 = 75 ! Maximum number of levels
      LP4 = 1500
      LP6 = 32
      LP7 = lp0 - 4900 ! Start of collision coefficients in ZETA
      LP8 = LP3*28 + 1
      LP9 = lp0 - LP3*28
      LP10 = 600
      LP11 = LP8 - 1
      LP12 = 365 ! Maximum number of steps of omega (dimension of ADB, SH, CH)
      LP13 = LP9 + 1
      LP14 = 4900 ! Maximum number of collision coefficients

      JZB = 5

C---- gosia2 changes start
C     Read target/projectile switch and first option from standard input
      READ (*,*) IBPS
      READ (*,99001) op1 , op2

C     Open files
      CALL OPENF1
      nawr = 0

C     Use input unit 25 for target and 26 for projectile
      IBPS = IBPS - 1
      JZB = 25
      IF ( IBPS.EQ.1 ) JZB = 26
C---- gosia2 changes end

C     Initialize normalization to 1.
      DO i = 1 , LP3 ! LP3 = 75
        DO j = 1 , LP6 ! LP6 = 32
            CNOR(j,i) = 1.
            CNOR1(j,i) = 1. ! Added for gosia2
            CNOR2(j,i) = 1. ! Added for gosia2
         ENDDO
      ENDDO

C---- gosia2 changes start
      IUNIT3 = 33
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
      BEQ = -983872.
      ipinf = 0
      iyr = 0
      pi = 3.141592654
      INNR = 0
      itno = 0
      chisq = 0.
      chilo = 0.
      IWF = 1
      ifm = 0 ! Fast minimisation switch off by default
      IPS1 = 11
      ifwd = -1
      INTR = 0
      LNY = 0
      JENTR = 0
      ICS = 0

      DO i = 1 , LP1 ! LP1 = 50
         jpin(i) = 0
         iecd(i) = 0
      ENDDO

      txx = 0.
      SGW = 3.
      SUBCH1 = 0.
      SUBCH2 = 0.
      ITS = 0 ! Create tape 18 flag
      iosr = 0
      LOCKS = 0
      DLOCK = 1.1
      kerf = 0
      IFBFL = 0
      NLOCK = 0
      LOCKF = 0
      DO i = 1 , LP4 ! LP4 = 1500
         DO j = 1 , LP6 ! LP6 = 32
            CORF(i,j) = 1.
         ENDDO
      ENDDO
      DO i = 1 , 20
         IPRM(i) = 1
         DO j = 1 , 5
            CC(i,j) = 0.
         ENDDO
      ENDDO
      IPRM(2) = 0 ! changed from 1 to 0 for gosia2
      IPRM(4) = 0 ! changed from -2 to 0 for gosia2
      IPRM(5) = 1 ! changed from 11111 to 1 for gosia2
      IPRM(6) = 11111
      IPRM(7) = 0
      IPRM(11) = 0 ! changed from 1 to 0 for gosia2
      IPRM(12) = 0 ! changed from 1 to 0 for gosia2
      IPRM(14) = 0 ! changed from 1 to 0 for gosia2
      IPRM(16) = 0
      IPRM(17) = 0
      IPRM(18) = 0
      IPRM(19) = 0
      IPRM(20) = 0
      DO i = 1 , LP1 ! LP1 = 50
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
      DO k = 1 , LP1 ! LP1 = 50
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
      NDIM = LP3 ! LP3 = 75
      ISO = 1
      B(1) = 1.
      DO i = 2 , 20
         B(i) = B(i-1)*(i-1)
      ENDDO
      LMAXE = 0
      CALL FAKP
      CALL FHIP
      NCM = 2 ! Default final spin for kinematics calculation (OP,CONT NCM,)
      DO ijx = 1 , LP1 ! LP1 = 50
         INTERV(ijx) = 1
      ENDDO
      la = 0
      ipo3 = 1
      indx = 0
      ACCUR = .00001
      icg = 1
      ient = 1
      jphd = 1
      DIPOL = 0.005
      MAGEXC = 0
      LAMMAX = 0
      DO lam = 1 , 8
         DO lexp = 1 , LP3 ! LP3 = 75
            LDNUM(lam,lexp) = 0
         ENDDO
         MULTI(lam) = 0
         LAMDA(lam) = 0
      ENDDO
      DO j = 1 , LP2 ! LP2 = 500
         EXPO(j) = (1.,0.)
         KVAR(j) = 1
         ELM(j) = 0.
      ENDDO
      DO j = 1 , LP1 ! LP1 = 50
         JSKIP(j) = 1
         ISKIN(j) = 0
      ENDDO
      DO j = 1 , LP3 ! LP3 = 75
         ISEX(j) = 1111
      ENDDO
      ISEX(1) = 0
      ACCA = .00001
      oph = '    '
      nmemx = LP2 + 9 ! LP2 = 500
      IEXP = 1
      IMIN = 0
      i122 = 0
      DO j = 1 , LP2 ! LP2 = 500
         DO k = 1 , 2
            DO l = 1 , 7
               QAPR(j,k,l) = 0.
            ENDDO
         ENDDO
      ENDDO
      ERR = .FALSE.
      intend = 0 ! End of initialization

C     Start reading input file.
 100  READ (JZB,99001) op1 , op2
99001 FORMAT (1A3,1A4)

      IF ( op1.EQ.'OP, ' ) THEN
         IF ( op2.EQ.'GOSI' ) oph = op2
         IF ( op2.EQ.'GOSI' ) opcja = op2

C        Treat OP,FILE (attach files to fortran units)
         IF ( op2.EQ.'FILE' ) THEN
            CALL OPENF
            GOTO 100 ! Back to input loop
         ENDIF
         IF ( jphd.EQ.1 ) WRITE (22,99101)
99101    FORMAT ('1'/1X,125('*')/1X,125('*')/1X,50('*'),25X,50('*')/1X,
     &           50('*'),10X,'GOSIA2',10X,50('*')/1X,50('*'),25X,50('*')
     &           /1X,125('*')/1X,125('*')////)
         IF ( jphd.EQ.1 ) WRITE (22,99002)
99002    FORMAT (1X/20X,'ROCHESTER COULOMB EXCITATION DATA ANALYSIS ',
     &           'CODE BY T.CZOSNYKA,D.CLINE AND C.Y.WU'/50X,
     &           'LATEST REVISION- JUNE  2006'//////)
         jphd = 0

C        Handle OP,GDET (germanium detectors)
         IF ( op2.EQ.'GDET' ) THEN
            nl = 7
            READ (JZB,*) nfdd ! number of physical detectors

            nfd = ABS(nfdd) ! Negative value means graded absorber
            IF ( nfdd.LE.0 ) THEN
               REWIND 8
               DO i = 1 , nl
                  WRITE (8,*) (tau2(l,i),l=1,10)
               ENDDO
               WRITE (8,*) (eng(l),l=1,10)
            ENDIF

C           Write file for gamma-ray energy dependence of Ge solid-angle
C           attenuation coefficients
            REWIND 9
            WRITE (9,*) nfd
            DO i = 1 , nfd ! For each detector
               READ (JZB,*) (DIX(k),k=1,4) ! radius of core, outer radius, length, distance
               READ (JZB,*) (xl1(k),k=1,nl) ! thicknesses of 7 kinds of absorber
               IF ( DIX(1).LE.0. ) DIX(1) = .01
               WRITE (9,*) DIX(4)
               IF ( nfdd.LE.0 ) WRITE (8,*) (xl1(k),k=1,nl)
               ind = 1
               IF ( xl1(5).GT.0. ) ind = 3
               IF ( xl1(6).GT.0. ) ind = 4
               IF ( xl1(7).GT.0. ) ind = 5
               WRITE (9,*) eng(ind)
               CALL QFIT(qui,tau1,tau2,eng,xl1,cf,nl,ind)
               WRITE (22,99003) i
99003          FORMAT (10X,'DETECTOR',1X,1I2)
               DO k = 1 , 8
                  WRITE (22,99004) k , cf(k,1) , cf(k,2)
99004             FORMAT (1X,//5X,'K=',1I1,2X,'C1=',1E14.6,2X,'C2=',
     &                    1E14.6/5X,'ENERGY(MEV)',5X,'FITTED QK',5X,
     &                    'CALC.QK',5X,'PC.DIFF.'/)
                  WRITE (9,*) cf(k,1) , cf(k,2) , qui(k,ind)
                  DO l = 1 , 10
                     arg = (eng(l)-eng(ind))**2
                     qc = (qui(k,ind)*cf(k,2)+cf(k,1)*arg)/(cf(k,2)+arg)
                     WRITE (22,99005) eng(l) , qc , qui(k,l) , 
     &                                100.*(qc-qui(k,l))/qui(k,l)
99005                FORMAT (8X,1F4.2,6X,1F9.4,5X,1F9.4,3X,1E10.2)
                  ENDDO
               ENDDO
            ENDDO
            GOTO 100 ! Back to input loop

C         Treat OP,RAND (randomise matrix elements)
         ELSEIF ( op2.EQ.'RAND' ) THEN
            READ (JZB,*) SE ! Seed for random number generator
            CALL MIXUP
            WRITE (22,99006)
99006       FORMAT (1X///5X,'MATRIX ELEMENTS RANDOMIZED...'///)
            CALL PRELM(2)
            GOTO 100 ! End of OP,RAND

C        Treat OP,TROU (troubleshooting)
         ELSEIF ( op2.EQ.'TROU' ) THEN
            ITS = 1
            READ (JZB,*) kmat , rlr
            GOTO 100 ! Back to input loop

C        Treat OP,REST (restart)
         ELSEIF ( op2.EQ.'REST' ) THEN
            irix = 12
            IF ( IBPS.EQ.1 ) irix = 32 ! unit 12 for target, 32 for beam
            REWIND irix
            memax1 = MEMAX + 1
            DO lkj = 1 , MEMAX
               READ (irix,*) ELM(lkj)
            ENDDO
            DO lkj = 1 , memax1
               READ (JZB,*) lkj1 , xlk
               IF ( lkj1.EQ.0 ) GOTO 220
               IF ( mres1.EQ.0 ) ELM(lkj1) = xlk ! condition added for gosia2
            ENDDO
C---- gosia2 changes start
 220        DO lkj = 1 , MEMAX
               IF ( mres1.EQ.0 .AND. JZB.EQ.25 ) ELM25(lkj) = ELM(lkj)
               IF ( mres1.EQ.0 .AND. JZB.EQ.26 ) ELM26(lkj) = ELM(lkj)
            ENDDO
            IF ( JZB.EQ.26 ) mres1 = 1
            WRITE (22,99007)
C---- gosia2 changes end
99007       FORMAT (1X///5X,'*****',2X,
     &              'RESTART-MATRIX ELEMENTS OVERWRITTEN',2X,'*****'///)
            DO kk = 1 , MEMAX
               la = mlt(kk)
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
            GOTO 100 ! End of OP,REST

C        Treat other options
         ELSE

C           Treat OP,RE,A (release A)
            IF ( op2.EQ.'RE,A' ) GOTO 1000
           
C           Treat OP,RE,F (release F)
            IF ( op2.EQ.'RE,F' ) GOTO 1000

C           Treat OP,ERRO (calculate errors)
            IF ( op2.EQ.'ERRO' ) THEN
               READ (JZB,*) idf , ms , mend , irep , ifc , remax
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
               JENTR = 1
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
                  ENDIF
               ENDIF
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
               IF ( ERR ) GOTO 2000
               IF ( IMIN.EQ.0 ) GOTO 2400
               GOTO 500 ! End of OP,ERRO

C           Treat OP,RE,C (release C)
            ELSEIF ( op2.EQ.'RE,C' ) THEN
               jfre = 1
               irfix = 0
               GOTO 1100 ! End of OP,RE,C

C           Treat OP,TITL (title)
            ELSEIF ( op2.EQ.'TITL' ) THEN
               READ (JZB,99008) (title(k),k=1,20)
99008          FORMAT (20A4)
               WRITE (22,99009) (title(k),k=1,20)
99009          FORMAT (10X,20A4/10X,100('-'))
               GOTO 100 ! End of OP,TITL

            ELSE
C              Treat OP,GOSI
               IF ( op2.EQ.'GOSI' ) GOTO 300

C              Treat OP,COUL
               IF ( op2.EQ.'COUL' ) GOTO 300

C              Treat OP,EXIT
               IF ( op2.EQ.'EXIT' ) GOTO 3000 ! End of OP,EXIT

C              Treat OP,MINI
               IF ( op2.EQ.'MINI' ) THEN
                  READ (JZB,*) imode , nptl , chiok , conu , xtest , 
     &                         LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
                  op2 = opcja
                  IMIN = IMIN + 1
                  IF ( IMIN.EQ.1 ) GOTO 1300
                  GOTO 1400 ! End of OP,MINI

C              Treat OP,THEO
               ELSEIF ( op2.EQ.'THEO' ) THEN
                  REWIND (12)
                  ibaf = 1
                  DO jb = 1 , LP1 ! LP1 = 50
                     DO lb = 1 , 2
                        xlevb(jb,lb) = 0
                     ENDDO
                  ENDDO
                  READ (JZB,*) nbands ! Number of bands
                  IF ( nbands.LE.0 ) ibaf = 0
                  nbands = ABS(nbands)
                  DO nl = 1 , 8
                     DO jb = 1 , nbands
                        DO jl = 1 , nbands
                           DO kl = 1 , 3
                              bm(nl,jb,jl,kl) = 0.
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
                  DO jb = 1 , nbands
                     READ (JZB,*) bk , ilevls ! K of band, number of levels in band
                     READ (JZB,*) (levl(ib),ib=1,ilevls) ! Level list for band
                     DO kb = 1 , ilevls
                        inva = levl(kb)
                        xlevb(inva,2) = bk
                        xlevb(inva,1) = DBLE(jb)
                     ENDDO
                  ENDDO
                  DO nl = 1 , 8
                     READ (JZB,*) nnl ! Multipolarity
 222                 IF ( nnl.LE.0 ) GOTO 225
                     READ (JZB,*) jb1 , jb2 ! band indices
                     IF ( jb1.NE.0 ) THEN
                        READ (JZB,*) (bm(nnl,jb1,jb2,j),j=1,3) ! intrinsic moments
                        DO j = 1 , 3
                           bm(nnl,jb2,jb1,j) = bm(nnl,jb1,jb2,j)
                        ENDDO
                        GOTO 222
                     ENDIF
                  ENDDO
 225              DO kb = 1 , MEMAX
                     IF ( ibaf.NE.0 ) THEN
                        ind1 = LEAD(1,kb)
                        ind2 = LEAD(2,kb)
                        xi1 = SPIN(ind1)
                        xi2 = SPIN(ind2)
                        lamd = mlt(kb)
                        nb1 = INT(xlevb(ind1,1)+.1)
                        nb2 = INT(xlevb(ind2,1)+.1)
                        xk1 = xlevb(ind1,2)
                        xk2 = xlevb(ind2,2)
                        xm1 = bm(lamd,nb1,nb2,1)
                        xm2 = bm(lamd,nb1,nb2,2)
                        xm3 = bm(lamd,nb1,nb2,3)
                        ELM(kb) = ELMT(xi1,xi2,lamd,nb1,nb2,xk1,xk2,xm1,
     &                            xm2,xm3)
                        IF ( ABS(ELM(kb)).LT.1E-6 ) ELM(kb) = 1.E-6
                        irix = 12
                        IF ( IBPS.EQ.1 ) irix = 32 ! unit 12 for target, 32 for beam
                        WRITE (irix,*) ELM(kb)
                     ENDIF
                  ENDDO
                  GOTO 100 ! End of OP,THEO

C              Treat OP,YIEL
               ELSEIF ( op2.EQ.'YIEL' ) THEN
                  CALL ADHOC(oph,idr,nfd,ntap,iyr)
                  GOTO 100 ! End of OP,YIEL

C              Treat OP,INTG
               ELSEIF ( op2.EQ.'INTG' ) THEN
                  REWIND 14
                  lfagg = 1
                  IF ( SPIN(1).LT..25 ) ISO = 0
                  DO lx = 1 , NEXPT
                     lpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) lpin = jpin(lx)
                     ENDIF
                     IEXP = lx
                     tth = TLBDG(lx)
                     enh = EP(lx)
                     DO mpin = 1 , lpin
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           READ (JZB,*) ne , ntt , emn , emx , wth , 
     &                                  wph , wthh
                           mfla = 1
                           CALL COORD(wth,wph,wthh,ntt,0,pfi,wpi,tth,lx,
     &                                tmn,tmx)
                        ELSE
                           READ (JZB,*) ne , ntt , emn , emx , tmn , tmx
                           mfla = 0
                           IF ( ntt.LT.0 ) mfla = 1
                        ENDIF
                        ntt = ABS(ntt)
                        jan = NANG(lx)
                        jan1 = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) jan1 = jan
                        IF ( iecd(lx).EQ.1 ) THEN ! Circular detector
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , wth , wph , wthh
                        ELSE
                           WRITE (14,*) ne , ntt , emn , emx , tmn , 
     &                                  tmx , jan1 , tmx , tmx , tmx
                        ENDIF
                        READ (JZB,*) (XV(i),i=1,ne)
                        IF ( iecd(lx).NE.1 ) READ (JZB,*)
     &                       (YV(i),i=1,ntt)
                        IF ( tth.LT.0. ) ELMH(2*lx-1) = YV(1)
                        IF ( tth.LT.0. ) ELMH(2*lx) = YV(ntt)
                        DO kloop = 1 , ne
                           enb = XV(kloop)
                           EP(lx) = enb
                           DO ktt = 1 , ntt
                              tta = SIGN(YV(ktt),tth)
                              IF ( IAXS(lx).NE.0 ) THEN
                                 IF ( iecd(lx).NE.1 ) THEN
                                    IF ( kloop.EQ.1 ) THEN
                                       READ (JZB,*) nfi
                                       READ (JZB,*)
     &                                    (fiex1(ktt,jfi,1),fiex1(ktt,
     &                                    jfi,2),jfi=1,nfi)
                                       IF ( tth.LT.0. ) THEN
                                         DO jfi = 1 , nfi
                                         fiex1(ktt,jfi,1)
     &                                      = fiex1(ktt,jfi,1) + 180.
                                         fiex1(ktt,jfi,2)
     &                                      = fiex1(ktt,jfi,2) + 180.
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                 ENDIF
                              ENDIF
                              TLBDG(lx) = tta
                              IF ( kloop.EQ.1 ) THEN
                                 IF ( iecd(lx).NE.0 ) THEN
                                    nfi = 1
                                    fiex1(ktt,1,1) = wpi(ktt,1)
                                    fiex1(ktt,1,2) = wpi(ktt,2)
                                 ENDIF
                              ENDIF
                              CALL CMLAB(lx,dsig,tetrc)
                              IF ( ERR ) GOTO 2000
                              tting = TLBDG(lx)
                              IF ( ERR ) GOTO 1900
                              CALL LOAD(lx,1,1,0.D0,jj)
                              CALL ALLOC(ACCUR)
                              CALL SNAKE(lx,ZPOL)
                              CALL SETIN
                              DO j = 1 , LMAX
                                 polm = DBLE(j-1) - SPIN(1)
                                 CALL LOAD(lx,2,1,polm,jj)
                                 CALL STING(jj)
                                 CALL PATH(jj)
                                 CALL INTG(IEXP)
                                 CALL TENB(j,bten,LMAX)
                              ENDDO
                              CALL TENS(bten)
                              CALL DECAY(ccd,0,ccc)
                              DO j = 1 , LP2
                                 DO ijan = 1 , 20
                                    SUMCL(ijan,j) = 0.
                                 ENDDO
                              ENDDO
                              ija0 = 0
                              DO ijan = 1 , jan
                                 IF ( IAXS(lx).EQ.0 ) nfi = 1
                                 DO jyi = 1 , idr
                                    GRAD(jyi) = 0.
                                 ENDDO
                                 todfi = 0.
                                 DO jfi = 1 , nfi
                                    fi0 = fiex1(ktt,jfi,1)/57.2957795
                                    fi1 = fiex1(ktt,jfi,2)/57.2957795
                                    gth = AGELI(IEXP,ijan,1)
                                    fm = (fi0+fi1)/2.
                                    figl = AGELI(IEXP,ijan,2)
                                    CALL ANGULA(YGN,idr,1,fi0,fi1,tetrc,
     &                                 gth,figl,ijan,op2)
                                    IF ( IFMO.NE.0 ) THEN
                                       id = ITMA(IEXP,ijan)
                                       d = ODL(id)
                                       rx = d*SIN(gth)*COS(figl-fm)
     &                                    - .25*SIN(tetrc)*COS(fm)
                                       ry = d*SIN(gth)*SIN(figl-fm)
     &                                    - .25*SIN(tetrc)*SIN(fm)
                                       rz = d*COS(gth) - .25*COS(tetrc)
                                       rl = SQRT(rx*rx+ry*ry+rz*rz)
                                       sf = d*d/rl/rl
                                       thc = TACOS(rz/rl)
                                       fic = ATAN2(ry,rx)
                                       CALL ANGULA(YGP,idr,1,fi0,fi1,
     &                                    tetrc,thc,fic,ijan,op2)
                                       DO ixl = 1 , idr
                                         ixm = KSEQ(ixl,3)
                                         tfac = TAU(ixm)
                                         YGN(ixl) = YGN(ixl)
     &                                      + .01199182*tfac*BETAR(IEXP)
     &                                      *(sf*YGP(ixl)-YGN(ixl))
                                       ENDDO
                                    ENDIF
                                    IF ( IRAWEX(lx).NE.0 ) THEN
                                       ipd = ITMA(lx,ijan)
                                       DO jyi = 1 , idr
                                         ni = KSEQ(jyi,3)
                                         nf = KSEQ(jyi,4)
                                         decen = EN(ni) - EN(nf)
                                         cocos = SIN(tetrc)*SIN(gth)
     &                                      *COS(fm-figl) + COS(tetrc)
     &                                      *COS(gth)
                                         decen = decen*(1.+BETAR(lx)
     &                                      *cocos)
                                         CALL EFFIX(ipd,decen,effi)
                                         YGN(jyi) = YGN(jyi)*effi
                                       ENDDO
                                       inclus = ICLUST(lx,ijan)
                                       IF ( inclus.NE.0 ) THEN
                                         DO jyi = 1 , idr
                                         SUMCL(inclus,jyi)
     &                                      = SUMCL(inclus,jyi)
     &                                      + YGN(jyi)
                                         ENDDO
                                         IF ( ijan.NE.LASTCL(lx,inclus)
     &                                      ) GOTO 226
                                         DO jyi = 1 , idr
                                         YGN(jyi) = SUMCL(inclus,jyi)
                                         ENDDO
                                       ENDIF
                                    ENDIF
                                    IF ( jfi.EQ.1 ) ija0 = ija0 + 1
                                    DO jyi = 1 , idr
                                       GRAD(jyi) = GRAD(jyi) + YGN(jyi)
                                    ENDDO
                                    todfi = todfi + ABS(fi1-fi0)
                                 ENDDO
                                 IF ( IAXS(lx).EQ.0 ) todfi = 6.283185
                                 ax = 1.
                                 IF ( mfla.EQ.1 ) ax = 1./todfi
                                 dsx = dsig
                                 IF ( mfla.NE.1 ) dsx = dsig*todfi
                                 dsxm(mpin,kloop,ktt) = dsx
                                 WRITE (17,*) lx , mpin , kloop , ktt , 
     &                                  dsx
                                 WRITE (14,*) lx , enb , tting , ija0 , 
     &                                  dsx , 
     &                                  (GRAD(jyi)*dsig*ax,jyi=1,idr)
                                 IF ( IPRM(11).EQ.1 ) THEN
                                    WRITE (22,99048) lx , ija0 , enb , 
     &                                 tta
                                    IF ( tta.LT.0. ) WRITE (22,99010)
     &                                 tting
99010                               FORMAT (5X,
     &                             'RESPECTIVE TARGET SCATTERING ANGLE='
     &                             ,1F7.3,1X,'DEG'/)
                                    DO jyi = 1 , idr
                                       ni = KSEQ(jyi,3)
                                       nf = KSEQ(jyi,4)
                                       WRITE (22,99049) ni , nf , 
     &                                    SPIN(ni) , SPIN(nf) , 
     &                                    GRAD(jyi)*dsig*ax , GRAD(jyi)
     &                                    /GRAD(IDRN)
                                    ENDDO
                                 ENDIF
 226                          ENDDO
                           ENDDO
                        ENDDO
                     ENDDO
                     EP(lx) = enh
                     TLBDG(lx) = tth
                  ENDDO
                  REWIND 14
                  REWIND 15
                  iske = 0
                  DO na = 1 , LP6
                     ILE(na) = 1
                  ENDDO
                  ilx = 0
                  DO lx = 1 , NEXPT
                     REWIND 17
                     DO ijaja = 1 , 300000
                        READ (17,*,END=228) jjlx , jmpin , jkloo , 
     &                        jktt , dsx
                        IF ( jjlx.EQ.lx ) dsxm(jmpin,jkloo,jktt) = dsx
                     ENDDO
 228                 na = NANG(lx)
                     IF ( lx.NE.1 ) THEN
                        DO na1 = 1 , LP6
                           ILE(na1) = ILE(na1) + NYLDE(lx-1,na1)
                        ENDDO
                     ENDIF
                     READ (JZB,*) nptx
                     IF ( nptx.NE.0 ) THEN
                        READ (JZB,*) (esp(i),i=1,nptx)
                        READ (JZB,*) (dedx(i),i=1,nptx)
                        npt = nptx
                     ENDIF
                     READ (JZB,*) npce , npct
                     mfla = 0
                     IF ( npct.LT.0 ) mfla = 1
                     IF ( iecd(lx).EQ.1 ) mfla = 1
                     npct = ABS(npct)
                     npce = npce + MOD(npce,2)
                     npct = npct + MOD(npct,2)
                     mpin = 1
                     IF ( ipinf.NE.0 ) THEN
                        IF ( jpin(lx).NE.0 ) mpin = jpin(lx)
                     ENDIF
                     dst = 0.
                     DO lpin = 1 , mpin
                        ilx = ilx + 1
                        IF ( ilx.NE.1 )
     &                       CALL TAPMA(lx,iske,isko,iskf,nflr,idr,0,
     &                       nft,enb)
                        READ (14,*) ne , ntt , emn , emx , tmn , tmx , 
     &                              jan , wth , wph , wthh
                        iocc = (ne+ntt)*idr
                        IF ( iocc.GT.izcap ) GOTO 1800
                        hen = (emx-emn)/npce
                        npce1 = npce + 1
                        het = (tmx-tmn)/npct
                        npct1 = npct + 1
                        IF ( iecd(lx).EQ.1 ) ! Circular detector
     &                       CALL COORD(wth,wph,wthh,npct1,1,pfi,wpi,
     &                       TLBDG(lx),lx,tmn,tmx)
                        IF ( iecd(lx).NE.1 ) THEN
                           IF ( mfla.EQ.1 ) READ (JZB,*)
     &                          (pfi(j),j=1,npct1)
                        ENDIF
                        het = het/57.2957795
                        DO j = 1 , npce1
                           xx = (j-1)*hen + emn
                           CALL LAGRAN(esp,dedx,npt,1,xx,yy,3,1)
                           HLMLM(j) = 1./yy
                        ENDDO
                        naa = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) naa = NANG(lx)
                        iskf = naa - 1
                        DO ja = 1 , naa
                           icll = 3
                           DO je = 1 , ne
                              lu = ILE(ja)
                              isko = (je-1)*naa*ntt + ja - 1
                              CALL TAPMA(lx,iske,isko,iskf,ntt,idr,1,
     &                           nft,enb)
                              IF ( nft.EQ.1 ) GOTO 1900
                              DO jd = 1 , idr
                                 DO jtp = 1 , ntt
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                                 DSG(jtp) = dsxm(lpin,je,jtp)
                                    jyv = (jtp-1)*idr + jd
                                    YV(jtp) = ZETA(jyv)
                                 ENDDO
                                 DO jt = 1 , npct1
                                    xx = (jt-1)*het + tmn/57.2957795
                                    CALL LAGRAN(XV,YV,ntt,jt,xx,yy,2,
     &                                 icll)
                                    CALL LAGRAN(XV,DSG,ntt,jt,xx,zz,2,
     &                                 icll)
                                    IF ( mfla.EQ.1 ) yy = yy*pfi(jt)
     &                                 /57.2957795
                                    IF ( yy.LE.0. ) yy = 1.E-15
                                    IF ( mfla.EQ.1 ) zz = zz*pfi(jt)
     &                                 /57.2957795
                                    XI(jt) = yy*SIN(xx)
                                    IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                                 = zz*SIN(xx)
                                 ENDDO
                                 icll = 4
                                 locat = ntt*idr + (je-1)*idr + jd
                                 ZETA(locat) = SIMIN(npct1,het,XI)
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) DSE(je)
     &                                = SIMIN(npct1,het,HLM)
                                 ZV(je) = enb
                              ENDDO
                           ENDDO
                           icll = 3
                           DO jd = 1 , idr
                              DO jtp = 1 , ne
                                 jyv = (jtp-1)*idr + jd + ntt*idr
                                 YV(jtp) = ZETA(jyv)
                              ENDDO
                              DO jt = 1 , npce1
                                 xx = (jt-1)*hen + emn
                                 CALL LAGRAN(ZV,YV,ne,jt,xx,yy,2,icll)
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                                CALL LAGRAN(ZV,DSE,ne,jt,xx,zz,2,
     &                                icll)
                                 IF ( jd.EQ.1 .AND. ja.EQ.1 ) HLM(jt)
     &                                = zz*HLMLM(jt)
                                 XI(jt) = yy*HLMLM(jt)
                              ENDDO
                              icll = 4
                              IF ( jd.EQ.1 .AND. ja.EQ.1 )
     &                             DS = SIMIN(npce1,hen,HLM)
                              GRAD(jd) = SIMIN(npce1,hen,XI)
                           ENDDO
                           IF ( ja.EQ.1 ) dst = dst + DS
                           IF ( ja.EQ.1 ) WRITE (22,99111) DS , lx
99111                      FORMAT (1X/////5X,
     &                            'INTEGRATED RUTHERFORD CROSS SECTION='
     &                            ,1E9.4,2X,'FOR EXP.',1I2///)
                           WRITE (22,99112) lx , ja , emn , emx , tmn , 
     &                            tmx
99112                      FORMAT (1X,//50X,'INTEGRATED YIELDS'//5X,
     &                             'EXPERIMENT ',1I2,2X,'DETECTOR ',
     &                             1I2/5X,'ENERGY RANGE ',1F8.3,'---',
     &                             1F8.3,1X,'MEV',3X,
     &                             'SCATTERING ANGLE RANGE ',1F7.3,
     &                             '---',1F7.3,1X,'DEG'//5X,'NI',5X,
     &                             'NF',5X,'II',5X,'IF',5X,'YIELD',5X,
     &                             'NORMALIZED YIELD'/)
                           DO jd = 1 , idr
                              WRITE (15,*) GRAD(jd)
                           ENDDO
                           DO jd = 1 , idr
                              ni = KSEQ(jd,3)
                              nf = KSEQ(jd,4)
                              WRITE (22,99049) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , GRAD(jd) , GRAD(jd)
     &                               /GRAD(IDRN)
                           ENDDO
                        ENDDO
                        IF ( iecd(lx).EQ.1 ) THEN
                           IF ( jpin(lx).EQ.0 ) THEN
                              CALL COORD(wth,wph,wthh,1,2,pfi,wpi,
     &                           TLBDG(lx),lx,txx,txx)
                              WRITE (22,99113) FIEX(lx,1)*57.2957795 , 
     &                               FIEX(lx,2)*57.2957795 , lx
99113                         FORMAT (//5X,
     &                          'WARNING: THE PHI ANGLE WAS REPLACED BY'
     &                          ,1X,F8.3,1X,'TO',F8.3,3X,
     &                          'FOR EXPERIMENT',2X,I3)
                              IF ( TLBDG(lx).LT.0 ) THEN
                                 FIEX(lx,1) = FIEX(lx,1) + 3.14159265
                                 FIEX(lx,2) = FIEX(lx,2) + 3.14159265
                              ENDIF
                           ENDIF
                        ENDIF
                        iske = iske + ne*ntt*naa
                     ENDDO
                     IF ( mpin.GT.1 ) WRITE (22,99114) dst , lx
99114                FORMAT (1x//2x,
     &                      'Total integrated Rutherford cross section='
     &                      ,1E8.3,' for exp. ',1I2/)
                  ENDDO
                  IF ( ipinf.NE.0 ) THEN
                     ngpr = 0
                     DO lx = 1 , NEXPT
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        IF ( lx.NE.1 ) ngpr = ngpr + idr*jpin(lx-1)
     &                       *NDST(lx-1)
                        lpin = jpin(lx)
                        IF ( lpin.EQ.0 ) lpin = 1
                        DO jgd = 1 , nged
                           DO jd = 1 , idr
                              GRAD(jd) = 0.
                           ENDDO
                           DO mpin = 1 , lpin
                              REWIND 15
                              ndum = ngpr + (jgd-1)*idr + (mpin-1)
     &                               *jgd*idr
                              IF ( ndum.NE.0 ) THEN
                                 DO jd = 1 , ndum
                                    READ (15,*) xx
                                 ENDDO
                              ENDIF
                              DO jd = 1 , idr
                                 READ (15,*) xx
                                 GRAD(jd) = GRAD(jd) + xx
                              ENDDO
                           ENDDO
                           WRITE (17,*) (GRAD(jd),jd=1,idr)
                        ENDDO
                     ENDDO
                     REWIND 15
                     REWIND 17
                     DO lx = 1 , NEXPT
                        nged = NDST(lx)
                        IF ( IRAWEX(lx).EQ.0 ) nged = NANG(lx)
                        DO ija0 = 1 , nged
                           READ (17,*) (GRAD(jdy),jdy=1,idr)
                           DO jd = 1 , idr
                              WRITE (15,*) GRAD(jd)
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDIF
                  GOTO 100
               ELSEIF ( op2.EQ.'CORR' ) THEN
                  CALL READY(idr,ntap,0)
                  REWIND 3
                  REWIND 15
                  REWIND 4
                  GOTO 1300
               ELSE

C                 Treat OP,POIN
                  IF ( op2.EQ.'POIN' ) GOTO 1300

C                 Treat OP,MAP
                  IF ( op2.EQ.'MAP ' ) iobl = 1

C                 Treat OP,STAR
                  IF ( op2.EQ.'STAR' ) GOTO 1300

C                 Treat OP,SIXJ
                  IF ( op2.EQ.'SIXJ' ) THEN
                     DO k = 1 , 2
                        l = 4*k
                        DO j = 1 , 80
                           ixj = j - 1
                           DO ms = 1 , 5
                              mend = 2*(ms-3) + ixj
                              WRITE (14,*) WSIXJ(l,4,4,ixj,mend,ixj-4) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj-2) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj+2) , 
     &                               WSIXJ(l,4,4,ixj,mend,ixj+4)
                           ENDDO
                        ENDDO
                     ENDDO
                     GOTO 2000 ! End of OP,SIXJ

C                 Treat OP,RAW (raw uncorrected gamma yields)
                  ELSEIF ( op2.EQ.'RAW ' ) THEN
C                    Read absorber coefficients from unit 8
                     REWIND 8
                     DO l = 1 , 8
                        READ (8,*) (ABC(l,j),j=1,10) ! Absorption coefficients
                        DO j = 1 , 10
                           ABC(l,j) = LOG(ABC(l,j))
                        ENDDO
                     ENDDO
                     DO l = 1 , nfd
                        READ (8,*) (THICK(l,j),j=1,7) ! thickness of absorbers
                     ENDDO
                     DO l = 1 , LP1 ! LP1 = 50
                        DO j = 1 , 200
                           ICLUST(l,j) = 0
                        ENDDO
                        DO j = 1 , 20
                           LASTCL(l,j) = 0
                        ENDDO
                        IRAWEX(l) = 0
                     ENDDO

C                    Read input from standard input
                     DO l = 1 , LP1 ! LP1 = 50
                        READ (JZB,*) mexl ! experiment number
                        IF ( mexl.EQ.0 ) GOTO 100
                        IRAWEX(mexl) = 1
                        n = NANG(mexl)
                        DO j = 1 , n
                           jj = ITMA(mexl,j)
                           READ (JZB,*) (AKAVKA(k,jj),k=1,8) ! efficiency curve parameters
                        ENDDO
                        READ (JZB,*) kclust
                        IF ( kclust.NE.0 ) THEN
                           DO j = 1 , kclust
                              READ (JZB,*) numcl ! Number of detectors for this cluster
                              READ (JZB,*) (liscl(k),k=1,numcl) ! Indices of logical detectors
                              LASTCL(l,j) = liscl(numcl)
                              DO k = 1 , numcl
                                 kk = liscl(k)
                                 ICLUST(l,kk) = j
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                     GOTO 100 ! End of OP,RAW

C                 Treat OP,MAP
                  ELSEIF ( op2.EQ.'MAP ' ) THEN
                     GOTO 1300 ! End of OP,MAP 
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF ! End of if (op1.eq."OP, ") if statement

      WRITE (22,99115) op1 , op2
99115 FORMAT (5X,'UNRECOGNIZED OPTION',1X,1A3,1A4)
      GOTO 2000

C     Treat suboptions of OP,COUL and OP,GOSI
 300  READ (JZB,99116) op1 ! Read the suboption
99116 FORMAT (1A4)
      IF ( op1.EQ.'    ' ) GOTO 100

C     Treat suboption LEVE (levels)
      IF ( op1.EQ.'LEVE' ) THEN
         NMAX = 0
         IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99017)
99017    FORMAT (1X/40X,'LEVELS',//5X,'INDEX',5X,'PARITY',9X,'SPIN',11X,
     &           'ENERGY(MEV)')
         ndima = NDIM + 1
         DO k = 1 , ndima
            READ (JZB,*) ipo1 , ipo2 , po2 , po1
            IF ( ipo1.EQ.0 ) GOTO 300
            IF ( ipo1.EQ.1 .AND. ABS(po2).LT.1.E-6 ) ISO = 0
            NMAX = NMAX + 1
            SPIN(ipo1) = po2
            IF ( k.EQ.1 ) iph = ipo2
            iprc = ipo2 - iph
            IF ( iprc.NE.0 ) iprc = 1
            IFAC(ipo1) = (-1)**(iprc-INT(po2-SPIN(1)))
            EN(ipo1) = po1
            prp = '+'
            IF ( ipo2.EQ.-1 ) prp = '-'
            IF ( ABS(IPRM(1)).EQ.1 ) WRITE (22,99018) ipo1 , prp , 
     &           SPIN(ipo1) , EN(ipo1)
99018       FORMAT (6X,1I2,11X,1A1,10X,1F4.1,8X,1F10.4)
         ENDDO
      ELSEIF ( op1.EQ.'ME  ' ) THEN
         DO k = 1 , nmemx
            IF ( op2.EQ.'GOSI' ) THEN
               READ (JZB,*) ipo1 , ipo2 , po1 , bl , bu
               iopri = 2
               icg = 2
            ELSE
               iopri = 1
               READ (JZB,*) ipo1 , ipo2 , po1
            ENDIF
            IF ( ipo1.NE.0 ) THEN
               IF ( ipo2.EQ.0 ) THEN
                  IF ( ipo1.LE.la ) GOTO 1600
                  LAMMAX = LAMMAX + 1
                  LAMDA(LAMMAX) = ipo1
                  ipo3 = 0
                  IF ( indx.EQ.0 ) GOTO 320
               ELSE
                  MULTI(la) = MULTI(la) + 1
                  indx = indx + 1
                  IF ( ipo1.GT.ABS(ipo2) ) GOTO 1500
                  IF ( ipo1.NE.ipo3 ) THEN
                     IF ( ipo1.LT.ipo3 ) GOTO 1700
                     ipo3 = ipo1
                  ENDIF
                  ELM(indx) = po1
                  mlt(indx) = la
                  LEAD(1,indx) = ipo1
                  LEAD(2,indx) = ABS(ipo2)
                  LDNUM(la,ipo1) = LDNUM(la,ipo1) + 1
                  IF ( op2.EQ.'GOSI' ) THEN
                     IF ( ipo2.LT.0 ) THEN
                        IVAR(indx) = 10000*INT(bl) + INT(bu)
                     ELSE
                        ELMU(indx) = bu
                        ELML(indx) = bl
                        IF ( ABS(bl-bu).LT.1.E-6 ) THEN
                           IVAR(indx) = 0
                        ELSE
                           IVAR(indx) = 2
                           IF ( la.GT.4 ) IVAR(indx) = 1
                        ENDIF
                     ENDIF
                     isip = ISEX(ipo1) + 1
                     ISEX(ABS(ipo2)) = MIN(isip,ISEX(ABS(ipo2)))
                  ENDIF
                  GOTO 350
               ENDIF
            ENDIF
            DO kk = 1 , indx
               IF ( ABS(ELM(kk)).LE.1.E-6 ) ELM(kk) = 1.E-6
               IF ( IVAR(kk).GE.10000 ) THEN
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
            IF ( ipo1.EQ.0 ) GOTO 400
 320        la = ipo1
            IF ( la.GT.LMAXE .AND. la.LE.6 ) LMAXE = la
 350     ENDDO
 400     MEMAX = indx
         IF ( JZB.EQ.25 ) m25 = MEMAX
         IF ( JZB.EQ.26 ) m26 = MEMAX
         IF ( la.GT.6 ) MAGEXC = 1
         memx4 = MULTI(1) + MULTI(2) + MULTI(3) + MULTI(4)
         MEMX6 = memx4 + MULTI(5) + MULTI(6)
         IF ( ABS(IPRM(1)).EQ.1 ) CALL PRELM(iopri)
         DO kh = 1 , NMAX
            IF ( ISEX(kh).EQ.1111 ) ISEX(kh) = 1
         ENDDO
         DO kh = 1 , MEMAX
            ivarh(kh) = IVAR(kh)
         ENDDO
      ELSEIF ( op1.EQ.'CONT' ) THEN
 450     READ (JZB,99019) op1 , fipo1
99019    FORMAT (1A4,1F7.1)
         ipo1 = INT(fipo1)
         IF ( op1.EQ.'ACP,' ) ACCA = 10.**(-fipo1)
         IF ( op1.EQ.'SEL,' ) ITS = 2
         IF ( op1.EQ.'SMR,' ) iosr = 1
         IF ( op1.EQ.'FMI,' ) ifm = 1
         IF ( op1.EQ.'TEN,' ) itno = 1
         IF ( op1.EQ.'NCM,' ) NCM = ipo1
         IF ( op1.EQ.'WRN,' ) SGW = fipo1
         IF ( op1.EQ.'INT,' ) THEN
            DO jjx = 1 , ipo1
               READ (JZB,*) ipo2 , ijx
               INTERV(ipo2) = ijx
            ENDDO
         ELSE
            IF ( op1.EQ.'VAC,' ) THEN
               DO jjx = 1 , 7
                  READ (JZB,*) ijx , val
                  IF ( ijx.EQ.0 ) GOTO 450
                  G(ijx) = val
               ENDDO
            ELSE
               IF ( op1.EQ.'DIP,' ) DIPOL = 0.001*fipo1
               IF ( op1.EQ.'ACC,' ) ACCUR = 10.**(-fipo1)
               IF ( op1.EQ.'PRT,' ) THEN
                  DO jjx = 1 , 20
                     READ (JZB,*) inm1 , inm2
                     IF ( inm1.EQ.0 ) GOTO 450
                     IPRM(inm1) = inm2
                  ENDDO
                  GOTO 450
               ELSEIF ( op1.NE.'FIX,' ) THEN
                  IF ( op1.EQ.'SKP,' ) THEN
                     DO jjx = 1 , ipo1
                        READ (JZB,*) ijx
                        JSKIP(ijx) = 0
                     ENDDO
                     GOTO 450
                  ELSE
                     IF ( op1.EQ.'CRF,' ) ICS = 1
                     IF ( op1.EQ.'LCK,' ) THEN
 452                    READ (JZB,*) lck1 , lck2
                        IF ( lck1.EQ.0 ) GOTO 450
                        DO jjx = lck1 , lck2
                           ivarh(jjx) = 0
                           IVAR(jjx) = 0
                        ENDDO
                        GOTO 452
                     ELSE
                        IF ( op1.EQ.'INR,' ) INNR = 1
                        IF ( op1.EQ.'CRD,' ) THEN
                           DO jjx = 1 , ipo1
                              READ (JZB,*) ipo2
                              iecd(ipo2) = 1
                           ENDDO
                           GOTO 450
                        ELSE
                           IF ( op1.EQ.'CCF,' ) IPS1 = ipo1
                           IF ( op1.EQ.'PIN,' ) ipine = ipo1
                           IF ( op1.EQ.'PIN,' ) ipinf = 1
                           IF ( op1.EQ.'PIN,' ) THEN
                              DO ipp = 1 , ipine
                                 READ (JZB,*) ig1 , ig2
                                 jpin(ig1) = ig2
                              ENDDO
                              GOTO 450
                           ELSE
                              IF ( op1.EQ.'END,' ) GOTO 300
                              GOTO 450
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
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
         ENDIF
         GOTO 450
      ELSEIF ( op1.EQ.'EXPT' ) THEN
         READ (JZB,*) NEXPT , IZ , XA
         G(1) = 3.
         G(2) = .02
         G(3) = .0345
         G(4) = 3.5
         G(5) = DBLE(IZ)/XA
         G(6) = 6.E-06
         G(7) = .6
         DO k = 1 , NEXPT
            READ (JZB,*) IZ1(k) , XA1(k) , EP(k) , TLBDG(k) , EMMA(k) , 
     &                   MAGA(k) , IAXS(k) , fi0 , fi1 , ISKIN(k) , 
     &                   LNORM(k)
            ITTE(k) = 0
            IF ( XA1(k).LT.0. ) ITTE(k) = 1
            XA1(k) = ABS(XA1(k))
            FIEX(k,1) = fi0/57.2957795
            FIEX(k,2) = fi1/57.2957795
            IF ( TLBDG(k).LT.0. ) THEN
               FIEX(k,1) = FIEX(k,1) + 3.14159265
               FIEX(k,2) = FIEX(k,2) + 3.14159265
            ENDIF
         ENDDO
      ELSE
         WRITE (22,99020) op1
99020    FORMAT (5X,'UNRECOGNIZED SUBOPTION',1X,1A4)
         GOTO 2000
      ENDIF
      GOTO 300
 500  IF ( ICS.EQ.1 ) THEN
         REWIND 11
         DO kh1 = 1 , LP4
            READ (11) (CORF(kh1,kh2),kh2=1,LP6)
         ENDDO
      ELSE
         CALL FTBM(0,chiss,idr,0,chilo,bten)
         REWIND 11
         DO kh1 = 1 , LP4
            WRITE (11) (CORF(kh1,kh2),kh2=1,LP6)
         ENDDO
      ENDIF
      CALL FTBM(3,chiss,idr,1,chilo,bten)
      chis0 = chiss
      WRITE (22,99021) chis0
99021 FORMAT (1X///10X,'***** CENTRAL CHISQ=',1E12.4,1X,'*****'//)
      INHB = 1
      chisl = chiss
      DO kh = 1 , MEMAX
         HLM(kh) = ELM(kh)
      ENDDO
      IF ( idf.EQ.1 ) THEN
         IFBFL = 1
         IF ( irep.NE.2 ) GOTO 800
         IF ( iosr.EQ.0 ) GOTO 800
         REWIND 33
         READ (33,*) ll , mm , kk , inn
         DO inn = 1 , ll
            READ (33,*) mm , yyy , zz
         ENDDO
         DO inn = 1 , MEMAX
            READ (33,*) mm , ll , kk
         ENDDO
         DO inn = 1 , MEMAX
            READ (33,*) mm , yyy
         ENDDO
 550     READ (33,*) mm , ll
         IF ( mm.EQ.0 ) THEN
            BACKSPACE 33
            GOTO 800
         ELSE
            READ (33,*) kk , ll , yyy
            READ (33,*) (SA(mm),mm=1,MEMAX)
            GOTO 550
         ENDIF
      ELSE
         naxfl = 0
         IF ( ms.EQ.0 ) mend = MEMAX
         IF ( ms.EQ.0 ) ms = 1
         DO kh = ms , mend
            DO ij = 1 , 2
               pv = (ELMU(kh)-ELML(kh))/100.
               IF ( ij.NE.1 .OR. (ELM(kh)-ELML(kh)).GE.pv ) THEN
                  IF ( ij.NE.2 .OR. (ELMU(kh)-ELM(kh)).GE.pv ) THEN
                     DO kh1 = 1 , MEMAX
                        SA(kh1) = 0.
                     ENDDO
                     IF ( IVAR(kh).EQ.0 ) GOTO 600
                     SA(kh) = 1.*(-1)**ij
                     kh1 = kh
                     CALL KONTUR(idr,chis0,chisl,ifbp,-1,kh1,sh,bten,
     &                           rem)
                     ELM(kh) = HLM(kh)
                  ENDIF
               ENDIF
            ENDDO
            REWIND 15
            WRITE (15,*) (DEVD(ij),DEVU(ij),ij=1,MEMAX)
 600     ENDDO
      ENDIF
 700  IF ( ifbp.EQ.1 ) THEN
         REWIND 17
         DO lkj = 1 , MEMAX
            READ (17,*) ELM(lkj)
         ENDDO
         WRITE (22,99022)
99022    FORMAT (1X///20X,'*** BEST POINT FOUND (TAPE17) ***'///)
         CALL PRELM(3)
      ENDIF
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      WRITE (22,99023)
99023 FORMAT (40X,'ESTIMATED ERRORS'//5X,'INDEX',5X,'NI',5X,'NF',5X,
     &        'ME AND ERRORS'//)
      DO kh1 = 1 , MEMAX
         IF ( IVAR(kh1).NE.0 .AND. IVAR(kh1).LE.999 ) THEN
            WRITE (22,99024) kh1 , LEAD(1,kh1) , LEAD(2,kh1) , HLM(kh1)
     &                       , DEVD(kh1) , DEVU(kh1) , DEVD(kh1)
     &                       *100./ABS(HLM(kh1)) , DEVU(kh1)
     &                       *100./ABS(HLM(kh1))
99024       FORMAT (6X,1I3,6X,1I2,5X,1I2,5X,1F9.5,2X,'(',1F9.5,' ,',
     &              1F9.5,')','......',1F7.1,' ,',1F7.1,1X,'PC')
         ENDIF
      ENDDO
      IF ( naxfl.NE.0 ) WRITE (22,99050)
      IF ( naxfl.EQ.0 ) WRITE (22,99051)
      WRITE (22,99025)
99025 FORMAT (40X,'ESTIMATED ERRORS',//5X,'INDEX',5X,'NI',5X,'NF',5X,
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
               IF ( ABS(be2a-be2c).LT.1.E-6 ) be2a = be2c
               IF ( be2a/HLM(kh2).LE.0. .OR. be2b/HLM(kh2).LE.0. )
     &              be2a = 0.
               be2a = be2a**2/sbe
               be2b = be2b**2/sbe
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , be2 , 
     &                          be2a - be2 , be2b - be2
            ELSE
               ispb = INT(SPIN(ispa))*2
               qfac = 3.170662*WTHREJ(ispb,4,ispb,-ispb,0,ispb)
               WRITE (22,99052) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , 
     &                          HLM(kh2)*qfac , DEVD(kh2)*qfac , 
     &                          DEVU(kh2)*qfac
            ENDIF
         ENDIF
      ENDDO
      GOTO 2000
 800  irea = 0
      IF ( ms.LT.0 ) irea = 1
      IF ( ms.EQ.0 ) mend = MEMAX
      IF ( ms.EQ.0 ) ms = 1
 900  naxfl = 1
      IF ( irea.EQ.1 ) READ (JZB,*) ms , mend
      IF ( ms.NE.0 ) THEN
         DO kh = ms , mend
            IF ( ifc.NE.1 ) THEN
               REWIND 18
               DO kh1 = 1 , kh
                  READ (18,*) (KVAR(jyi),jyi=1,MEMAX)
               ENDDO
               DO kh1 = 1 , MEMAX
                  ivrh = IVAR(kh1)
                  IF ( KVAR(kh1).EQ.0 ) IVAR(kh1) = 0
                  KVAR(kh1) = ivrh
               ENDDO
            ENDIF
            DO ij = 1 , 2
               sh = DEVU(kh)
               IF ( ij.EQ.1 ) sh = DEVD(kh)
               IF ( ABS(sh).LT.1.E-6 ) sh = (-1)**ij*ABS(HLM(kh))/10.
               ELM(kh) = HLM(kh) + 1.5*sh
               mm = 0
               DO kh1 = 1 , MEMAX
                  IF ( ifc.EQ.1 ) KVAR(kh1) = IVAR(kh1)
                  mm = mm + IVAR(kh1)
               ENDDO
               IF ( mm.EQ.0 ) WRITE (22,99026) kh
99026          FORMAT (10X,'ME=',1I3,5X,'NO FREE MATRIX ELEMENTS')
               IF ( mm.NE.0 ) THEN
                  KFERR = 1
                  IF ( iosr.EQ.1 ) WRITE (33,*) kh , kh
                  IF ( iosr.EQ.1 ) WRITE (33,*) kh , ij , ELM(kh)
                  LOCKS = 1
                  DLOCK = .05
                  CALL MINI(chiss,-1.D0,2,.0001D0,1000,idr,100000.D0,0,
     &                      iosr,kh,bten)
                  DO kh1 = 1 , MEMAX
                     SA(kh1) = (ELM(kh1)-HLM(kh1))/ABS(sh)
                  ENDDO
                  CALL KONTUR(idr,chis0,chisl,ifbp,inpo,kh,sh,bten,rem)
               ENDIF
               DO kh1 = 1 , MEMAX
                  IF ( ifc.EQ.1 ) IVAR(kh1) = KVAR(kh1)
                  ELM(kh1) = HLM(kh1)
               ENDDO
            ENDDO
            IF ( ifc.NE.1 ) THEN
               DO kh1 = 1 , MEMAX
                  IVAR(kh1) = KVAR(kh1)
               ENDDO
            ENDIF
            REWIND 15
            WRITE (15,*) (DEVD(kh1),DEVU(kh1),kh1=1,MEMAX)
         ENDDO
         IF ( irea.EQ.1 ) GOTO 900
      ENDIF
      IF ( iosr.NE.0 ) THEN
         im = 0
         WRITE (33,*) im , im
      ENDIF
      GOTO 700

 1000 jfre = 0
      irfix = 0
      IF ( op2.EQ.'RE,F' ) irfix = 1
 1100 DO jrls = 1 , MEMAX
         IF ( IVAR(jrls).NE.0 .OR. irfix.NE.1 ) THEN
            IF ( IVAR(jrls).GT.999 ) THEN
               IF ( jfre.EQ.1 ) GOTO 1200
            ENDIF
            IVAR(jrls) = 2
            ELML(jrls) = -ABS(ELML(jrls))
            ELMU(jrls) = ABS(ELMU(jrls))
            IF ( jrls.GT.MEMX6 ) IVAR(jrls) = 1
         ENDIF
 1200 ENDDO
      DO jrls = 1 , MEMAX
         ivarh(jrls) = IVAR(jrls)
      ENDDO
      GOTO 100

 1300 CALL CMLAB(0,dsig,ttttt)
      IF ( ERR ) GOTO 2000
      IF ( op2.EQ.'POIN' ) READ (JZB,*) ifwd , slim
      ient = 1
      icg = 1
      IF ( SPIN(1).LT.1.E-6 ) ISO = 0
      IF ( iobl.LT.1 ) THEN
         IF ( op2.NE.'GOSI' ) THEN
            iapx = 0
            DO ii = 1 , LP6
               ILE(ii) = 1
            ENDDO
            nch = 0
            DO jexp = 1 , NEXPT
               IEXP = jexp
               ttttt = TREP(IEXP)
               dsig = DSIGS(IEXP)
               IF ( op2.NE.'STAR' ) THEN
                  jmm = IEXP
                  IF ( IEXP.NE.1 ) THEN
                     DO lli = 1 , LP6
                        ILE(lli) = ILE(lli) + NYLDE(IEXP-1,lli)
                     ENDDO
                  ENDIF
               ENDIF
               fi0 = FIEX(IEXP,1)
               fi1 = FIEX(IEXP,2)
               CALL LOAD(IEXP,1,icg,0.D0,jj)
               CALL ALLOC(ACCUR)
               CALL SNAKE(IEXP,ZPOL)
               CALL SETIN
               DO j = 1 , LMAX
                  polm = DBLE(j-1) - SPIN(1)
                  CALL LOAD(IEXP,2,icg,polm,jj)
                  CALL STING(jj)
                  CALL PATH(jj)
                  CALL INTG(IEXP)
                  CALL TENB(j,bten,LMAX)
                  pr = 0.
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99027) (DBLE(j)-1.-SPIN(1)) , IEXP
99027             FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',
     &                    1F5.1,5X,'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,
     &                    'SPIN',2X,'M',5X,'REAL AMPLITUDE',2X,
     &                    'IMAGINARY AMPLITUDE'//)
                  DO k = 1 , ISMAX
                     pr = pr + DBLE(ARM(k,5))**2 + IMAG(ARM(k,5))**2
                     IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                    WRITE (22,99028) INT(CAT(k,1)) , CAT(k,2) , 
     &                    CAT(k,3) , DBLE(ARM(k,5)) , IMAG(ARM(k,5))
99028                FORMAT (7X,1I2,3X,1F4.1,2X,1F5.1,2X,1E14.6,2X,
     &                       1E14.6)
                  ENDDO
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99029) pr
99029             FORMAT (1X/5X,'SUM OF PROBABILITIES=',1E14.6)
               ENDDO
               CALL TENS(bten)
               IF ( itno.NE.0 ) THEN
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
                  nogeli = NANG(IEXP)
                  jgl1 = 0
                  DO js = 1 , LP2
                     DO jgl = 1 , 20
                        SUMCL(jgl,js) = 0.
                     ENDDO
                  ENDDO
                  DO jgl = 1 , nogeli
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                       WRITE (23,99030) IEXP , jgl , EP(IEXP) , 
     &                       TLBDG(IEXP)
99030                   FORMAT (1x//50x,'CALCULATED YIELDS'//5x,
     &                          'EXPERIMENT ',1I2,2x,'DETECTOR ',1I2/5x,
     &                          'ENERGY ',1F10.3,1x,'MEV',2x,'THETA ',
     &                          1F7.3,1x,'DEG'//5x,'NI',5x,'NF',5x,'II',
     &                          5x,'IF',5x,'E(MeV)',5x,'EFFICIENCY'/)
                     ENDIF
                     gth = AGELI(IEXP,jgl,1)
                     figl = AGELI(IEXP,jgl,2)
                     fm = (fi0+fi1)/2.
                     CALL ANGULA(YGN,idr,1,fi0,fi1,ttttt,gth,figl,jgl,
     &                           op2)
                     IF ( IFMO.NE.0 ) THEN
                        id = ITMA(IEXP,jgl)
                        d = ODL(id)
                        rx = d*SIN(gth)*COS(figl-fm) - .25*SIN(ttttt)
     &                       *COS(fm)
                        ry = d*SIN(gth)*SIN(figl-fm) - .25*SIN(ttttt)
     &                       *SIN(fm)
                        rz = d*COS(gth) - .25*COS(ttttt)
                        rl = SQRT(rx*rx+ry*ry+rz*rz)
                        thc = TACOS(rz/rl)
                        sf = d*d/rl/rl
                        fic = ATAN2(ry,rx)
                        CALL ANGULA(YGP,idr,1,fi0,fi1,ttttt,thc,fic,jgl,
     &                              op2)
                        DO ixl = 1 , idr
                           ixm = KSEQ(ixl,3)
                           tfac = TAU(ixm)
                           YGN(ixl) = YGN(ixl)
     &                                + .01199182*tfac*BETAR(IEXP)
     &                                *(sf*YGP(ixl)-YGN(ixl))
                        ENDDO
                     ENDIF
                     IF ( IRAWEX(IEXP).NE.0 ) THEN
                        ipd = ITMA(IEXP,jgl)
                        DO jyi = 1 , idr
                           ni = KSEQ(jyi,3)
                           nf = KSEQ(jyi,4)
                           decen = EN(ni) - EN(nf)
                           cocos = SIN(ttttt)*SIN(gth)*COS(fm-figl)
     &                             + COS(ttttt)*COS(gth)
                           decen = decen*(1.+BETAR(IEXP)*cocos)
                           CALL EFFIX(ipd,decen,effi)
                           IF ( op2.EQ.'POIN' .AND. IPRM(20).EQ.1 )
     &                          WRITE (23,99049) ni , nf , SPIN(ni) , 
     &                                 SPIN(nf) , decen , effi
                           YGN(jyi) = YGN(jyi)*effi
                        ENDDO
                        inclus = ICLUST(IEXP,jgl)
                        IF ( inclus.NE.0 ) THEN
                           DO jyi = 1 , idr
                              SUMCL(inclus,jyi) = SUMCL(inclus,jyi)
     &                           + YGN(jyi)
                           ENDDO
                           IF ( jgl.NE.LASTCL(IEXP,inclus) ) GOTO 1305
                           DO jyi = 1 , idr
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
                              CORF(jmm,3) = YGN(jyi)
                              IF ( YGN(jyi).GE.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)/20.
                              IF ( YGN(jyi).LT.YGN(IDRN) ) CORF(jmm,4)
     &                             = CORF(jmm,3)
     &                             *(.05+.2*(1.-YGN(jyi)/YGN(IDRN)))
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
                                 IF ( jyi2.LE.0 ) GOTO 1302
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
 1302                ENDDO
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
 1305             ENDDO
                  IF ( op2.EQ.'CORR' ) THEN
                     jgl1 = 0
                     DO jgl = 1 , nogeli
                        IF ( IRAWEX(jexp).NE.0 ) THEN
                           inclus = ICLUST(jexp,jgl)
                           IF ( inclus.NE.0 ) THEN
                              IF ( jgl.NE.LASTCL(jexp,inclus) )
     &                             GOTO 1306
                           ENDIF
                        ENDIF
                        jgl1 = jgl1 + 1
                        READ (3,*) ne , na , zp , ap , xep , nval , waga
                        WRITE (4,*) ne , na , zp , ap , EP(IEXP) , 
     &                              nval , waga
                        WRITE (22,99031) IEXP , jgl1
99031                   FORMAT (///10X,'EXPERIMENT',1X,I2,8X,'DETECTOR',
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
                           ycorr = YEXP(jgl1,ile1+itp-1)
                           WRITE (4,*) ns1 , ns2 , ycorr , 
     &                                 DYEX(jgl1,ile1+itp-1)
                           WRITE (22,99032) ns1 , ns2 , 
     &                            CORF(ile1+itp-1,jgl1) , ycorr , 
     &                            ycorr/CORF(ile1+itp-1,jgl1)
99032                      FORMAT (5X,I4,5X,I4,3X,E8.3,4X,E8.3,4X,E8.3)
                        ENDDO ! Loop over itp
 1306                ENDDO ! Loop over jgl
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
            GOTO 100
         ENDIF ! if (op2 .NE. 'GOSI') if statement
      ENDIF ! if ( iobl.LT.1 ) if statement

 2400 IF ( iobl.GE.1 ) THEN
         ient = 1
         icg = 2
         nmaxh = NMAX
         lmax1 = LMAX
         sh1 = SPIN(1)
         sh2 = SPIN(2)
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
         DO kh = 1 , 8
            ihlm(kh) = MULTI(kh)
            ihlm(kh+24) = LDNUM(kh,2)
            ihlm(kh+8) = LAMDA(kh)
            ihlm(kh+16) = LDNUM(kh,1)
         ENDDO
         DO jexp = 1 , NEXPT
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
            DO jgs = 1 , LMAX
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
                     IF ( zmir(kk,1,IEXP).LT..5 ) zmir(kk,1,IEXP) = .5
                     IF ( zmir(kk,2,IEXP).LT..5 ) zmir(kk,2,IEXP) = .5
                  ENDIF
               ENDDO
            ENDDO
            DO kk = 1 , 6
               XIR(kk,IEXP) = XIR(kk,IEXP)*1.01
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
                  sz1 = MIN(zmir(kk,1,IEXP),10.)
                  sz2 = zmir(kk,2,IEXP)/50.
                  acof = 2.4009604E-3/zmir(kk,2,IEXP)
                  bcof = 8.163265E-4
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
                        IF ( jk.EQ.1 ) XI(1) = .02
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
                              ARM(1,5) = (.9999999,0.)
                              ARM(2,5) = (1.2E-6,0.)
                              ARM(1,6) = (.9999998,0.)
                              ARM(2,6) = (.9E-6,0.)
                              DO kh = 1 , 4
                                 ARM(1,kh) = (-1.E-6,0.)
                                 ARM(2,kh) = (1.E-6,0.)
                              ENDDO
                           ENDIF
                           CALL INTG(IEXP)
                           jp = 2
                           IF ( MAGA(IEXP).NE.0 .AND. jd.EQ.2 ) jp = 3
                           p = DBLE(ARM(1,5))
                           r = IMAG(ARM(1,5))
                           qr = DBLE(ARM(jp,5))
                           s = IMAG(ARM(jp,5))
                           test = p*p + r*r + qr*qr + s*s
                           p = p/SQRT(test)
                           s = ABS(r/s)
                           IF ( jk.EQ.1 ) THEN
                              IF ( MAGA(IEXP).EQ.0 ) THEN
                                 q1 = 0.
                                 GOTO 1402
                              ELSEIF ( jd.EQ.2 .OR. MAGA(IEXP).EQ.0 )
     &                                 THEN
                                 q1 = 0.
                                 GOTO 1402
                              ENDIF
                           ENDIF
                           q1 = ARCTG(s,ph1,pi)
                           ph1 = q1
 1402                      IF ( jk.EQ.1 ) THEN
                              IF ( jd.EQ.1 .AND. MAGA(IEXP).NE.0 ) THEN
                                 q2 = 0.
                                 GOTO 1404
                              ENDIF
                           ENDIF
                           q2 = ARCCOS(p,ph2,pi)
                           ph2 = q2
 1404                      q1 = q1/ZETA(jd)/2.
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
                           PARXM(IEXP,1,jk,kk) = acof*(2.*s12-51.*s11)
                           PARXM(IEXP,2,jk,kk) = bcof*(101.*s11-3.*s12)
                           PARXM(IEXP,3,jk,kk) = acof*(2.*s22-51.*s21)
                           PARXM(IEXP,4,jk,kk) = bcof*(101.*s21-3.*s22)
                        ENDIF
                     ENDDO ! Loop over jk
                  ENDDO ! Loop over jd
               ENDIF
            ENDDO ! Loop over kk
            EMMA(IEXP) = emhl1
            NMAX = nmaxh
            SPIN(1) = sh1
            SPIN(2) = sh2
            IFAC(1) = ih1
            IFAC(2) = ih2
            MAGEXC = magh
            ISO = isoh
            ELM(1) = eh1
            LEAD(1,1) = lh1
            LEAD(2,1) = lh2
            LAMMAX = lamh
            MEMAX = memh
            DO kh = 1 , 8
               LDNUM(kh,2) = ihlm(kh+24)
               MULTI(kh) = ihlm(kh)
               LAMDA(kh) = ihlm(kh+8)
               LDNUM(kh,1) = ihlm(kh+16)
            ENDDO
            INTERV(IEXP) = intvh
         ENDDO
         REWIND 7
         REWIND 27
         jkmx = 7
         IF ( IBPS.EQ.1 ) jkmx = 27
         DO iuy = 1 , 6
            WRITE (jkmx,*) (XIR(iuy,jj),jj=1,NEXPT)
            WRITE (jkmx,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  WRITE (jkmx,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               WRITE (jkmx,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
         jkmy = 12
         IF ( IBPS.EQ.1 ) jkmy = 32
         DO kuku = 1 , MEMAX
            WRITE (jkmy,*) ELM(kuku)
         ENDDO
         DO jj = 1 , 2
            DO jj1 = 1 , LP1
               IDIVE(jj1,jj) = 1
            ENDDO
         ENDDO
      ELSE
         REWIND 7
         REWIND 27
         mdupa = 7
         IF ( IBPS.EQ.1 ) mdupa = 27
         DO iuy = 1 , 6
            READ (mdupa,*) (XIR(iuy,jj),jj=1,NEXPT)
            READ (mdupa,*) (zmir(iuy,1,jj),zmir(iuy,2,jj),jj=1,NEXPT)
         ENDDO
         DO jj = 1 , NEXPT
            DO jk = 1 , 4
               DO kuku = 1 , 6
                  READ (mdupa,*) (PARXM(jj,jk,jl,kuku),jl=1,10)
               ENDDO
            ENDDO
            DO jk = 1 , 12
               READ (mdupa,*) (PARX(jj,jk,jl),jl=1,5)
            ENDDO
         ENDDO
         DO jgs = 1 , MEMAX
            DO jgr = 1 , 7
               QAPR(jgs,1,jgr) = 0.
            ENDDO
         ENDDO
      ENDIF
      IF ( IPRM(12).NE.0 .OR. op2.EQ.'MAP ' ) THEN
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
                     xxi = XIR(lex,jex)*(kex-1)/9.
                     WRITE (22,99055) xxi , 
     &                                (PARXM(jex,ilx,kex,lex),ilx=1,4)
                  ENDDO
                  IF ( MAGA(jex).NE.0 ) THEN
                     WRITE (22,99042) lex , zmir(lex,1,jex)
99042                FORMAT (1X//30X,'E',1I1,8X,'MI=+/-1',5X,
     &                       'MAX.ZETA=',1F6.3//)
                     WRITE (22,99054)
                     DO kex = 1 , 5
                        xxi = XIR(lex,jex)*(kex-1)/4.
                        u = 0.
                        WRITE (22,99055) xxi , u , PARX(jex,2*lex-1,kex)
     &                         , u , PARX(jex,2*lex,kex)
                     ENDDO ! Loop on kex
                  ENDIF ! if maga(jex).ne.0
               ENDIF ! if multi(lex).ne.0
            ENDDO ! Loop on lex
         ENDDO ! Loop on jex
      ENDIF ! IPRM(12).ne.0 .or. op2.eq.'MAP '

      IF ( op2.NE.'GOSI' .AND. op2.NE.'ERRO' ) GOTO 100
      IF ( op2.EQ.'ERRO' ) GOTO 500

 1400 mret = ABS(mret) ! Added for gosia2
      DO kh1 = 1 , MEMAX
         HLM(kh1) = ELM(kh1)
      ENDDO
      lfagg = 0
      DO kh1 = 1 , MEMAX
         IVAR(kh1) = ivarh(kh1)
      ENDDO
      IF ( nawr.EQ.1 ) THEN
         MCFIX = 0 ! Calculate chisq using CNOR1 not CNOR
         CALL MINI(chisq,chiok,nptl,conu,imode,idr,xtest,0,0,0,bten)
         IF ( IBPS.EQ.0 ) JZB = 25
         IF ( IBPS.EQ.0 ) mrepf = 1
         IF ( IBPS.NE.0 ) THEN
            IBPS = IBPS + 1
            IF ( IBPS.EQ.1 ) JZB = 26 ! But IBPS must be 2 here! Ooops!
            IF ( IBPS.EQ.1 ) GOTO 2200
         ENDIF
         IF ( IPS1.EQ.0 ) GOTO 2000
         IMIN = IMIN + 1
         DO iva = 1 , LP1 ! LP1 = 50
            JSKIP(iva) = 1
         ENDDO
         irix = 12
         IF ( IBPS.EQ.2 ) irix = 32 ! unit 12 for target, 32 for beam
         REWIND irix
         DO lkj = 1 , MEMAX
            WRITE (irix,*) ELM(lkj)
         ENDDO
         IF ( ifm.EQ.1 ) CALL PRELM(3)
         IF ( ifm.EQ.1 ) GOTO 2000 ! changed for gosia2
         GOTO 3000
      ELSE
         MCFIX = 1
         CALL MINI(chisq,1.D+38,nptl,conu,imode,idr,xtest,0,0,0,bten)
         mawr = mawr + 1
         IF ( JZB.EQ.25 ) ccch1 = chisq
         IF ( JZB.EQ.26 ) ccch2 = chisq
         chir = chir + chisq
         IF ( JZB.EQ.26 ) cht = ABS(chir-chp)
         IF ( JZB.EQ.26 .AND. cht.LT.0.1 ) mret = 0
         IF ( JZB.EQ.26 .AND. mawr.GE.20 ) mret = 0
         IF ( JZB.EQ.26 ) nawr = 1
         DO kh1 = 1 , 32
            DO kh2 = 1 , LP3
               IF ( JZB.EQ.26 ) CNOR2(kh1,kh2) = CNOR(kh1,kh2)
               IF ( JZB.EQ.25 ) CNOR1(kh1,kh2) = CNOR(kh1,kh2)
            ENDDO
         ENDDO
         IF ( JZB.EQ.26 ) THEN
            DO kh1 = 1 , 32
               DO kh2 = 1 , 75
c	cnor1(kh1,kh2)=cnor1(kh1,kh2)*cnor2(kh1,kh2)*2/(cnor1(kh1,kh2)
c     *+cnor2(kh1,kh2))
c      cnor1(kh1,kh2)=(cnor1(kh1,kh2)*ccch1+cnor2(kh1,kh2)*ccch2)/
c     *(ccch1+ccch2)
                  CNOR1(kh1,kh2) = (CNOR1(kh1,kh2)+CNOR2(kh1,kh2))/2.
               ENDDO
            ENDDO
         ENDIF
         IF ( IBPS.EQ.0 .AND. JZB.EQ.25 ) JZB = 26
         IF ( IBPS.NE.0 ) JZB = 25
         GOTO 2200
      ENDIF

 1500 WRITE (22,99043)
99043 FORMAT (5X,'ERROR-M.E. DOES NOT BELONG TO THE UPPER TRIANGLE')
      GOTO 1900

 1600 WRITE (22,99044)
99044 FORMAT (5X,'ERROR-WRONG SEQUENCE OF MULTIPOLARITIES')
      GOTO 1900

 1700 WRITE (22,99045)
99045 FORMAT (5X,'ERROR-REPEATED APPEARANCE OF THE STATE')
      GOTO 1900

 1800 WRITE (22,99046)
99046 FORMAT (1X///10X,'ERROR-INSUFFICIENT SPACE FOR E-THETA INTEGR ',
     &        'ATION')
      GOTO 1900

 3000 IF ( mret.EQ.1 .AND. JZB.EQ.26 ) nawr = 0
      IF ( mret.EQ.1 ) GOTO 2300
      IF ( IPRM(18).NE.0 ) CALL PTICC(idr)
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
                  DO iva1 = 1 , 10
                     IF ( LIFCT(iva1).EQ.iva ) GOTO 122
                  ENDDO
                  WRITE (22,99012) iva , TAU(iva)
99012             FORMAT (7X,1I2,7X,1E10.4)
                  GOTO 124
 122              WRITE (22,99013) iva , TAU(iva) , TIMEL(1,iva1) , 
     &                             TIMEL(2,iva1)
99013             FORMAT (7X,1I2,7X,1E10.4,5X,1E10.4,4X,1E10.4)
 124              IF ( iva.EQ.NMAX ) THEN
                     IF ( NAMX.GE.1 ) THEN
                        WRITE (22,99014)
99014                   FORMAT (5x,//,
     &                     'CALCULATED AND EXPERIMENTAL MATRIX ELEMENTS'
     &                     ,//)
                        WRITE (22,99015)
99015                   FORMAT (5x,'NI ','NF ',' EXP. ME   ',
     &                          'CURRENT ME','   SIGMA')
                        DO kq = 1 , NAMX
                           ni = IAMY(kq,1)
                           nf = IAMY(kq,2)
                           ind = IAMX(kq)
                           ess = ELM(ind)
                           esd = EAMX(kq,1)
                           dsd = EAMX(kq,2)
                           WRITE (22,99016) ni , nf , esd , ess , 
     &                            (ess-esd)/dsd
99016                      FORMAT (5x,1I2,1x,1I2,1x,1F9.4,1x,1F9.4,1x,
     &                             1F9.4)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF ( IMIN.NE.0 ) CALL PRELM(3)
         ENDIF
      ENDIF

 1900 IF ( ITS.NE.0 ) THEN
         iva = 0
         WRITE (18,*) iva , iva , iva , chisq
         IF ( ITS.NE.2 ) THEN
            WRITE (15,*) iva , chisq , chisq , chisq , chisq
            CALL KLOPOT(kmat,rlr)
         ENDIF
      ENDIF

 2000 WRITE (22,99047)
99047 FORMAT (15X,'********* END OF EXECUTION **********')

C     Decide if we have to loop again for beam/projectile
 2300 IF ( mrepf.NE.1 ) THEN
         IF ( mret.EQ.1 ) JZB = 25
         IF ( mret.EQ.1 ) GOTO 2200
         DO mmmm = 1 , 32
            DO kkkk = 1 , 50
               WRITE (13,*) CNOR1(mmmm,kkkk)
            ENDDO
         ENDDO
      ELSE
         mrepf = 2
         JZB = 26
         GOTO 2200
      ENDIF

99048 FORMAT (1X//50X,'CALCULATED YIELDS'//5X,'EXPERIMENT ',1I2,2X,
     &        'DETECTOR ',1I2/5X,'ENERGY ',1F10.3,1X,'MEV',2X,'THETA ',
     &        1F7.3,1X,'DEG'//5X,'NI',5X,'NF',5X,'II',5X,'IF',5X,
     &        'YIELD',5X,'NORMALIZED YIELD'/)
99049 FORMAT (5X,1I2,5X,1I2,3X,1F4.1,3X,1F4.1,3X,1E11.5,3X,1E11.5)
99050 FORMAT (1X///44X,'OVERALL')
99051 FORMAT (1X///43X,'DIAGONAL')
99052 FORMAT (6X,1I3,6X,1I2,5X,1I2,5X,1F10.5,2X,'(',1F10.5,' ,',1F10.5,
     &        ')')
99053 FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
99054 FORMAT (5X,'XI',13X,'Q1',22X,'Q2'///13X,'SLOPE',2X,'INTERCEPT',7X,
     &        'SLOPE',5X,'INTERCEPT'//)
99055 FORMAT (2X,1F6.4,3X,1E8.2,2X,1E8.2,6X,1E8.2,2X,1E8.2)
      END
