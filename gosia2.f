C
C       PROGRAM GOSIA2(INPUT,OUTPUT,TAPE6=OUTPUT,TAPE3,TAPE14,TAPE1,
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
C            DOUGLAS CLINE,NSRL,THE UNIVERSITY OF ROCHESTER
C            ROCHESTER,NY14627,USA------PHONE(716)275-4943
C
C            REF.----   UR/NSRL REPORT 308/1986
C
C************************  VERSION FROM SEPT, 2005  *******************
C
C**********************************************************************
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
     &          ispa , ispb , ITMA , itno , itp , ITS , ITTE , iuy , 
     &          iva , iva1 , IVAR , ivarh , ivari
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
      COMMON /SWITCH/ JZB , IBPS
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
      COMMON /RESC  / ELM25(500) , ELM26(500)
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
      DATA (eng(k),k=1,10)/.05 , .06 , .08 , .1 , .15 , .2 , .3 , .5 , 
     &      1. , 1.5/
      DATA (tau1(k),k=1,10)/17.656 , 10.726 , 5.076 , 2.931 , 1.3065 , 
     &      .8828 , .5959 , .4357 , .3041 , .2472/
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
      READ * , IBPS
      READ 99055 , op1 , op2
      CALL OPENF1
      IBPS = IBPS - 1
      JZB = 25
      nawr = 0
      IF ( IBPS.EQ.1 ) JZB = 26
      DO i = 1 , 75
         DO j = 1 , 32
            CNOR1(j,i) = 1.
            CNOR2(j,i) = 1.
            CNOR(j,i) = 1.
         ENDDO
      ENDDO
      MCFIX = 1
      mres1 = 0
      mrepf = 0
      mret = -1
      chir = 0.
      mawr = 0
 100  DO i = 1 , 50
         LIFCT(i) = 0
      ENDDO
      chp = chir
      IF ( mret.EQ.1 .AND. JZB.EQ.26 ) chir = 0.
      REWIND 25
      REWIND 26
      READ (JZB,*) IBPS
      IBPS = IBPS - 1
      JZB = 25
      IF ( IBPS.EQ.1 ) JZB = 26
      IBYP = 0
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
      ifm = 0
      IPS1 = 11
      ifwd = -1
      INTR = 0
      LNY = 0
      JENTR = 0
      lp0 = 50000
      ICS = 0
      LP1 = 50
      LP2 = 500
      LP3 = 75
      LP4 = 1500
      LP6 = 32
      LP7 = lp0 - 4900
      LP8 = LP3*28 + 1
      LP9 = lp0 - LP3*28
      LP10 = 600
      LP11 = LP8 - 1
      LP12 = 365
      LP13 = LP9 + 1
      LP14 = 4900
      DO i = 1 , LP1
         jpin(i) = 0
         iecd(i) = 0
      ENDDO
      txx = 0.
      SGW = 3.
      SUBCH1 = 0.
      SUBCH2 = 0.
      ITS = 0
      iosr = 0
      LOCKS = 0
      DLOCK = 1.1
      kerf = 0
      IFBFL = 0
      NLOCK = 0
      LOCKF = 0
      DO i = 1 , LP4
         DO j = 1 , LP6
            CORF(i,j) = 1.
         ENDDO
      ENDDO
      DO i = 1 , 20
         IPRM(i) = 1
         DO j = 1 , 5
            CC(i,j) = 0.
         ENDDO
      ENDDO
      IPRM(2) = 0
      IPRM(4) = 0
      IPRM(5) = 1
      IPRM(6) = 11111
      IPRM(7) = 0
      IPRM(16) = 0
      IPRM(17) = 0
      IPRM(18) = 0
      IPRM(19) = 0
      IPRM(20) = 0
      IPRM(11) = 0
      IPRM(12) = 0
      IPRM(14) = 0
      DO i = 1 , LP1
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
      DO k = 1 , LP1
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
      NDIM = LP3
      ISO = 1
      B(1) = 1.
      DO i = 2 , 20
         B(i) = B(i-1)*(i-1)
      ENDDO
      LMAXE = 0
      CALL FAKP
      CALL FHIP
      NCM = 2
      DO ijx = 1 , LP1
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
         DO lexp = 1 , LP3
            LDNUM(lam,lexp) = 0
         ENDDO
         MULTI(lam) = 0
         LAMDA(lam) = 0
      ENDDO
      DO j = 1 , LP2
         EXPO(j) = (1.,0.)
         KVAR(j) = 1
         ELM(j) = 0.
      ENDDO
      DO j = 1 , LP1
         JSKIP(j) = 1
         ISKIN(j) = 0
      ENDDO
      DO j = 1 , LP3
         ISEX(j) = 1111
      ENDDO
      ISEX(1) = 0
      ACCA = .00001
      oph = '    '
      nmemx = LP2 + 9
      IEXP = 1
      IMIN = 0
      i122 = 0
      DO j = 1 , LP2
         DO k = 1 , 2
            DO l = 1 , 7
               QAPR(j,k,l) = 0.
            ENDDO
         ENDDO
      ENDDO
      ERR = .FALSE.
      intend = 0
 200  READ (JZB,99055) op1 , op2
      IF ( op1.EQ.'OP, ' ) THEN
         IF ( op2.EQ.'GOSI' ) oph = op2
         IF ( op2.EQ.'GOSI' ) opcja = op2
         IF ( op2.EQ.'FILE' ) CALL OPENF
         IF ( op2.EQ.'FILE' ) GOTO 200
         IF ( jphd.EQ.1 ) WRITE (22,99001)
99001    FORMAT ('1'/1X,125('*')/1X,125('*')/1X,50('*'),25X,50('*')/1X,
     &           50('*'),10X,'GOSIA2',10X,50('*')/1X,50('*'),25X,50('*')
     &           /1X,125('*')/1X,125('*')////)
         IF ( jphd.EQ.1 ) WRITE (22,99002)
99002    FORMAT (1X/20X,'ROCHESTER COULOMB EXCITATION DATA ANALYSIS ',
     &           'CODE BY T.CZOSNYKA,D.CLINE AND C.Y.WU'/50X,
     &           'LATEST REVISION- JUNE  2006'//////)
         jphd = 0
         IF ( op2.EQ.'GDET' ) THEN
            nl = 7
            READ (JZB,*) nfdd
            nfd = IABS(nfdd)
            IF ( nfdd.LE.0 ) THEN
               REWIND 8
               DO i = 1 , nl
                  WRITE (8,*) (tau2(l,i),l=1,10)
               ENDDO
               WRITE (8,*) (eng(l),l=1,10)
            ENDIF
            REWIND 9
            WRITE (9,*) nfd
            DO i = 1 , nfd
               READ (JZB,*) (DIX(k),k=1,4)
               READ (JZB,*) (xl1(k),k=1,nl)
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
            GOTO 200
         ELSEIF ( op2.EQ.'RAND' ) THEN
            READ (JZB,*) SE
            CALL MIXUP
            WRITE (22,99006)
99006       FORMAT (1X///5X,'MATRIX ELEMENTS RANDOMIZED...'///)
            CALL PRELM(2)
            GOTO 200
         ELSEIF ( op2.EQ.'TROU' ) THEN
            ITS = 1
            READ (JZB,*) kmat , rlr
            GOTO 200
         ELSEIF ( op2.EQ.'REST' ) THEN
            irix = 12
            IF ( IBPS.EQ.1 ) irix = 32
            REWIND irix
            memax1 = MEMAX + 1
            DO lkj = 1 , MEMAX
               READ (irix,*) ELM(lkj)
            ENDDO
            DO lkj = 1 , memax1
               READ (JZB,*) lkj1 , xlk
               IF ( lkj1.EQ.0 ) GOTO 220
               IF ( mres1.EQ.0 ) ELM(lkj1) = xlk
            ENDDO
 220        DO lkj = 1 , MEMAX
               IF ( mres1.EQ.0 .AND. JZB.EQ.25 ) ELM25(lkj) = ELM(lkj)
               IF ( mres1.EQ.0 .AND. JZB.EQ.26 ) ELM26(lkj) = ELM(lkj)
            ENDDO
            IF ( JZB.EQ.26 ) mres1 = 1
            WRITE (22,99007)
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
            CALL PRELM(4)
            GOTO 200
         ELSE
            IF ( op2.EQ.'RE,A' ) GOTO 1000
            IF ( op2.EQ.'RE,F' ) GOTO 1000
            IF ( op2.EQ.'ERRO' ) THEN
               READ (JZB,*) idf , ms , mend , irep , ifc , remax
               MCFIX = 0
               DO mmmm = 1 , 32
                  DO kkkk = 1 , 50
                     READ (13,*) CNOR1(mmmm,kkkk)
                  ENDDO
               ENDDO
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
                     WRITE (33,*) NMAX , MEMAX , inpo , inko
                     DO inn = 1 , NMAX
                        WRITE (33,*) inn , SPIN(inn) , EN(inn)
                     ENDDO
                     DO inn = 1 , MEMAX
                        WRITE (33,*) inn , LEAD(1,inn) , LEAD(2,inn)
                     ENDDO
                     DO inn = 1 , MEMAX
                        WRITE (33,*) inn , ELM(inn)
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
               IF ( ERR ) GOTO 2200
               IF ( IMIN.NE.0 ) GOTO 500
               GOTO 1400
            ELSEIF ( op2.EQ.'RE,C' ) THEN
               jfre = 1
               irfix = 0
               GOTO 1100
            ELSEIF ( op2.EQ.'TITL' ) THEN
               READ (JZB,99008) (title(k),k=1,20)
99008          FORMAT (20A4)
               WRITE (22,99009) (title(k),k=1,20)
99009          FORMAT (10X,20A4/10X,100('-'))
               GOTO 200
            ELSE
               IF ( op2.EQ.'GOSI' ) GOTO 300
               IF ( op2.EQ.'COUL' ) GOTO 300
               IF ( op2.EQ.'EXIT' ) GOTO 2000
               IF ( op2.EQ.'MINI' ) THEN
                  READ (JZB,*) imode , nptl , chiok , conu , xtest , 
     &                         LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
                  op2 = opcja
                  IMIN = IMIN + 1
                  IF ( IMIN.NE.1 ) GOTO 1500
                  GOTO 1300
               ELSEIF ( op2.EQ.'THEO' ) THEN
                  REWIND (12)
                  ibaf = 1
                  DO jb = 1 , LP1
                     DO lb = 1 , 2
                        xlevb(jb,lb) = 0
                     ENDDO
                  ENDDO
                  READ (JZB,*) nbands
                  IF ( nbands.LE.0 ) ibaf = 0
                  nbands = IABS(nbands)
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
                     READ (JZB,*) bk , ilevls
                     READ (JZB,*) (levl(ib),ib=1,ilevls)
                     DO kb = 1 , ilevls
                        inva = levl(kb)
                        xlevb(inva,2) = bk
                        xlevb(inva,1) = FLOAT(jb)
                     ENDDO
                  ENDDO
                  DO nl = 1 , 8
                     READ (JZB,*) nnl
 222                 IF ( nnl.LE.0 ) GOTO 225
                     READ (JZB,*) jb1 , jb2
                     IF ( jb1.NE.0 ) THEN
                        READ (JZB,*) (bm(nnl,jb1,jb2,j),j=1,3)
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
                        IF ( IBPS.EQ.1 ) irix = 32
                        WRITE (irix,*) ELM(kb)
                     ENDIF
                  ENDDO
                  GOTO 200
               ELSEIF ( op2.EQ.'YIEL' ) THEN
                  CALL ADHOC(oph,idr,nfd,ntap,iyr)
                  GOTO 200
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
                        IF ( iecd(lx).EQ.1 ) THEN
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
                        IF ( iecd(lx).EQ.1 ) THEN
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
                              IF ( ERR ) GOTO 2200
                              tting = TLBDG(lx)
                              IF ( ERR ) GOTO 2100
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
                                    WRITE (22,99047) lx , ija0 , enb , 
     &                                 tta
                                    IF ( tta.LT.0. ) WRITE (22,99010)
     &                                 tting
99010                               FORMAT (5X,
     &                             'RESPECTIVE TARGET SCATTERING ANGLE='
     &                             ,1F7.3,1X,'DEG'/)
                                    DO jyi = 1 , idr
                                       ni = KSEQ(jyi,3)
                                       nf = KSEQ(jyi,4)
                                       WRITE (22,99048) ni , nf , 
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
                        IF ( iocc.GT.izcap ) GOTO 1900
                        hen = (emx-emn)/npce
                        npce1 = npce + 1
                        het = (tmx-tmn)/npct
                        npct1 = npct + 1
                        IF ( iecd(lx).EQ.1 )
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
                              IF ( nft.EQ.1 ) GOTO 2100
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
                           IF ( ja.EQ.1 ) WRITE (22,99011) DS , lx
99011                      FORMAT (1X/////5X,
     &                            'INTEGRATED RUTHERFORD CROSS SECTION='
     &                            ,1E9.4,2X,'FOR EXP.',1I2///)
                           WRITE (22,99012) lx , ja , emn , emx , tmn , 
     &                            tmx
99012                      FORMAT (1X,//50X,'INTEGRATED YIELDS'//5X,
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
                              WRITE (22,99048) ni , nf , SPIN(ni) , 
     &                               SPIN(nf) , GRAD(jd) , GRAD(jd)
     &                               /GRAD(IDRN)
                           ENDDO
                        ENDDO
                        IF ( iecd(lx).EQ.1 ) THEN
                           IF ( jpin(lx).EQ.0 ) THEN
                              CALL COORD(wth,wph,wthh,1,2,pfi,wpi,
     &                           TLBDG(lx),lx,txx,txx)
                              WRITE (22,99013) FIEX(lx,1)*57.2957795 , 
     &                               FIEX(lx,2)*57.2957795 , lx
99013                         FORMAT (//5X,
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
                     IF ( mpin.GT.1 ) WRITE (22,99014) dst , lx
99014                FORMAT (1x//2x,
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
                  GOTO 200
               ELSEIF ( op2.EQ.'CORR' ) THEN
                  CALL READY(idr,ntap,0)
                  REWIND 3
                  REWIND 15
                  REWIND 4
                  GOTO 1300
               ELSE
                  IF ( op2.EQ.'POIN' ) GOTO 1300
                  IF ( op2.EQ.'MAP ' ) iobl = 1
                  IF ( op2.EQ.'STAR' ) GOTO 1300
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
                     GOTO 2200
                  ELSEIF ( op2.EQ.'RAW ' ) THEN
                     REWIND 8
                     DO l = 1 , 8
                        READ (8,*) (ABC(l,j),j=1,10)
                        DO j = 1 , 10
                           ABC(l,j) = LOG(ABC(l,j))
                        ENDDO
                     ENDDO
                     DO l = 1 , nfd
                        READ (8,*) (THICK(l,j),j=1,7)
                     ENDDO
                     DO l = 1 , LP1
                        DO j = 1 , 200
                           ICLUST(l,j) = 0
                        ENDDO
                        DO j = 1 , 20
                           LASTCL(l,j) = 0
                        ENDDO
                        IRAWEX(l) = 0
                     ENDDO
                     DO l = 1 , LP1
                        READ (JZB,*) mexl
                        IF ( mexl.EQ.0 ) GOTO 200
                        IRAWEX(mexl) = 1
                        n = NANG(mexl)
                        DO j = 1 , n
                           jj = ITMA(mexl,j)
                           READ (JZB,*) (AKAVKA(k,jj),k=1,8)
                        ENDDO
                        READ (JZB,*) kclust
                        IF ( kclust.NE.0 ) THEN
                           DO j = 1 , kclust
                              READ (JZB,*) numcl
                              READ (JZB,*) (liscl(k),k=1,numcl)
                              LASTCL(l,j) = liscl(numcl)
                              DO k = 1 , numcl
                                 kk = liscl(k)
                                 ICLUST(l,kk) = j
                              ENDDO
                           ENDDO
                        ENDIF
                     ENDDO
                     GOTO 200
                  ELSEIF ( op2.EQ.'MAP ' ) THEN
                     GOTO 1300
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      WRITE (22,99015) op1 , op2
99015 FORMAT (5X,'UNRECOGNIZED OPTION',1X,1A3,1A4)
      GOTO 2200
 300  READ (JZB,99016) op1
99016 FORMAT (1A4)
      IF ( op1.EQ.'    ' ) GOTO 200
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
         GOTO 300
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
                  IF ( ipo1.LE.la ) GOTO 1700
                  LAMMAX = LAMMAX + 1
                  LAMDA(LAMMAX) = ipo1
                  ipo3 = 0
                  IF ( indx.EQ.0 ) GOTO 320
               ELSE
                  MULTI(la) = MULTI(la) + 1
                  indx = indx + 1
                  IF ( ipo1.GT.ABS(ipo2) ) GOTO 1600
                  IF ( ipo1.NE.ipo3 ) THEN
                     IF ( ipo1.LT.ipo3 ) GOTO 1800
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
         GOTO 300
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
                              IF ( op1.NE.'END,' ) GOTO 450
                              GOTO 300
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
         GOTO 300
      ELSE
         WRITE (22,99020) op1
99020    FORMAT (5X,'UNRECOGNIZED SUBOPTION',1X,1A4)
         GOTO 2200
      ENDIF
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
      IF ( naxfl.EQ.0 ) WRITE (22,99050)
      IF ( naxfl.NE.0 ) WRITE (22,99049)
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
      IF ( naxfl.NE.0 ) WRITE (22,99049)
      IF ( naxfl.EQ.0 ) WRITE (22,99050)
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
               WRITE (22,99051) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , be2 , 
     &                          be2a - be2 , be2b - be2
            ELSE
               ispb = INT(SPIN(ispa))*2
               qfac = 3.170662*WTHREJ(ispb,4,ispb,-ispb,0,ispb)
               WRITE (22,99051) kh2 , LEAD(2,kh2) , LEAD(1,kh2) , 
     &                          HLM(kh2)*qfac , DEVD(kh2)*qfac , 
     &                          DEVU(kh2)*qfac
            ENDIF
         ENDIF
      ENDDO
      GOTO 2200
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
                  CALL MINI(chiss,-1.,2,.0001,1000,idr,100000.,0,iosr,
     &                      kh,bten)
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
      GOTO 200
 1300 CALL CMLAB(0,dsig,ttttt)
      IF ( ERR ) GOTO 2200
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
               IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 ) WRITE (22,99052)
     &              jgl , pop1
               DO jgl = 2 , NMAX
                  loct = (jgl-1)*28 + 1
                  IF ( op2.EQ.'STAR' .OR. IPRM(19).EQ.1 )
     &                 WRITE (22,99052) jgl , ZETA(loct)
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
     &                          WRITE (23,99048) ni , nf , SPIN(ni) , 
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
     &                    WRITE (22,99047) IEXP , jgl1 , EP(IEXP) , 
     &                    TLBDG(IEXP)
                     jmm = 0
                     DO jyi = 1 , idr
                        ni = KSEQ(jyi,3)
                        nf = KSEQ(jyi,4)
                        IF ( op2.EQ.'POIN' .OR. IPRM(11).EQ.1 )
     &                       WRITE (22,99048) ni , nf , SPIN(ni) , 
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
                        ENDDO
 1306                ENDDO
                  ENDIF
               ENDIF
            ENDDO
            IF ( op2.EQ.'STAR' ) oph = op2
            IF ( op2.NE.'STAR' ) THEN
               IF ( op2.EQ.'CORR' ) THEN
                  ntap = 4
                  CALL READY(idr,ntap,ipri)
                  REWIND ntap
               ENDIF
            ENDIF
            GOTO 200
         ENDIF
      ENDIF
 1400 IF ( iobl.GE.1 ) THEN
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
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
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
                  WRITE (22,99033) jex , XIR(lex,jex)
99033             FORMAT (1X//30X,'EXPERIMENT',1X,1I2,10X,'MAX.XI=',
     &                    1F6.4)
                  WRITE (22,99034) lex , zmir(lex,2,jex)
99034             FORMAT (1X/30X,'E',1I1,8X,'MI=0',5X,'MAX.ZETA=',
     &                    1F6.3//)
                  WRITE (22,99053)
                  DO kex = 1 , 10
                     xxi = XIR(lex,jex)*(kex-1)/9.
                     WRITE (22,99054) xxi , 
     &                                (PARXM(jex,ilx,kex,lex),ilx=1,4)
                  ENDDO
                  IF ( MAGA(jex).NE.0 ) THEN
                     WRITE (22,99035) lex , zmir(lex,1,jex)
99035                FORMAT (1X//30X,'E',1I1,8X,'MI=+/-1',5X,
     &                       'MAX.ZETA=',1F6.3//)
                     WRITE (22,99053)
                     DO kex = 1 , 5
                        xxi = XIR(lex,jex)*(kex-1)/4.
                        u = 0.
                        WRITE (22,99054) xxi , u , PARX(jex,2*lex-1,kex)
     &                         , u , PARX(jex,2*lex,kex)
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      IF ( op2.NE.'GOSI' .AND. op2.NE.'ERRO' ) GOTO 200
      IF ( op2.EQ.'ERRO' ) GOTO 500
 1500 mret = IABS(mret)
      DO kh1 = 1 , MEMAX
         HLM(kh1) = ELM(kh1)
      ENDDO
      lfagg = 0
      DO kh1 = 1 , MEMAX
         IVAR(kh1) = ivarh(kh1)
      ENDDO
      IF ( nawr.EQ.1 ) THEN
         MCFIX = 0
         CALL MINI(chisq,chiok,nptl,conu,imode,idr,xtest,0,0,0,bten)
         IF ( IBPS.EQ.0 ) JZB = 25
         IF ( IBPS.EQ.0 ) mrepf = 1
         IF ( IBPS.NE.0 ) THEN
            IBPS = IBPS + 1
            IF ( IBPS.EQ.1 ) JZB = 26
            IF ( IBPS.EQ.1 ) GOTO 100
         ENDIF
         IF ( IPS1.EQ.0 ) GOTO 2200
         IMIN = IMIN + 1
         DO iva = 1 , LP1
            JSKIP(iva) = 1
         ENDDO
         irix = 12
         IF ( IBPS.EQ.2 ) irix = 32
         REWIND irix
         DO lkj = 1 , MEMAX
            WRITE (irix,*) ELM(lkj)
         ENDDO
         IF ( ifm.EQ.1 ) CALL PRELM(3)
         IF ( ifm.NE.1 ) GOTO 2000
         GOTO 2200
      ELSE
         MCFIX = 1
         CALL MINI(chisq,1.E+38,nptl,conu,imode,idr,xtest,0,0,0,bten)
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
         GOTO 100
      ENDIF
 1600 WRITE (22,99036)
99036 FORMAT (5X,'ERROR-M.E. DOES NOT BELONG TO THE UPPER TRIANGLE')
      GOTO 2100
 1700 WRITE (22,99037)
99037 FORMAT (5X,'ERROR-WRONG SEQUENCE OF MULTIPOLARITIES')
      GOTO 2100
 1800 WRITE (22,99038)
99038 FORMAT (5X,'ERROR-REPEATED APPEARANCE OF THE STATE')
      GOTO 2100
 1900 WRITE (22,99039)
99039 FORMAT (1X///10X,'ERROR-INSUFFICIENT SPACE FOR E-THETA INTEGR ',
     &        'ATION')
      GOTO 2100
 2000 IF ( mret.EQ.1 .AND. JZB.EQ.26 ) nawr = 0
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
               WRITE (22,99040)
99040          FORMAT (1X//20X,'CALCULATED LIFETIMES'//5X,'LEVEL',5X,
     &                 'LIFETIME(PSEC)',5X,'EXP',8X,'ERROR'/)
               DO iva = 2 , NMAX
                  DO iva1 = 1 , 10
                     IF ( LIFCT(iva1).EQ.iva ) GOTO 2005
                  ENDDO
                  WRITE (22,99041) iva , TAU(iva)
99041             FORMAT (7X,1I2,7X,1E10.4)
                  GOTO 2010
 2005             WRITE (22,99042) iva , TAU(iva) , TIMEL(1,iva1) , 
     &                             TIMEL(2,iva1)
99042             FORMAT (7X,1I2,7X,1E10.4,5X,1E10.4,4X,1E10.4)
 2010             IF ( iva.EQ.NMAX ) THEN
                     IF ( NAMX.GE.1 ) THEN
                        WRITE (22,99043)
99043                   FORMAT (5x,//,
     &                     'CALCULATED AND EXPERIMENTAL MATRIX ELEMENTS'
     &                     ,//)
                        WRITE (22,99044)
99044                   FORMAT (5x,'NI ','NF ',' EXP. ME   ',
     &                          'CURRENT ME','   SIGMA')
                        DO kq = 1 , NAMX
                           ni = IAMY(kq,1)
                           nf = IAMY(kq,2)
                           ind = IAMX(kq)
                           ess = ELM(ind)
                           esd = EAMX(kq,1)
                           dsd = EAMX(kq,2)
                           WRITE (22,99045) ni , nf , esd , ess , 
     &                            (ess-esd)/dsd
99045                      FORMAT (5x,1I2,1x,1I2,1x,1F9.4,1x,1F9.4,1x,
     &                             1F9.4)
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF ( IMIN.NE.0 ) CALL PRELM(3)
         ENDIF
      ENDIF
 2100 IF ( ITS.NE.0 ) THEN
         iva = 0
         WRITE (18,*) iva , iva , iva , chisq
         IF ( ITS.NE.2 ) THEN
            WRITE (15,*) iva , chisq , chisq , chisq , chisq
            CALL KLOPOT(kmat,rlr)
         ENDIF
      ENDIF
 2200 WRITE (22,99046)
99046 FORMAT (15X,'********* END OF EXECUTION **********')
 2300 IF ( mrepf.NE.1 ) THEN
         IF ( mret.EQ.1 ) JZB = 25
         IF ( mret.EQ.1 ) GOTO 100
         DO mmmm = 1 , 32
            DO kkkk = 1 , 50
               WRITE (13,*) CNOR1(mmmm,kkkk)
            ENDDO
         ENDDO
      ELSE
         mrepf = 2
         JZB = 26
         GOTO 100
      ENDIF
99047 FORMAT (1X//50X,'CALCULATED YIELDS'//5X,'EXPERIMENT ',1I2,2X,
     &        'DETECTOR ',1I2/5X,'ENERGY ',1F10.3,1X,'MEV',2X,'THETA ',
     &        1F7.3,1X,'DEG'//5X,'NI',5X,'NF',5X,'II',5X,'IF',5X,
     &        'YIELD',5X,'NORMALIZED YIELD'/)
99048 FORMAT (5X,1I2,5X,1I2,3X,1F4.1,3X,1F4.1,3X,1E11.5,3X,1E11.5)
99049 FORMAT (1X///44X,'OVERALL')
99050 FORMAT (1X///43X,'DIAGONAL')
99051 FORMAT (6X,1I3,6X,1I2,5X,1I2,5X,1F10.5,2X,'(',1F10.5,' ,',1F10.5,
     &        ')')
99052 FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
99053 FORMAT (5X,'XI',13X,'Q1',22X,'Q2'///13X,'SLOPE',2X,'INTERCEPT',7X,
     &        'SLOPE',5X,'INTERCEPT'//)
99054 FORMAT (2X,1F6.4,3X,1E8.2,2X,1E8.2,6X,1E8.2,2X,1E8.2)
99055 FORMAT (1A3,1A4)
      END
      FUNCTION ARCCOS(A,F,Pi)
      IMPLICIT NONE
      REAL*8 A , an , ARCCOS , F , Pi , q , qa , qap , TACOS
      INTEGER*4 ie , j , k
      q = TACOS(A)
      qa = q
      qap = q
      IF ( q.LE.F ) THEN
         DO j = 1 , 20
            an = 2*j*Pi
            DO k = 1 , 2
               qap = qa
               ie = (-1)**k
               qa = an + ie*q
               IF ( qa.GT.F ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
 100  ARCCOS = qa
      IF ( (qa-F).GT.Pi/2. ) ARCCOS = qap
      END
      FUNCTION ARCTG(A,F,Pi)
      IMPLICIT NONE
      REAL*8 A , an , ARCTG , F , Pi , q , qa , qap
      INTEGER*4 ie , j , k
      q = ATAN(A)
      qa = q
      qap = q
      IF ( q.LE.F ) THEN
         DO j = 1 , 40
            an = j*Pi
            DO k = 1 , 2
               qap = qa
               ie = (-1)**k
               qa = an + ie*q
               IF ( qa.GT.F ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
 100  ARCTG = qa
      IF ( (qa-F).GT.Pi/4. ) ARCTG = qap
      END
      SUBROUTINE LOAD(Iexp,Ient,Icg,Polm,Joj)
      IMPLICIT NONE
      REAL*8 a1 , a2 , aaz2 , aaz3 , aazz , ACCA , ACCUR , ah , CAT , 
     &       cpsi , dep , DIPOL , EMMA , EN , EP , eta , etan , Polm , 
     &       pp1 , pp2
      REAL*8 ppp , PSI , QAPR , rlam , SPIN , ssqrt , szet , TLBDG , 
     &       VINF , wrt , wrtm , XA , XA1 , XI , z1 , z2 , zet , ZETA , 
     &       ZPOL , zsqa
      INTEGER*4 i , i1 , i2 , i3 , IAPR , Icg , Ient , Iexp , IPATH , 
     &          ir , is , ISEX , ISHA , ISMAX , ISO , ispi , ispo , 
     &          IVAR , IZ , IZ1
      INTEGER*4 jj , jjj , Joj , la , lam , lam1 , LAMDA , LAMMAX , ld , 
     &          LDNUM , LEAD , LMAX , LMAXE , LP1 , LP10 , LP11 , LP12 , 
     &          LP13 , LP14 , LP2
      INTEGER*4 LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , LZETA , m , m1 , 
     &          m2 , MAGA , MAGEXC , MEMAX , MEMX6 , mstop , MULTI , n , 
     &          n2 , n3 , NCM
      INTEGER*4 NDIM , NEXPT , NMAX , NMAX1 , nn , NSTART , NSTOP , nz
      LOGICAL ERR
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /PSPIN / ISHA
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /PCOM  / PSI(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CLM   / LMAX
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /CLCOM9/ ERR
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /CXI   / XI(500)
      COMMON /CAUX0 / EMMA(75) , NCM
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      DIMENSION etan(75) , cpsi(8)
      LMAX = INT(SPIN(1)+1.1)
      IF ( Ient.EQ.1 ) THEN
         ISHA = 0
         ispi = INT(SPIN(1)+.51)
         ispo = INT(SPIN(1)+.49)
         IF ( ispi.NE.ispo ) ISHA = 1
         z1 = DBLE(ABS(IZ1(Iexp)))
         z2 = DBLE(IZ)
         a1 = XA1(Iexp)
         a2 = XA
         ZPOL = DIPOL*EP(Iexp)*a2/(z2*z2*(1.+a1/a2))
         IF ( IZ1(Iexp).LT.0 ) ZPOL = DIPOL*EP(Iexp)
     &                                *a1/(z1*z1*(1.+a2/a1))
         IF ( IZ1(Iexp).LE.0 ) THEN
            ah = a1
            a1 = a2
            a2 = ah
         ENDIF
         eta = z1*z2*SQRT(a1/EP(Iexp))/6.349770
         DO m = 1 , NMAX
            dep = (1.0+a1/a2)*EN(m)
            zet = dep/EP(Iexp)
            szet = SQRT(1.0-zet)
            etan(m) = eta/szet
         ENDDO
         DO n = 1 , MEMAX
            i1 = LEAD(1,n)
            i2 = LEAD(2,n)
            XI(n) = etan(i1) - etan(i2)
         ENDDO
         aazz = 1./(1.+a1/a2)/z1/z2
         cpsi(1) = 5.169286*aazz
         IF ( LMAXE.NE.1 ) THEN
            aaz2 = aazz*aazz
            cpsi(2) = 14.359366*aaz2
            IF ( LMAXE.NE.2 ) THEN
               aaz3 = aazz*aaz2
               cpsi(3) = 56.982577*aaz3
               IF ( LMAXE.NE.3 ) THEN
                  aazz = aaz2*aaz2
                  cpsi(4) = 263.812653*aazz
                  IF ( LMAXE.NE.4 ) THEN
                     aaz2 = aaz3*aaz2
                     cpsi(5) = 1332.409500*aaz2
                     IF ( LMAXE.NE.5 ) THEN
                        aazz = aaz3*aaz3
                        cpsi(6) = 7117.691577*aazz
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF ( MAGEXC.NE.0 ) THEN
            aazz = VINF(Iexp)/95.0981942
            cpsi(7) = aazz*cpsi(1)
            IF ( LAMMAX.NE.8 ) cpsi(8) = aazz*cpsi(2)
         ENDIF
         zsqa = z1*SQRT(a1)
         i3 = 1
         ppp = 1. + a1/a2
         DO i1 = 1 , LAMMAX
            lam = LAMDA(i1)
            lam1 = lam
            IF ( lam.GT.6 ) lam1 = lam - 6
            DO n2 = 1 , NMAX
               nn = LDNUM(lam,n2)
               IF ( nn.NE.0 ) THEN
                  n3 = LEAD(1,i3)
                  pp1 = EP(Iexp) - ppp*EN(n3)
                  DO m1 = 1 , nn
                     m2 = LEAD(2,i3)
                     i2 = i3
                     i3 = i3 + 1
                     pp2 = EP(Iexp) - ppp*EN(m2)
                     PSI(i2) = cpsi(lam)*zsqa*(pp1*pp2)
     &                         **((2.*DBLE(lam1)-1.)/4.)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         IF ( Ient.EQ.1 ) RETURN
      ENDIF
      DO n = 1 , NMAX
         NSTART(n) = 0
         NSTOP(n) = 0
      ENDDO
      is = 1
      NSTART(1) = 1
      DO n = 1 , NMAX
         wrt = Polm - EMMA(Iexp)
         wrtm = Polm + EMMA(Iexp)
         IF ( Icg.EQ.2 ) wrt = Polm - DBLE(MAGA(Iexp))
         IF ( Icg.EQ.2 ) wrtm = Polm + DBLE(MAGA(Iexp))
         IF ( wrtm.LT.-SPIN(n) ) THEN
            NSTART(n) = 0
         ELSE
            IF ( ABS(wrt).GT.SPIN(n) ) wrt = -SPIN(n)
            IF ( wrtm.GT.SPIN(n) ) wrtm = SPIN(n)
            mstop = INT(wrtm-wrt+1.01)
            DO i = 1 , mstop
               CAT(is,1) = n
               CAT(is,2) = SPIN(n)
               CAT(is,3) = wrt + DBLE(i-1)
               IF ( n.EQ.1 .AND. ABS(CAT(is,3)-Polm).LT.1.E-6 ) Joj = is
               is = is + 1
            ENDDO
         ENDIF
         NSTART(n+1) = is
         NSTOP(n) = is - 1
      ENDDO
      ISMAX = is - 1
      IF ( ISMAX.LE.LP10 ) THEN
         IF ( Ient.EQ.3 ) RETURN
         nz = 0
         DO jj = 1 , 7
            DO jjj = 1 , MEMAX
               QAPR(jjj,1,jj) = 0.
               QAPR(jjj,2,jj) = 0.
            ENDDO
         ENDDO
         DO i = 1 , 8
            LZETA(i) = 0
         ENDDO
         DO i1 = 1 , LAMMAX
            lam = LAMDA(i1)
            IF ( Icg.NE.2 .OR. lam.LE.6 ) THEN
               la = lam
               IF ( lam.GT.6 ) lam = lam - 6
               rlam = DBLE(lam)
               ssqrt = SQRT(2.*rlam+1.)
               LZETA(la) = nz
               ir = 0
 10            ir = ir + 1
               IF ( ir.LE.ISMAX ) THEN
                  n = CAT(ir,1)
                  IF ( Icg.NE.1 ) THEN
                     IF ( MAGA(Iexp).EQ.0 .AND. ir.NE.IPATH(n) ) GOTO 10
                     IF ( ABS(ir-IPATH(n)).GT.1 ) GOTO 10
                  ENDIF
                  ld = LDNUM(la,n)
                  IF ( ld.EQ.0 ) THEN
                     ir = ir + NSTOP(n) - NSTART(n)
                  ELSE
                     CALL LSLOOP(ir,n,nz,ld,lam,la,ssqrt,Icg,Iexp)
                  ENDIF
                  GOTO 10
               ENDIF
            ENDIF
         ENDDO
         IF ( nz.GT.LP7 ) THEN
            WRITE (22,99001) LP7
99001       FORMAT (1x,
     &              'ERROR - NUMBER OF ELEMENTS IN ZETA ARRAY EXCEEDS',
     &              'ZEMAX',5X,'(ZEMAX =',I6,')')
         ELSE
            RETURN
         ENDIF
      ELSE
         WRITE (22,99002) LP10
99002    FORMAT (' ERROR-ISMAX EXCEEDS MAGMAX',5X,'(MAGMAX =',I4,')')
      ENDIF
      ERR = .TRUE.
      RETURN
      END
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
      FUNCTION LEADF(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMDA , LAMMAX , LDNUM , LEAD , LEADF , lsum , 
     &          MULTI , N1 , n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      lsum = 0
      n3m = N3 - 1
      IF ( n3m.NE.0 ) THEN
         DO k = 1 , n3m
            lsum = lsum + MULTI(k)
         ENDDO
      ENDIF
      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m
            lsum = lsum + LDNUM(N3,k)
         ENDDO
      ENDIF
      n1m = lsum + N2
      LEADF = LEAD(2,n1m)
      END
      FUNCTION MEM(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMDA , LAMMAX , LDNUM , LEAD , MEM , msum , MULTI , 
     &          N1 , n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      msum = 0
      IF ( N3.NE.1 ) THEN
         n3m = N3 - 1
         DO k = 1 , n3m
            msum = msum + MULTI(k)
         ENDDO
      ENDIF
      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m
            msum = msum + LDNUM(N3,k)
         ENDDO
      ENDIF
      n1m = msum + 1
      n3m = n1m + LDNUM(N3,N1)
      DO k = n1m , n3m
         msum = msum + 1
         IF ( LEAD(2,k).EQ.N2 ) GOTO 100
      ENDDO
 100  MEM = msum
      END
      SUBROUTINE CMLAB(Ii,Dsig,Tetrn)
      IMPLICIT NONE
      REAL*8 a1 , a2 , ACCA , ACCUR , ared , BETAR , d2a , DIPOL , 
     &       dista , dists , Dsig , DSIGS , emax , EMMA , EN , EP , 
     &       epmin , EPS , EROOT , FIEX
      INTEGER*4 IAXS , IEXP , iflaa , Ii , IPRM , ISKIN , ISO , IZ , 
     &          IZ1 , lexp , lexp0 , lexp1 , n , NCM , NDIM , NEXPT , 
     &          NMAX , NMAX1
      REAL*8 r3 , SPIN , TASIN , tau , taup , tcmdg , tcmrad , TETACM , 
     &       Tetrn , TLBDG , tlbrad , tmxdg , TREP , VINF , XA , XA1 , 
     &       z1 , z2 , zcmdg , zcmrad
      REAL*8 zlbrad , ZPOL
      LOGICAL ERR
      COMMON /CLCOM9/ ERR
      COMMON /SECK  / ISKIN(50)
      COMMON /PRT   / IPRM(20)
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      COMMON /BREC  / BETAR(50)
      COMMON /CAUX0 / EMMA(75) , NCM
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      lexp0 = 1
      lexp1 = NEXPT
      IF ( Ii.NE.0 ) lexp0 = Ii
      IF ( Ii.NE.0 ) lexp1 = Ii
      DO lexp = lexp0 , lexp1
         iflaa = 0
         IF ( TLBDG(lexp).LT.0 ) iflaa = 1
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99001) lexp
99001       FORMAT (1X,///10X,'** EXPERIMENT',1X,1I2,1X,'**'//)
         ENDIF
         TLBDG(lexp) = ABS(TLBDG(lexp))
         a1 = XA1(lexp)
         IF ( IZ1(lexp).LT.0 ) a1 = XA
         a2 = XA
         IF ( IZ1(lexp).LT.0 ) a2 = XA1(lexp)
         z1 = DBLE(ABS(IZ1(lexp)))
         z2 = DBLE(IZ)
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( IZ1(lexp).LT.0 .AND. (Ii.EQ.0 .AND. IPRM(10).EQ.1) )
     &           WRITE (22,99002) IZ , XA , ABS(IZ1(lexp)) , XA1(lexp)
99002       FORMAT (5X,'PROJECTILE EXCITATION OF(',1I3,',',1F7.3,
     &              ') ON(',1I3,',',1F7.3,')')
            IF ( IZ1(lexp).GT.0 .AND. (Ii.EQ.0 .AND. IPRM(10).EQ.1) )
     &           WRITE (22,99003) IZ , XA , IZ1(lexp) , XA1(lexp)
99003       FORMAT (5X,'TARGET EXCITATION OF(',1I3,',',1F7.3,') BY(',
     &              1I3,',',1F7.3,')')
         ENDIF
         dists = 1.44*(a1+a2)*z1*z2/((a1**.33333+a2**.33333)*1.25+5.)/a2
         dista = 0.0719949*(1.0+a1/a2)*z1*z2/EP(lexp)
         d2a = 20.0*dista
         VINF(lexp) = 0.0463365*SQRT(EP(lexp)/a1)
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99004) EP(lexp)
     &           , VINF(lexp)
99004       FORMAT (5X,'ENERGY',1X,1F10.3,1X,'MEV',5X,'BETA',1X,1E14.6)
            IF ( EP(lexp).GT.dists .AND. (Ii.EQ.0 .AND. IPRM(10).EQ.1) )
     &           WRITE (22,99005) (EP(lexp)/dists-1.)*100.
99005       FORMAT (5X,'***** ','BE CAREFUL-ACCORDING',
     &              ' TO D.CLINE BOMBARDING ENERGY',1X,1F6.2,1X,'PC',1X,
     &              ' TOO HIGH FOR HEAD-ON COLLISIONS! *****')
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99006) d2a
99006       FORMAT (5X,
     &             'DISTANCE OF CLOSEST APPROACH FOR HEAD-ON COLLISIONS'
     &             ,1X,1F10.4,1X,'FM')
         ENDIF
         tlbrad = TLBDG(lexp)/57.2957795
         ared = 1.0 + a1/a2
         emax = EP(lexp)/ared
         DO n = 1 , NMAX
            IF ( EN(n).GT.emax ) GOTO 50
         ENDDO
         epmin = EP(lexp) - EN(NCM)*ared
         taup = SQRT(EP(lexp)/epmin)
         tau = taup*a1/a2
         IF ( tau.LE.1.0 ) GOTO 100
         tmxdg = TASIN(1.0/tau)*57.2957795
         IF ( tmxdg.GE.TLBDG(lexp) ) GOTO 100
         WRITE (22,99007) tmxdg , lexp
99007    FORMAT (1X,'ERROR- MAXIMUM SCATTERING ANGLE IS ',F7.2,
     &           ' DEGREES',' FOR EXPERIMENT ',1I2)
         GOTO 200
 50      WRITE (22,99008) emax , lexp
99008    FORMAT (1X,'ERROR- MAXIMUM EXCITATION ENERGY IS ',F8.4,' MEV',
     &           ' FOR EXPERIMENT ',1I2)
         GOTO 200
 100     tcmrad = tlbrad + TASIN(tau*SIN(tlbrad))
         tcmdg = tcmrad*57.2957795
         IF ( tau.GT.1.0 ) THEN
            IF ( IPRM(1).EQ.1 ) THEN
               IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99009)
     &              tcmdg , lexp
99009          FORMAT (5X,'SECOND POSSIBLE CM SCATTERING ANGLE IS',F7.2,
     &                 ' DEGREES FOR EXPERIMENT ',1I2)
            ENDIF
            IF ( ISKIN(lexp).NE.1 ) THEN
               tcmdg = 180. + 2.*TLBDG(lexp) - tcmdg
               tcmrad = tcmdg/57.2957795
            ENDIF
         ENDIF
         EPS(lexp) = 1./SIN(tcmrad/2.)
         TETACM(lexp) = tcmrad
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99010) tcmdg , 
     &           EPS(lexp)
99010       FORMAT (5X,'CM SCATTERING ANGLE',1X,1F10.3,1X,'DEG',5X,
     &              'EPSILON',1X,1F10.4)
         ENDIF
         IF ( IZ1(lexp).GT.0 ) BETAR(lexp) = a1*a2/(a1+a2)
     &        **2*(1.+taup*taup-2.*taup*COS(tcmrad))*epmin
         IF ( IZ1(lexp).LT.0 ) BETAR(lexp) = (a2/(a1+a2))
     &        **2*(1.+tau*tau+2.*tau*COS(tcmrad))*epmin
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99011)
     &           BETAR(lexp)
99011       FORMAT (5X,'RECOIL ENERGY(MEV)',2X,1F10.4)
         ENDIF
         BETAR(lexp) = .0463365*SQRT(BETAR(lexp)/XA)
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99012)
     &           BETAR(lexp)
99012       FORMAT (5X,'RECOIL BETA',2X,1E14.6)
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99013) EP(lexp)
     &           /(dists*.5*(1.+EPS(lexp)))
99013       FORMAT (5X,'BOMBARDING ENERGY=',1F10.3,1X,
     &              'OF SAFE BOMBARDING ENERGY AT THIS ANGLE')
         ENDIF
         IF ( iflaa.NE.1 ) THEN
            IF ( ABS(tcmdg-180.).LT.1.E-5 ) THEN
               r3 = (1.-tau)**2
            ELSE
               r3 = SIN(tlbrad)/SIN(tcmrad)
               r3 = r3*r3*ABS(COS(tcmrad-tlbrad))
               r3 = 1./r3
            ENDIF
         ENDIF
         zcmdg = 180. - tcmdg
         zcmrad = zcmdg/57.2957795
         zlbrad = ATAN(SIN(zcmrad)/(COS(zcmrad)+taup))
         IF ( iflaa.NE.0 ) THEN
            IF ( ABS(tcmdg-180.).LT.1.E-5 ) THEN
               r3 = (1.+taup)**2
               TLBDG(lexp) = 0.
            ELSE
               r3 = SIN(zlbrad)/SIN(zcmrad)
               r3 = r3*r3
               r3 = r3*ABS(COS(zcmrad-zlbrad))
               r3 = 1./r3
               TLBDG(lexp) = zlbrad*57.2955795
            ENDIF
         ENDIF
         Dsig = 250.*r3*SQRT(EP(lexp)/(EP(lexp)-ared*EN(NCM)))
     &          *dista*dista*(EPS(lexp))**4
         EROOT(lexp) = SQRT(EPS(lexp)*EPS(lexp)-1.)
         DSIGS(lexp) = Dsig
         Tetrn = zlbrad
         IF ( IZ1(lexp).LT.0. ) Tetrn = tlbrad
         TREP(lexp) = Tetrn
      ENDDO
      IPRM(10) = 0
      RETURN
 200  ERR = .TRUE.
      RETURN
      END
      SUBROUTINE QE(C,D,B2,C2,D2,B4,B6,D3,B8,C4,D4,B10,D5,B12,D6,Lmda,
     &              Pol,Cq)
      IMPLICIT NONE
      REAL*8 B10 , B12 , B2 , B4 , B6 , B8 , C , C2 , C4 , Cq , D , D2 , 
     &       D3 , D4 , D5 , D6 , Pol
      INTEGER*4 Lmda
      DIMENSION Cq(7)
      IF ( Lmda.EQ.2 ) THEN
         Cq(1) = 0.75*(2.0*C2-D2)/B4*Pol
         Cq(2) = -1.83711730*C*D/B4*Pol
         Cq(3) = -0.91855865*D2/B4*Pol
         RETURN
      ELSEIF ( Lmda.EQ.3 ) THEN
         Cq(1) = 1.875*C*(2.0*C2-3.0*D2)/B6
         Cq(2) = -1.62379763*(4.0*C2-D2)*D/B6
         Cq(3) = -5.13489890*C*D2/B6
         Cq(4) = 2.09631373*D3/B6
         RETURN
      ELSEIF ( Lmda.EQ.4 ) THEN
         Cq(1) = 1.09375000*(8.0*C4-24.0*C2*D2+3.0*D4)/B8
         Cq(2) = -4.89139867*C*(4.0*C2-3.0*D2)*D/B8
         Cq(3) = -3.45874113*(6.0*C2-D2)*D2/B8
         Cq(4) = 12.9414244*C*D3/B8
         Cq(5) = 4.57548440*D4/B8
         RETURN
      ELSEIF ( Lmda.EQ.5 ) THEN
         Cq(1) = 1.230468*C*(-14.*C2*(9.*D2+B2)+30.*B4)/B10
         Cq(2) = -1.347911*D*(35.*C2*(-3.*D2+B2)+5.*B4)/B10
         Cq(3) = -35.662372*D2*C*(-3.*D2+2.*B2)/B10
         Cq(4) = 7.279552*D3*(9.*C2-B2)/B10
         Cq(5) = 30.884521*D4*C/B10
         Cq(6) = -9.766543*D5/B10
         RETURN
      ELSEIF ( Lmda.EQ.6 ) THEN
         Cq(1) = 2.707031*(21.*C2*(-C2*(11.*D2+4.*B2)+5.*B4)-5.*B6)/B12
         Cq(2) = -17.543567*D*C*(3.*C2*(-11.*D2+B2)+5.*B4)/B12
         Cq(3) = -13.869408*D2*(3.*C2*(-11.*D2+5.*B2)+B4)/B12
         Cq(4) = -27.738815*D3*C*(-11.*D2+8.*B2)/B12
         Cq(5) = 15.193177*D4*(11.*C2-B2)/B12
         Cq(6) = -71.262308*D5*C/B12
         Cq(7) = -20.571656*D6/B12
         GOTO 99999
      ENDIF
      Cq(1) = 0.5*C/B2
      Cq(2) = -0.35355339*D/B2
      RETURN
99999 END
      SUBROUTINE QM(C,D,B2,B4,Ert,Lmda,Cq)
      IMPLICIT NONE
      REAL*8 B2 , B4 , C , Cq , D , Ert
      INTEGER*4 Lmda
      DIMENSION Cq(7)
      IF ( Lmda.EQ.8 ) THEN
         Cq(1) = -.9185586536*C*Ert/B4
         Cq(2) = -Cq(1)*D/C
         GOTO 99999
      ENDIF
      Cq(1) = -.3535533905*Ert/B2
      RETURN
99999 END
      SUBROUTINE SNAKE(Nexp,Zpol)
      IMPLICIT NONE
      REAL*8 b10 , b12 , b2 , b4 , b6 , b8 , c , c2 , c4 , c6 , CH , 
     &       chi , cq , d , d2 , d3 , d4 , d5 , d6 , EPS
      REAL*8 EROOT , ert , FIEX , pol , SH , shi , ZETA , Zpol
      INTEGER*4 IAXS , ibm , icm , icnt , idm , IEXP , irl , j , k , 
     &          lloc , lmd , lmda , LOCQ , LP1 , LP10 , LP11 , LP12 , 
     &          LP13 , LP14 , LP2
      INTEGER*4 LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , LZETA , mimx , 
     &          Nexp , nind , nlm
      DIMENSION lloc(8) , cq(7) , irl(8)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /HIPER / SH(365) , CH(365)
      icnt = 0
 100  icnt = icnt + 1
      CALL QRANGE(icnt,nlm,lloc,ibm,icm,idm,irl)
      IF ( nlm.EQ.0 ) RETURN
      chi = CH(icnt)
      shi = SH(icnt)
      b2 = EPS(Nexp)*chi + 1.
      pol = 1. - Zpol/b2
      b2 = b2*b2
      IF ( ibm.NE.2 ) THEN
         b4 = b2*b2
         IF ( ibm.NE.4 ) THEN
            b6 = b4*b2
            IF ( ibm.NE.6 ) THEN
               b8 = b4*b4
               IF ( ibm.NE.8 ) THEN
                  b10 = b6*b4
                  IF ( ibm.NE.10 ) b12 = b6*b6
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF ( icm.NE.0 ) THEN
         c = chi + EPS(Nexp)
         IF ( icm.NE.1 ) THEN
            c2 = c*c
            IF ( icm.NE.2 ) THEN
               c4 = c2*c2
               IF ( icm.NE.4 ) c6 = c2*c4
            ENDIF
         ENDIF
      ENDIF
      IF ( idm.NE.0 ) THEN
         d = EROOT(Nexp)*shi
         IF ( idm.NE.1 ) THEN
            d2 = d*d
            IF ( idm.NE.2 ) THEN
               d3 = d*d2
               IF ( idm.NE.3 ) THEN
                  d4 = d2*d2
                  IF ( idm.NE.4 ) THEN
                     d5 = d3*d2
                     IF ( idm.NE.5 ) d6 = d3*d3
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      DO j = 1 , nlm
         lmda = lloc(j)
         IF ( lmda.GT.6 ) THEN
            lmd = lmda
            lmda = lmda - 6
            ert = EROOT(Nexp)
            CALL QM(c,d,b2,b4,ert,lmda,cq)
            mimx = lmda
            DO k = 1 , mimx
               nind = LOCQ(lmd,k) + icnt
               ZETA(nind+LP7) = cq(k)
            ENDDO
         ELSE
            CALL QE(c,d,b2,c2,d2,b4,b6,d3,b8,c4,d4,b10,d5,b12,d6,lmda,
     &              pol,cq)
            mimx = lmda + 1
            DO k = 1 , mimx
               nind = LOCQ(lmda,k) + icnt
               ZETA(nind+LP7) = cq(k)
            ENDDO
         ENDIF
      ENDDO
      GOTO 100
      END
      SUBROUTINE FHIP
      IMPLICIT NONE
      REAL*8 CH , er , ex , SH , w
      INTEGER*4 j , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , 
     &          LP4 , LP6 , LP7 , LP8 , LP9
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /HIPER / SH(365) , CH(365)
      w = -.03
      DO j = 1 , LP12
         w = w + .03
         ex = EXP(w)
         er = 1./ex
         SH(j) = (ex-er)/2.
         CH(j) = (ex+er)/2.
      ENDDO
      END
      SUBROUTINE ALLOC(Accur)
      IMPLICIT NONE
      REAL*8 Accur , u , v
      INTEGER*4 iflag , IRA , j , k , k1 , load , LOCQ , LP1 , LP10 , 
     &          LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , LP6 , 
     &          LP7 , LP8 , LP9
      INTEGER*4 MAXLA
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      CALL RANGEL(Accur)
      load = 0
      iflag = 0
      DO j = 1 , 8
         DO k = 1 , 7
            LOCQ(j,k) = 0
         ENDDO
      ENDDO
      DO k = 1 , 6
         k1 = k + 1
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k)
         ENDDO
      ENDDO
      DO k = 7 , 8
         k1 = k - 6
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k)
         ENDDO
      ENDDO
      IF ( load.LE.LP14 ) RETURN
      WRITE (22,99001)
99001 FORMAT (5X,'NO SPACE FOR Q FUNCTIONS TABULATION'//5X,
     &        'SORRY,JOB WILL BE BRUTALLY TERMINATED!')
      v = -1.
      u = LOG10(v)
      u = SIN(u)
      END
      SUBROUTINE RANGEL(Acc1)
      IMPLICIT NONE
      REAL*8 Acc1 , ACC50 , acl , w
      INTEGER*4 i , IRA , LAMDA , LAMMAX , LDNUM , LEAD , MAXLA , MULTI
      COMMON /A50   / ACC50
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      acl = -LOG(Acc1)
      ACC50 = Acc1/50.
      DO i = 1 , 8
         IF ( MULTI(i).NE.0 ) THEN
            IF ( i.EQ.2 .OR. i.EQ.7 ) THEN
               w = acl/2. + .203
            ELSEIF ( i.EQ.3 .OR. i.EQ.8 ) THEN
               w = acl/3. + .536
            ELSEIF ( i.EQ.4 ) THEN
               w = acl/4. + .716
            ELSEIF ( i.EQ.5 ) THEN
               w = acl/5. + .829
            ELSEIF ( i.EQ.6 ) THEN
               w = acl/6. + .962
            ELSE
               w = acl - .693
            ENDIF
            w = w/.03
            IRA(i) = INT(w+1.5)
         ELSE
            IRA(i) = 0
         ENDIF
      ENDDO
      IF ( IRA(7).NE.0 ) IRA(7) = IRA(7) + 1
      IF ( IRA(8).NE.0 ) IRA(8) = IRA(8) + 1
      END
      SUBROUTINE QRANGE(Icnt,Nlm,Lloc,Ibm,Icm,Idm,Irl)
      IMPLICIT NONE
      INTEGER*4 Ibm , Icm , Icnt , Idm , IRA , Irl , is , k , ke , km , 
     &          l , LAMDA , LAMMAX , ld , LDNUM , LEAD , Lloc , ls , 
     &          MAXLA , MULTI
      INTEGER*4 nlend , Nlm
      DIMENSION Lloc(8) , Irl(8)
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      IF ( Icnt.EQ.1 ) THEN
         Nlm = 0
         DO l = 1 , 8
            Lloc(l) = 0
            Irl(l) = 0
         ENDDO
         DO k = 1 , 6
            ke = 7 - k
            km = 13 - k
            IF ( km.LE.8 ) THEN
               IF ( MULTI(km).NE.0 ) THEN
                  Nlm = Nlm + 1
                  Lloc(Nlm) = km
                  Irl(Nlm) = IRA(km)
               ENDIF
            ENDIF
            IF ( MULTI(ke).NE.0 ) THEN
               Nlm = Nlm + 1
               Lloc(Nlm) = ke
               Irl(Nlm) = IRA(ke)
            ENDIF
         ENDDO
         nlend = INT((DBLE(Nlm)+1.1)/2.)
         DO k = 1 , nlend
            ke = Nlm - k + 1
            ls = Lloc(ke)
            is = Irl(ke)
            Lloc(ke) = Lloc(k)
            Irl(ke) = Irl(k)
            Lloc(k) = ls
            Irl(k) = is
         ENDDO
         l = 0
         DO k = 1 , 6
            IF ( MULTI(k).NE.0 ) l = k
         ENDDO
         Icm = MIN(4,l)
         Ibm = 2*l
         Idm = l
         l = 0
         DO k = 7 , 8
            ke = k - 6
            IF ( MULTI(k).NE.0 ) l = ke
         ENDDO
         Ibm = MAX(Ibm,2*l)
         Idm = MAX(Idm,l)
         IF ( Icm.EQ.1 .AND. l.GT.1 ) Icm = 2
         MAXLA = Lloc(1)
         RETURN
      ELSE
         IF ( Irl(Nlm).GE.Icnt ) RETURN
         ld = Lloc(Nlm)
         Lloc(Nlm) = 0
         Nlm = Nlm - 1
         IF ( Nlm.EQ.0 ) RETURN
         IF ( ld.GT.6 ) RETURN
         l = Lloc(Nlm)
         IF ( l.GT.6 ) l = l - 6
         Icm = MIN(2,l)
         Ibm = 2*l
         Idm = l
      ENDIF
      END
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
               IF ( NSTART(nhold).NE.0 ) GOTO 20
               GOTO 30
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
                        IF ( n.LE.NMAX ) GOTO 42
                        GOTO 100
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
      SUBROUTINE LAISUM(Ir,N,Rsg,Lam,Ld,Nz,I57)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , CAT , D2W , DIPOL , ELM , ELML , ELMU , EN , 
     &       q , rmir , rmis , rmu , Rsg , SA , SPIN , z , ZETA , ZPOL
      INTEGER*4 i2 , i3 , I57 , iii , indq , indx , Ir , irs , is , 
     &          is1 , is2 , ISG , ISG1 , ISHA , ISMAX , ismin , ISO , 
     &          isplus , ISSTAR , ISSTO
      INTEGER*4 KDIV , la , Lam , LAMR , Ld , LOCQ , LP1 , LP10 , LP11 , 
     &          LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , 
     &          LP9 , LZETA
      INTEGER*4 m , mrange , MSTORE , mua , N , NDIV , NPT , NSTART , 
     &          NSTOP , NSW , Nz
      COMPLEX*16 ARM , FAZA , pamp , EXPO , pamp1
      COMMON /PSPIN / ISHA
      COMMON /AZ    / ARM(600,7)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      COMMON /ADBXI / EXPO(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      rmir = CAT(Ir,3)
      iii = 0
      IF ( Lam.GT.6 ) iii = 1
      la = Lam
      IF ( Lam.GT.6 ) Lam = Lam - 6
      DO i2 = 1 , Ld
         pamp = (0.,0.)
         m = MSTORE(1,i2)
         indx = MSTORE(2,i2)
         ismin = 0
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
                  IF ( ISO.NE.0 .OR. rmir.LE..1 .OR. rmis.LE..1 ) THEN
                     rmu = rmis - rmir
                     mua = ABS(rmu) + 1.1
                     IF ( la.LE.6 .OR. mua.NE.1 ) THEN
                        indq = LOCQ(Lam,mua) + NPT
                        Nz = Nz + 1
                        z = ZETA(Nz)
                        q = ZETA(indq+LP7)
                        IF ( NDIV.NE.0 ) q = ZETA(indq+LP7) + DBLE(KDIV)
     &                       *(ZETA(indq+LP7+ISG1)-ZETA(indq+LP7))
     &                       /DBLE(NDIV)
                        pamp1 = FAZA(la,mua,rmu,Rsg)*q*z
                        IF ( ISO.NE.0 .OR. rmir.LE..1 ) THEN
                           pamp = pamp1*ARM(is,I57) + pamp
                           IF ( ISO.EQ.0 .AND. rmis.GT..1 ) GOTO 10
                        ENDIF
                        IF ( N.NE.m ) THEN
                           irs = (-1)**(INT(rmir+rmis)-ISHA+iii)
                           ARM(is,6) = ARM(is,6) + irs*pamp1*ARM(Ir,I57)
                           ISSTAR(i2) = MIN(is,ISSTAR(i2))
                           ISSTO(i2) = MAX(is,ISSTO(i2))
                        ENDIF
                     ENDIF
                  ENDIF
 10            ENDDO
               IF ( N.EQ.m ) THEN
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)
               ELSE
                  ARM(Ir,4) = ARM(Ir,4) + pamp*ELM(indx)*EXPO(indx)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      Lam = la
      END
      COMPLEX*16 FUNCTION EXPON(Inx,Npt,Isg,Isg1,Ndiv,Kdiv)
      IMPLICIT NONE
      REAL*8 ADB , XI
      INTEGER*4 Inx , Isg , Isg1 , Kdiv , Ndiv , Npt
      COMPLEX*16 expo1 , ci , expox , TCEXP
      COMMON /ADX   / ADB(365)
      COMMON /CXI   / XI(500)
      DATA ci/(0.,1.)/
      expox = TCEXP(ci*XI(Inx)*ADB(Npt)*Isg)
      EXPON = expox
      IF ( Ndiv.NE.0 ) THEN
         expo1 = TCEXP(ci*XI(Inx)*ADB(Npt+Isg1)*Isg)
         EXPON = expox + DBLE(Kdiv)*(expo1-expox)/DBLE(Ndiv)
      ENDIF
      END
      COMPLEX*16 FUNCTION FAZA(La,Mi,Rmu,Rsg)
      IMPLICIT NONE
      INTEGER*4 ieven , La , Mi
      REAL*8 Rmu , Rsg
      COMPLEX*16 ci
      DATA ci/(0.,1.)/
      IF ( La.GT.6 ) THEN
         FAZA = -ci
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         IF ( La.EQ.7 ) RETURN
         IF ( Mi.EQ.2 ) RETURN
         FAZA = CMPLX(Rsg,0.)
         IF ( Rmu.LT.0. ) FAZA = -FAZA
         GOTO 99999
      ELSE
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            FAZA = -ci
            RETURN
         ENDIF
      ENDIF
      FAZA = CMPLX(Rsg,0.)
      RETURN
99999 END
      SUBROUTINE SETIN
      IMPLICIT NONE
      REAL*8 ADB , CH , EPS , EROOT , FIEX , SH
      INTEGER*4 IAXS , IEXP , k , LP1 , LP10 , LP11 , LP12 , LP13 , 
     &          LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /HIPER / SH(365) , CH(365)
      COMMON /ADX   / ADB(365)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      DO k = 1 , LP12
         ADB(k) = EPS(IEXP)*SH(k) + .03*(k-1)
      ENDDO
      END
      SUBROUTINE STING(Irld)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELML , ELMU , rsg , SA , w0 , ZETA
      INTEGER*4 i , i57 , ibg , iend , IFLG , indx , IRA , Irld , is2 , 
     &          ISG , ISG1 , ISMAX , ISSTAR , ISSTO , j , j1 , jj , 
     &          KDIV , lam , LAMDA
      INTEGER*4 LAMMAX , LAMR , ld , LDNUM , LEAD , LZETA , maxh , 
     &          MAXLA , mm , MSTORE , MULTI , n , NDIV , NPT , NSW , nz
      COMPLEX*16 ARM , EXPO
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /AZ    / ARM(600,7)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /ADBXI / EXPO(500)
      COMMON /FLA   / IFLG
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /RNG   / IRA(8) , MAXLA
      maxh = MAXLA
 100  ISG = -1
      n = 1
      rsg = -1.
      IFLG = 1
      w0 = IRA(MAXLA)*.03 + .03
      DO j = 1 , ISMAX
         DO jj = 1 , 6
            ARM(j,jj) = (0.,0.)
         ENDDO
      ENDDO
      ARM(Irld,5) = (1.,0.)
      DO j = 1 , 8
         LAMR(j) = 0
      ENDDO
      LAMR(MAXLA) = 1
      NPT = IRA(MAXLA) + 1
      IF ( MAXLA.EQ.7 .AND. IRA(2).NE.0 ) THEN
         LAMR(2) = 1
         NPT = NPT - 1
         w0 = w0 - .03
      ENDIF
      NDIV = 0
      KDIV = 0
      DO j = 1 , 4
         NPT = NPT - 1
         DO j1 = 1 , LAMMAX
            lam = LAMDA(j1)
            IF ( LAMR(lam).NE.0 ) THEN
               CALL NEWLV(n,ld,lam)
               IF ( ld.NE.0 ) THEN
                  nz = LZETA(lam)
                  ld = LDNUM(lam,1)
                  i57 = 5
                  CALL LAISUM(Irld,n,rsg,lam,ld,nz,i57)
                  DO mm = 1 , ld
                     indx = MSTORE(2,mm)
                     ibg = ISSTAR(mm)
                     iend = ISSTO(mm)
                     DO is2 = ibg , iend
                        ARM(is2,4) = ARM(is2,4) + ARM(is2,6)*ELM(indx)
     &                               /EXPO(indx)
                        ARM(is2,6) = (0.,0.)
                     ENDDO
                  ENDDO
               ELSEIF ( j1.EQ.MAXLA ) THEN
                  IRA(MAXLA) = -IRA(MAXLA)
                  DO jj = 1 , LAMMAX
                     lam = LAMDA(jj)
                     IF ( IRA(lam).GT.0 ) GOTO 105
                  ENDDO
 105              MAXLA = LAMDA(jj)
                  GOTO 100
               ENDIF
            ENDIF
         ENDDO
         IF ( j.EQ.4 ) GOTO 200
         DO i = 1 , ISMAX
            ARM(i,j) = ARM(i,4)
            ARM(i,4) = (0.,0.)
         ENDDO
      ENDDO
 200  CALL LAIAMP(Irld,w0)
      MAXLA = maxh
      DO jj = 1 , 8
         IRA(jj) = ABS(IRA(jj))
      ENDDO
      END
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
      rmir = CAT(Ir,3)
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
                           rmis = CAT(is,3)
                           rmu = rmis - rmir
                           mua = ABS(rmu) + 1.1
                           IF ( lam.LE.6 .OR. mua.NE.1 ) THEN
                              CALL FAZA1(la,mua,rmir,rmis,dis,rmu)
                              pm = ELM(indx)*z
                              uhuj = STAMP(epsi,errt,xiv,.03,W0,lam,mua)
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
      SUBROUTINE FAZA1(La,Mi,Rmir,Rmis,Dis,Rmu)
      IMPLICIT NONE
      INTEGER*4 ieven , irs , La , Mi
      REAL*8 Rmir , Rmis , Rmu
      COMPLEX*16 Dis , ci
      DATA ci/(0.,1.)/
      irs = (-1)**INT(Rmir+Rmis)
      IF ( La.EQ.7 ) THEN
         Dis = -ci*irs
         IF ( Rmu.LT.0. ) Dis = -Dis
         GOTO 99999
      ELSE
         ieven = (-1)**Mi
         IF ( ieven.LE.0 ) THEN
            Dis = -ci*irs
            RETURN
         ENDIF
      ENDIF
      Dis = CMPLX(-DBLE(irs),0.)
      RETURN
99999 END
      SUBROUTINE TRINT(Arg,Si,Ci)
      IMPLICIT NONE
      REAL*8 a , Arg , c , Ci , f , g , POL4 , s , Si
      a = Arg*Arg
      IF ( Arg.LT.1. ) THEN
         Si = POL4(0.,2.83446712E-5,-1.66666667E-3,.055555555,-1.,a)
         Si = Si*Arg
         Ci = POL4(-3.100198413E-6,2.314814815E-4,-.0104166667,.25,0.,a)
         Ci = Ci - LOG(Arg)
         GOTO 99999
      ENDIF
      s = SIN(Arg)
      c = COS(Arg)
      f = POL4(1.,38.027246,265.187033,335.67732,38.102495,a)
      f = f/POL4(1.,40.021433,322.624911,570.23628,157.105423,a)/Arg
      g = POL4(1.,42.242855,302.757865,352.018498,21.821899,a)
      g = g/POL4(1.,48.196927,482.485984,1114.978885,449.690326,a)/a
      Si = f*c + g*s
      Ci = g*c - f*s
      RETURN
99999 END
      FUNCTION POL4(C0,C1,C2,C3,C4,A)
      IMPLICIT NONE
      REAL*8 A , C0 , C1 , C2 , C3 , C4 , POL4
      POL4 = 1.
      IF ( ABS(A).GT.1.E+9 ) RETURN
      POL4 = C4 + A*(C3+A*(C2+A*(C1+A*C0)))
      END
      COMPLEX*16 FUNCTION STAMP(Epsi,Errt,Xiv,Dw,W0,Lmda,Mua)
      IMPLICIT NONE
      REAL*8 a , axi , b , bic , bic2 , bis , bis2 , ca , cb , cia , 
     &       cib , cic , cis , Dw , dwi , Epsi , Errt , ex , exa , fct
      INTEGER*4 la , Lmda , mi , Mua
      REAL*8 sa , sb , sia , sib , W0 , Xiv
      mi = Mua - 1
      axi = ABS(Xiv)
      la = Lmda
      IF ( Lmda.EQ.7 ) la = 3
      IF ( axi.LT.1.E-5 ) THEN
         a = -2.*W0
         IF ( la.EQ.3 ) a = -W0
         exa = EXP(a)
         dwi = 3*Dw
         cic = exa*(EXP(dwi)-1.)
         STAMP = CMPLX(cic,0.)
         IF ( la.EQ.2 ) THEN
            IF ( mi.EQ.0 ) fct = 3.*(3.-Epsi*Epsi)/Epsi/Epsi/Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.837117307*Errt/Epsi/Epsi/Epsi/Epsi
            IF ( mi.EQ.2 ) fct = -3.674234613*Errt*Errt/Epsi/Epsi/Epsi/
     &                           Epsi
         ELSEIF ( la.EQ.3 ) THEN
            fct = -1.414213562*Errt/Epsi/Epsi
         ELSE
            IF ( mi.EQ.0 ) fct = 1./Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.414213562*Errt/Epsi/Epsi
         ENDIF
      ELSE
         ex = EXP(W0)/2.
         b = axi*(Epsi*ex+W0)
         CALL TRINT(b,sib,cib)
         sb = SIN(b)/b
         cb = COS(b)/b
         bis = sb + cib
         bic = cb - sib
         bis2 = -sb/b
         bic2 = -cb/b
         dwi = -3.*Dw
         exa = EXP(dwi)
         a = axi*(Epsi*ex*exa+W0+dwi)
         sa = SIN(a)/a
         ca = COS(a)/a
         CALL TRINT(a,sia,cia)
         cis = sa + cia - bis
         cic = ca - sia - bic
         IF ( la.EQ.1 ) THEN
            STAMP = CMPLX(cic,cis)
         ELSE
            dwi = (bic2-cis+ca/a)/2.
            exa = (bis2+cic+sa/a)/2.
            STAMP = CMPLX(dwi,exa)
         ENDIF
         IF ( la.EQ.2 ) THEN
            IF ( mi.EQ.0 ) fct = .75*(3.-Epsi*Epsi)*axi*axi/Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.837117307*Errt*axi*axi/Epsi/Epsi
            IF ( mi.EQ.2 ) fct = -.9185586535*Errt*Errt*axi*axi/Epsi/
     &                           Epsi
         ELSEIF ( la.EQ.3 ) THEN
            fct = -.3535533905*Errt*axi*axi
         ELSE
            IF ( mi.EQ.0 ) fct = .5*axi/Epsi
            IF ( mi.EQ.1 ) fct = .3535533907*Errt*axi/Epsi
         ENDIF
      ENDIF
      STAMP = STAMP*fct
      STAMP = CONJG(STAMP)
      END
      SUBROUTINE RESET(Iso)
      IMPLICIT NONE
      REAL*8 CAT
      INTEGER*4 ir , ISMAX , Iso , j , NDIM , NMAX , NMAX1 , NSTART , 
     &          NSTOP
      COMPLEX*16 ARM
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /AZ    / ARM(600,7)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX
            ir = NSTART(j) - 1
 20         ir = ir + 1
            ARM(ir,1) = ARM(ir,2)
            ARM(ir,2) = ARM(ir,3)
            ARM(ir,3) = ARM(ir,4)
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
      DO j = 1 , ISMAX
         ARM(j,1) = ARM(j,2)
         ARM(j,2) = ARM(j,3)
         ARM(j,3) = ARM(j,4)
      ENDDO
      RETURN
99999 END
      SUBROUTINE HALF(Iso)
      IMPLICIT NONE
      REAL*8 CAT
      INTEGER*4 ir , ISMAX , Iso , j , NDIM , NMAX , NMAX1 , NSTART , 
     &          NSTOP
      COMPLEX*16 ARM , fpom
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /AZ    / ARM(600,7)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX
            ir = NSTART(j) - 1
 20         ir = ir + 1
            fpom = ARM(ir,3)
            ARM(ir,1) = -.0625*(ARM(ir,1)+ARM(ir,4))
     &                  + .5625*(ARM(ir,2)+ARM(ir,3))
            ARM(ir,3) = ARM(ir,3)*.75 + .375*ARM(ir,4) - ARM(ir,2)/8.
            ARM(ir,2) = fpom
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
      DO j = 1 , ISMAX
         fpom = ARM(j,3)
         ARM(j,1) = -.0625*(ARM(j,4)+ARM(j,1))
     &              + .5625*(ARM(j,2)+ARM(j,3))
         ARM(j,3) = ARM(j,3)*.75 + .375*ARM(j,4) - ARM(j,2)/8.
         ARM(j,2) = fpom
      ENDDO
      RETURN
99999 END
      SUBROUTINE DOUBLE(Iso)
      IMPLICIT NONE
      REAL*8 CAT
      INTEGER*4 ir , ISMAX , Iso , j , NDIM , NMAX , NMAX1 , NSTART , 
     &          NSTOP
      COMPLEX*16 ARM , fpom
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /AZ    / ARM(600,7)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX
            ir = NSTART(j) - 1
 20         ir = ir + 1
            fpom = ARM(ir,2)
            ARM(ir,2) = -8.*ARM(ir,3) + 6.*ARM(ir,2) + 3.*ARM(ir,4)
            ARM(ir,1) = -16.*ARM(ir,1) + 9.*ARM(ir,2) + 9.*fpom - 
     &                  ARM(ir,4)
            ARM(ir,3) = fpom
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
      DO j = 1 , ISMAX
         fpom = ARM(j,2)
         ARM(j,2) = -8.*ARM(j,3) + 6.*ARM(j,2) + 3.*ARM(j,4)
         ARM(j,1) = -16.*ARM(j,1) + 9.*ARM(j,2) + 9.*fpom - ARM(j,4)
         ARM(j,3) = fpom
      ENDDO
      RETURN
99999 END
      SUBROUTINE PATH(Irld)
      IMPLICIT NONE
      REAL*8 CAT , spm , vl
      INTEGER*4 i , IPATH , Irld , ISMAX , isp , ist , j , MAGA , NDIM , 
     &          NMAX , NMAX1 , NSTART , NSTOP
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      spm = CAT(Irld,3)
      DO i = 2 , NMAX
         IPATH(i) = 0
         ist = NSTART(i)
         IF ( ist.NE.0 ) THEN
            isp = NSTOP(i)
            DO j = ist , isp
               vl = CAT(j,3)
               IF ( ABS(vl-spm).LT.1.E-6 ) GOTO 50
            ENDDO
         ENDIF
         GOTO 100
 50      IPATH(i) = j
 100  ENDDO
      IPATH(1) = Irld
      END
      SUBROUTINE INTG(Ien)
      IMPLICIT NONE
      REAL*8 ACC50 , ACCA , ACCUR , CAT , D2W , DIPOL , EN , f , rim , 
     &       rl , SPIN , srt , ZPOL
      INTEGER*4 i , i57 , Ien , IFAC , IFLG , ihold , intend , INTERV , 
     &          IPATH , ir , ir1 , IRA , ISG , ISG1 , ISMAX , ISO , k , 
     &          kast , KDIV , LAMR
      INTEGER*4 MAGA , MAXLA , mir , n , NDIM , NDIV , NMAX , NMAX1 , 
     &          NPT , NSTART , NSTOP , NSW
      COMPLEX*16 ARM , hold
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /AZ    / ARM(600,7)
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /A50   / ACC50
      COMMON /CLCOM0/ IFAC(75)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /FLA   / IFLG
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /CEXC9 / INTERV(50)
      intend = INTERV(Ien)
      D2W = .03
      NSW = 1
      kast = 0
      NDIV = 0
      KDIV = 0
 100  IF ( (NPT+NSW).GT.IRA(MAXLA) .AND. ISG.GT.0 ) RETURN
      DO i = 1 , 8
         LAMR(i) = 0
         IF ( (NPT+NSW).LT.IRA(i) ) LAMR(i) = 1
      ENDDO
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX
            ir = NSTART(n) - 1
 120        ir = ir + 1
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
            mir = CAT(ir,3)
            ir1 = ir - 2*mir
            ARM(ir1,7) = IFAC(n)*ARM(ir,7)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 120
         ENDDO
      ELSE
         DO ir = 1 , ISMAX
            ARM(ir,7) = ARM(ir,5)
     &                  + D2W/24.*(55.0*ARM(ir,4)-59.0*ARM(ir,3)
     &                  +37.0*ARM(ir,2)-9.0*ARM(ir,1))
         ENDDO
      ENDIF
      NPT = NPT + NSW*ISG
      IF ( NPT.GT.0 ) THEN
         IF ( NDIV.EQ.0 ) GOTO 200
         KDIV = KDIV + 1
         IF ( KDIV.LT.NDIV ) GOTO 200
         KDIV = 0
         NPT = NPT + ISG
         IF ( NPT.GT.0 ) GOTO 200
      ENDIF
      NPT = -NPT + 2
      ISG = 1
 200  CALL RESET(ISO)
      IFLG = 1
      i57 = 7
      CALL AMPDER(i57)
      IF ( ISO.EQ.0 ) THEN
         DO n = 1 , NMAX
            ir = NSTART(n) - 1
 220        ir = ir + 1
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
            mir = CAT(ir,3)
            ir1 = ir - 2*mir
            ARM(ir1,5) = IFAC(n)*ARM(ir,5)
            IF ( DBLE(mir).LT.-0.1 ) GOTO 220
         ENDDO
      ELSE
         DO ir = 1 , ISMAX
            ARM(ir,5) = ARM(ir,5)
     &                  + D2W/24.*(9.0*ARM(ir,4)+19.0*ARM(ir,3)
     &                  -5.0*ARM(ir,2)+ARM(ir,1))
         ENDDO
      ENDIF
      kast = kast + 1
      IFLG = 0
      i57 = 5
      CALL AMPDER(i57)
      IF ( (LAMR(2)+LAMR(3)).NE.0 ) THEN
         IF ( kast.GE.intend ) THEN
            kast = 0
            f = 0.
            DO k = 1 , NMAX
               ihold = IPATH(k)
               IF ( ihold.NE.0 ) THEN
                  hold = ARM(ihold,5) - ARM(ihold,7)
                  rl = DBLE(hold)
                  rim = IMAG(hold)
                  srt = rl*rl + rim*rim
                  f = MAX(f,srt)
               ENDIF
            ENDDO
            f = SQRT(f)/14.
            IF ( f.GT.ACCUR .OR. f.LT.ACC50 ) THEN
               IF ( f.LT.ACC50 ) THEN
                  CALL DOUBLE(ISO)
                  D2W = 2.*D2W
                  NSW = 2*NSW
                  intend = (DBLE(intend)+.01)/2.
                  IF ( intend.EQ.0 ) intend = 1
                  IF ( NSW.LT.1 ) THEN
                     NDIV = (DBLE(NDIV)+.01)/2.
                     IF ( NDIV.LT.2 ) THEN
                        NDIV = 0
                        NSW = 1
                     ENDIF
                  ENDIF
               ELSE
                  CALL HALF(ISO)
                  D2W = D2W/2.
                  NSW = (DBLE(NSW)+.01)/2.
                  intend = 2*intend
                  IF ( NSW.LT.1 ) THEN
                     NDIV = 2*NDIV
                     IF ( NDIV.EQ.0 ) NDIV = 2
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      GOTO 100
      END
      SUBROUTINE NEWLV(N,Ld,La)
      IMPLICIT NONE
      REAL*8 D2W
      INTEGER*4 i2 , IFLG , indx , ISG , ISG1 , ISSTAR , ISSTO , KDIV , 
     &          La , LAMDA , LAMMAX , LAMR , Ld , LDNUM , LEAD , LEADF , 
     &          m , MEM , MSTORE , MULTI
      INTEGER*4 N , NDIV , NPT , NSTART , NSTOP , NSW
      COMPLEX*16 EXPO , EXPON
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      COMMON /ADBXI / EXPO(500)
      COMMON /FLA   / IFLG
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      Ld = LDNUM(La,N)
      IF ( Ld.EQ.0 ) RETURN
      DO i2 = 1 , Ld
         m = LEADF(N,i2,La)
         ISSTAR(i2) = NSTOP(m)
         ISSTO(i2) = NSTART(m)
         MSTORE(1,i2) = m
         indx = MEM(N,m,La)
         MSTORE(2,i2) = indx
         IF ( IFLG.NE.0 ) THEN
            IF ( m.NE.N ) EXPO(indx) = EXPON(indx,NPT,ISG,ISG1,NDIV,KDIV
     &                                 )
         ENDIF
      ENDDO
      END
      SUBROUTINE CODE7(Ir,Is,N,Mt,Inqa,Indx)
      IMPLICIT NONE
      INTEGER*4 IAPR , idm , idn , Indx , Inqa , IPATH , Ir , Is , 
     &          ISEX , ism , MAGA , Mt , N
      REAL*8 QAPR
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      IAPR(Indx,1) = N
      IAPR(Indx,2) = Mt
      IF ( IPATH(N).EQ.0 .OR. IPATH(Mt).EQ.0 ) THEN
         Inqa = -1
         GOTO 99999
      ELSE
         idn = Ir - IPATH(N)
         idm = Is - IPATH(Mt)
         ism = idn + idm + 3
         IF ( ism.EQ.2 ) THEN
            Inqa = 2
            IF ( idn.GT.idm ) Inqa = 3
            RETURN
         ELSEIF ( ism.EQ.3 ) THEN
            Inqa = 4
            RETURN
         ELSEIF ( ism.EQ.4 ) THEN
            Inqa = 5
            IF ( idn.GT.idm ) Inqa = 6
            RETURN
         ELSEIF ( ism.NE.5 ) THEN
            Inqa = 1
            RETURN
         ENDIF
      ENDIF
      Inqa = 7
      RETURN
99999 END
      SUBROUTINE APRAM(Iexp,Inc,Indx,Irld,Acca)
      IMPLICIT NONE
      REAL*8 Acca , accah , ELM , ELML , ELMU , QAPR , SA , uwa
      INTEGER*4 i1 , i56 , i7 , IAPR , IDIVE , Iexp , img , Inc , Indx , 
     &          IPATH , Irld , ISEX , itm , IVAR , j , jidim , jj , k , 
     &          ktoto , l
      INTEGER*4 l1 , l2 , l3 , LERF , LMAXE , m , MAGA , MAGEXC , 
     &          MEMAX , MEMX6
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(600,7)
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /APRX  / LERF , IDIVE(50,2)
      LERF = 0
      accah = Acca
 100  i7 = 7
      itm = -1
      img = 3
      i1 = 1
      IF ( MAGA(Iexp).EQ.0 ) THEN
         i7 = 4
         i1 = 4
         img = 1
      ENDIF
      IF ( Inc.EQ.0 ) GOTO 300
      IF ( LERF.EQ.0 ) CALL NEWCAT(Iexp,jidim)
      IF ( LERF.EQ.0 ) CALL PODZIEL(3,Iexp)
      i56 = 5
      DO k = 1 , jidim
         ARM(k,2) = (0.,0.)
         ARM(k,5) = (0.,0.)
      ENDDO
      ARM(Irld+1,5) = (1.,0.)
 200  ktoto = 0
      LERF = 0
      l1 = IDIVE(Iexp,1)
      DO l3 = 1 , l1
         Acca = accah*l3/l1
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim)
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(1,Iexp)
            GOTO 100
         ENDIF
      ENDDO
      l2 = IDIVE(Iexp,2)
      DO l3 = 1 , l2
         Acca = accah + accah*l3/l2
         CALL POMNOZ(Acca,2,i56,ktoto,img,jidim)
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(2,Iexp)
            GOTO 100
         ENDIF
      ENDDO
      DO l = 1 , MEMX6
         DO m = i1 , i7
            QAPR(l,1,m) = -QAPR(l,1,m)
         ENDDO
      ENDDO
      DO l3 = 1 , l1
         Acca = accah*2. + accah*l3/l1
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim)
      ENDDO
      Acca = accah
      DO l = 1 , MEMX6
         DO m = i1 , i7
            QAPR(l,1,m) = -QAPR(l,1,m)
         ENDDO
      ENDDO
      IF ( Inc.NE.0 .OR. itm.NE.0 ) THEN
         IF ( Inc.EQ.0 ) THEN
            DO l = 1 , jidim
               ARM(l,6) = ARM(l,6) - ARM(l,7)
               ARM(l,6) = 50.*ARM(l,6)/ELM(Indx)
            ENDDO
            DO l = 1 , 2
               DO j = i1 , i7
                  QAPR(Indx,l,j) = QAPR(Indx,l,j)/.99
               ENDDO
            ENDDO
            DO jj = 2 , jidim
               ARM(jj-1,6) = ARM(jj,6)
            ENDDO
            GOTO 99999
         ELSE
            DO jj = 2 , jidim
               ARM(jj-1,5) = ARM(jj,5)
            ENDDO
            RETURN
         ENDIF
      ENDIF
 300  itm = itm + 1
      i56 = itm + 6
      DO k = 1 , jidim
         ARM(k,i56) = (0.,0.)
      ENDDO
      ARM(Irld+1,i56) = (1.,0.)
      uwa = -itm*.0298019802 + 1.01
      DO l = 1 , 2
         DO j = i1 , i7
            QAPR(Indx,l,j) = QAPR(Indx,l,j)*uwa
         ENDDO
      ENDDO
      DO j = 1 , jidim
         ARM(j,2) = (0.,0.)
      ENDDO
      GOTO 200
99999 END
      SUBROUTINE NEWCAT(Iexp,Jidim)
      IMPLICIT NONE
      REAL*8 a , b , FXIS1 , FXIS2 , PARX , PARXM , q1 , q2 , QAPR , 
     &       wg , wl , XI , XIR , xp , xx , zt
      INTEGER*4 IAPR , Iexp , IPATH , ISEX , ist , istop , Jidim , k , 
     &          kk , LAMDA , LAMMAX , LDNUM , LEAD , MAGA , MULTI , n , 
     &          NDIM , ng , nl , NMAX
      INTEGER*4 NMAX1
      COMMON /MAP   / PARX(50,12,5) , PARXM(50,4,10,6) , XIR(6,50)
      COMMON /CXI   / XI(500)
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      Jidim = NMAX + 1
      IF ( MAGA(Iexp).NE.0 ) Jidim = 3*NMAX + 1
      ist = 1
      DO kk = 1 , 6
         IF ( MULTI(kk).NE.0 ) THEN
            istop = MULTI(kk) - 1 + ist
            DO k = ist , istop
               xx = ABS(XI(k))
               xx = xx/XIR(kk,Iexp)
               DO n = 1 , 7 , 3
                  IF ( MAGA(Iexp).NE.0 .OR. n.EQ.4 ) THEN
                     zt = QAPR(k,1,n)
                     zt = ABS(zt)
                     xp = 9.*xx
                     nl = INT(xp) + 1
                     wg = xp - DBLE(nl-1)
                     ng = nl + 1
                     wl = DBLE(nl) - xp
                     a = wg*PARXM(Iexp,1,ng,kk) + wl*PARXM(Iexp,1,nl,kk)
                     b = wg*PARXM(Iexp,2,ng,kk) + wl*PARXM(Iexp,2,nl,kk)
                     q1 = a*zt + b
                     a = wg*PARXM(Iexp,3,ng,kk) + wl*PARXM(Iexp,3,nl,kk)
                     b = wg*PARXM(Iexp,4,ng,kk) + wl*PARXM(Iexp,4,nl,kk)
                     q2 = a*zt + b
                     QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n)
                     QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n)
                     IF ( IAPR(k,1).EQ.IAPR(k,2) ) THEN
                        QAPR(k,1,n) = 0.
                        QAPR(k,2,n) = QAPR(k,2,n)/2.
                     ENDIF
                  ENDIF
               ENDDO
               IF ( MAGA(Iexp).NE.0 ) THEN
                  DO n = 2 , 6
                     IF ( n.NE.4 ) THEN
                        zt = QAPR(k,1,n)
                        zt = ABS(zt)
                        xp = 4.*xx
                        nl = INT(xp) + 1
                        wg = xp - DBLE(nl-1)
                        ng = nl + 1
                        wl = DBLE(nl) - xp
                        q1 = wg*PARX(Iexp,2*kk-1,ng)
     &                       + wl*PARX(Iexp,2*kk-1,nl)
                        q2 = wg*PARX(Iexp,2*kk,ng)
     &                       + wl*PARX(Iexp,2*kk,nl)
                        QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n)
                        QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            ist = istop + 1
         ENDIF
      ENDDO
      END
      SUBROUTINE POMNOZ(Acca,L,Iw,Ktoto,Img,Jidim)
      IMPLICIT NONE
      REAL*8 Acca , QAPR , sig , TCABS , test , u
      INTEGER*4 IAPR , IDIVE , Img , INHB , IPATH , ISEX , IVAR , Iw , 
     &          Jidim , k , kk , Ktoto , L , LERF , LMAXE , m , MAGA , 
     &          MAGEXC , mc , mc1
      INTEGER*4 MEMAX , MEMX6 , mw , mw1
      COMPLEX*16 ARM , ci
      COMMON /INHI  / INHB
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /AZ    / ARM(600,7)
      COMMON /APRX  / LERF , IDIVE(50,2)
      DATA ci/(0.,-1.)/
      sig = 1.
      IF ( L.NE.2 ) sig = -1.
      DO kk = 1 , Jidim
         ARM(kk,1) = ARM(kk,Iw)
      ENDDO
      DO k = 1 , 100
         Ktoto = Ktoto + 1
         DO m = 1 , MEMX6
            mw1 = IAPR(m,1)
            mc1 = IAPR(m,2)
            IF ( IPATH(mw1).NE.0 .AND. IPATH(mc1).NE.0 ) THEN
               mw = IPATH(mw1) + 1
               mc = IPATH(mc1) + 1
               IF ( Ktoto.GE.ISEX(mc1) ) THEN
                  IF ( Img.EQ.1 ) THEN
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,4)*ARM(mw,1)
                  ELSE
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,4)*ARM(mw,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPR(m,L,2)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,2)*ARM(mw-1,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPR(m,L,1)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPR(m,L,1)
     &                             *ARM(mw-1,1)
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,3)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPR(m,L,3)
     &                             *ARM(mw,1)
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,5)*ARM(mc+1,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,6)*ARM(mw+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPR(m,L,5)
     &                             *ARM(mw,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPR(m,L,6)*ARM(mc,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPR(m,L,7)*ARM(mc+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPR(m,L,7)
     &                             *ARM(mw+1,1)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         test = 0.
         DO m = 1 , Jidim
            ARM(m,1) = ARM(m,2)/k
            ARM(m,2) = (0.,0.)
            IF ( L.NE.1 ) ARM(m,1) = ARM(m,1)*ci
            ARM(m,Iw) = ARM(m,Iw) + ARM(m,1)
            IF ( k.GT.5 ) THEN
               u = TCABS(ARM(m,Iw))
               test = test + u*u
            ENDIF
         ENDDO
         IF ( ABS(test-1.).LT.Acca ) GOTO 99999
      ENDDO
      IF ( INHB.NE.1 ) LERF = 1
99999 END
      SUBROUTINE TENB(Icl,Bten,Lmax)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , Bten , CAT , ce , DIPOL , EN , fc , si , 
     &       SPIN , WTHREJ , x , ZPOL
      INTEGER*4 i , Icl , iha , ila , ilg , ind , isi , ISMAX , ISO , 
     &          ite , jm , jmp , k , kk , kp , l , ll , Lmax , lp , m
      INTEGER*4 mm , mp , ms , msp , NDIM , NMAX , NMAX1 , NSTART , 
     &          NSTOP
      COMPLEX*16 ARM
      DIMENSION Bten(1200)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /AZ    / ARM(600,7)
      iha = (-1)**INT(2.*SPIN(1)+.01)
      IF ( Icl.EQ.1 ) THEN
         ms = 16*(NMAX-1)
         DO i = 1 , ms
            Bten(i) = 0.
         ENDDO
      ENDIF
      DO i = 2 , NMAX
         ms = NSTART(i)
         IF ( ms.NE.0 ) THEN
            msp = NSTOP(i)
            si = SPIN(i)
            isi = INT(2.*si+.01)
            ce = SQRT(2.*si+1.)
            DO kp = 1 , 7 , 2
               k = kp - 1
               kk = 2*k
               IF ( isi.GE.k ) THEN
                  ila = -1
                  DO lp = 1 , kp
                     ila = -ila
                     l = lp - 1
                     ll = 2*l
                     ind = k*k/4 + lp + (i-2)*16
                     DO m = ms , msp
                        mm = m
                        mp = m + l
                        jm = INT(2.01*CAT(mm,3))
                        IF ( mp.GT.NSTOP(i) ) GOTO 4
                        ilg = (-1)**INT(si-CAT(mp,3))
                        jmp = -INT(2.01*CAT(mp,3))
                        fc = WTHREJ(isi,kk,isi,jmp,ll,jm)
                        ite = 1
 2                      IF ( ila.EQ.1 ) x = DBLE(ARM(mp,5))
     &                       *DBLE(ARM(mm,5)) + IMAG(ARM(mp,5))
     &                       *IMAG(ARM(mm,5))
                        IF ( ila.NE.1 ) x = DBLE(ARM(mp,5))
     &                       *IMAG(ARM(mm,5)) - DBLE(ARM(mm,5))
     &                       *IMAG(ARM(mp,5))
                        Bten(ind) = Bten(ind) + x*fc*ilg
                        IF ( ite.EQ.2 ) GOTO 6
 4                      IF ( iha.NE.1 .OR. Icl.NE.Lmax ) THEN
                           ite = 2
                           mp = mp - 2*l
                           IF ( mp.GE.NSTART(i) ) THEN
                              jmp = INT(2.01*CAT(mp,3))
                              jm = -jm
                              fc = WTHREJ(isi,kk,isi,jmp,ll,jm)
                              ilg = (-1)**INT(si+CAT(mp,3))
                              GOTO 2
                           ENDIF
                        ENDIF
 6                   ENDDO
                     IF ( Icl.EQ.Lmax ) Bten(ind) = Bten(ind)
     &                    *ce/(2.*SPIN(1)+1.)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
      END
      SUBROUTINE TENS(Bten)
      IMPLICIT NONE
      REAL*8 arg , Bten , DJMM , DSIGS , EPS , EROOT , FIEX , TETACM , 
     &       TREP , ZETA
      INTEGER*4 i , IAXS , IEXP , ind , inz , iph , ix , k , k1 , kp , 
     &          l , lp , lpp , lx , lxx , LZETA , NDIM , NMAX , NMAX1
      DIMENSION Bten(1200)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      ix = NMAX*28
      arg = 1.570796327 + TETACM(IEXP)/2.
      DO i = 1 , ix
         ZETA(i) = 0.
      ENDDO
      DO i = 2 , NMAX
         DO kp = 1 , 7 , 2
            k = kp - 1
            k1 = INT(DBLE(k)/2.+.01)
            IF ( k.EQ.0 ) THEN
               ind = (i-2)*16 + 1
               inz = (i-1)*28 + 1
               ZETA(inz) = Bten(ind)
            ELSE
               DO lp = 1 , kp
                  IF ( IAXS(IEXP).NE.0 .OR. lp.EQ.1 ) THEN
                     inz = (i-1)*28 + k1*7 + lp
                     l = lp - 1
                     DO lpp = 1 , kp
                        ind = k*k/4 + lpp + (i-2)*16
                        lx = lpp - 1
                        lxx = lx
 2                      iph = (-1)**(l+INT(DBLE(lxx)/2.))
                        ZETA(inz) = ZETA(inz) + Bten(ind)
     &                              *iph*DJMM(arg,k,lx,l)
                        IF ( lpp.NE.1 ) THEN
                           IF ( lx.GE.0 ) THEN
                              lx = -lx
                              lxx = lx - 1
                              GOTO 2
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      END
      FUNCTION DJMM(Beta,K,Kpp,Kp)
      IMPLICIT NONE
      REAL*8 B , b1 , b2 , be , BEQ , Beta , cb , ctb , djm , DJMM , f , 
     &       g , sb , sk , ul
      INTEGER*4 iczy , ifla , ifza , ill , j , ja , jb , jc , jd , K , 
     &          Kp , Kpp , lca , loc , mas , mis
      DIMENSION djm(525) , iczy(525)
      COMMON /IDENT / BEQ
      COMMON /CB    / B(20)
      SAVE djm
      ifza = 1
      IF ( Beta.LT.0. ) ifza = (-1)**(Kp+Kpp)
      sk = DBLE(K)
      ul = sk*((sk-1.)*(4.*sk+7)/6.+1.)
      lca = INT(ul+.1)
      loc = lca + (2*K+1)*Kp + Kpp + K + 1
      IF ( ABS(BEQ-ABS(Beta)).GT.1.E-6 ) THEN
         BEQ = ABS(Beta)
         DO ill = 1 , 525
            iczy(ill) = 0
         ENDDO
      ELSEIF ( iczy(loc).EQ.1 ) THEN
         DJMM = djm(loc)*ifza
         GOTO 99999
      ENDIF
      be = BEQ/2.
      cb = COS(be)
      sb = SIN(be)
      ifla = 0
      IF ( BEQ.GT..01 .AND. ABS(BEQ-6.2832).GT..01 ) ifla = 1
      IF ( ifla.NE.1 ) THEN
         IF ( Kp.EQ.Kpp ) THEN
            sb = 1.
         ELSE
            DJMM = 0.
            RETURN
         ENDIF
      ENDIF
      ctb = cb*cb/sb/sb
      ja = K + Kp + 1
      jb = K - Kp + 1
      jc = K + Kpp + 1
      jd = K - Kpp + 1
      b1 = B(ja)*B(jb)*B(jc)*B(jd)
      ja = Kp + Kpp
      jb = 2*K - Kp - Kpp
      IF ( ABS(BEQ-3.141592654).LT..01 .AND. ja.LT.0 ) ifla = 3
      IF ( ifla.EQ.3 ) cb = 1.
      f = (-1)**(K-Kp)*(cb**ja)*(sb**jb)*SQRT(b1)
      mis = 0
      IF ( ja.LT.0 ) mis = -ja
      mas = K - Kpp
      IF ( Kpp.LT.Kp ) mas = K - Kp
      ja = Kp + Kpp + mis + 1
      jb = K - Kpp - mis + 1
      jc = K - Kp - mis + 1
      jd = mis + 1
      b2 = B(ja)*B(jb)*B(jc)*B(jd)
      IF ( ifla.NE.3 ) THEN
         g = (-ctb)**mis/b2
         DJMM = g
         ja = mis + 1
         IF ( mas.GE.ja ) THEN
            DO j = ja , mas
               g = -g*ctb*(K-Kpp-j+1)*(K-Kp-j+1)/(Kp+Kpp+j)/j
               DJMM = DJMM + g
            ENDDO
         ENDIF
         IF ( ifla.EQ.0 ) DJMM = g
         DJMM = DJMM*f*ifza
         djm(loc) = DJMM/ifza
         iczy(loc) = 1
         RETURN
      ENDIF
      DJMM = f*ifza/((-sb*sb)**mis)/b2
      djm(loc) = DJMM/ifza
      iczy(loc) = 1
      RETURN
99999 END
      SUBROUTINE FTBM(Icll,Chisq,Idr,Ncall,Chilo,Bten)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AGELI , aval , Bten , CAT , CC , Chilo , 
     &       chis1 , CHIS11 , chish , Chisq , chisx , chx , CORF , 
     &       DIPOL , DYEX , EG , ELM , ELML
      REAL*8 ELMU , EMH , EN , EP , EPS , EROOT , fc , FIEX , fx , 
     &       polm , pr , prop , Q , SA , SPIN , TAU , TLBDG , UPL , 
     &       val , VINF
      REAL*8 wz , XA , XA1 , YEXP , YNRM , ZETA , ZPOL
      INTEGER*4 i1 , i11 , iapx , IAXS , Icll , idec , Idr , IDRN , 
     &          IEXP , iflg , IGRD , ii , ILE , ile1 , ile2 , ile3 , 
     &          ilin , indx , inko , INM
      INTEGER*4 inp , inpo , inpx , INTR , inzz , inzzz , IPATH , IPRM , 
     &          IPS1 , ISMAX , ISO , issp , ITAK2 , itemp , IVAR , ixx , 
     &          IY , IZ , IZ1 , izzz
      INTEGER*4 j , jj , jjgg , jjj , jk , jkl , jm , jmf , jmt , jmte , 
     &          jpp , jpz , JSKIP , jy , k , karm , kk , kk6 , kkx , kmt
      INTEGER*4 knm , KSEQ , kx , larm , lcc , lcou , LFL , LFL1 , 
     &          LFL2 , licz , lix , llx , lm , LMAX , LMAXE , lmh , 
     &          LNY , loc , loch , loct
      INTEGER*4 lp , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , lpit , lput , lpx , 
     &          lpxd , ls , lst
      INTEGER*4 luu , lx , LZETA , MAGA , MAGEXC , MEMAX , MEMX6 , 
     &          NANG , Ncall , NDIM , NEXPT , NICC , NLIFT , nlin , 
     &          NMAX , NMAX1 , nowr , npoz , nrest , NSTART
      INTEGER*4 NSTOP , NWR , nwyr , NYLDE
      COMPLEX*16 ARM
      DIMENSION jmte(6) , prop(6) , Bten(1200)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      COMMON /ILEWY / NWR
      COMMON /CH1T  / CHIS11
      COMMON /IGRAD / IGRD
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      COMMON /UWAGA / ITAK2
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CLM   / LMAX
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /AZ    / ARM(600,7)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /PRT   / IPRM(20)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /SKP   / JSKIP(50)
      COMMON /LIFE  / NLIFT
      COMMON /LOGY  / LNY , INTR , IPS1
      issp = 0
      Chilo = 0.
      fx = 2.*SPIN(1) + 1.
      Chisq = 0.
      LFL = 0
      chis1 = 0.
      ixx = NDIM*MEMAX + LP11
      DO i1 = 1 , ixx
         ZETA(i1) = 0.
      ENDDO
      DO ii = 1 , LP6
         ILE(ii) = 1
      ENDDO
      itemp = 0
      NWR = 0
      iapx = 1
      DO jkl = 1 , NEXPT
         IEXP = jkl
         IGRD = 0
         LFL2 = 1
         IF ( ITAK2.EQ.-1 ) THEN
            DO larm = 1 , 4
               DO karm = 1 , LP10
                  ARM(karm,larm) = (0.,0.)
               ENDDO
            ENDDO
         ENDIF
         iflg = 0
         IF ( IEXP.NE.1 ) THEN
            kk = NANG(IEXP)
            DO jjj = 1 , LP6
               ILE(jjj) = ILE(jjj) + NYLDE(IEXP-1,jjj)
            ENDDO
         ENDIF
         lp = 3
         IF ( JSKIP(jkl).EQ.0 ) GOTO 200
         IF ( MAGA(IEXP).EQ.0 ) lp = 1
         IF ( Ncall.EQ.0 ) GOTO 150
         IF ( Icll.EQ.4 ) GOTO 100
 50      loch = LP3*(MEMAX-1) + NMAX + LP11
         DO k = 1 , loch
            ZETA(k) = 0.
         ENDDO
         CALL LOAD(IEXP,1,2,0.D0,jj)
         DO k = 1 , LMAX
            fc = 2.
            IF ( k.EQ.LMAX ) fc = 1.
            IF ( DBLE(INT(SPIN(1))).LT.SPIN(1) ) fc = 2.
            loc = 0
            polm = DBLE(k-1) - SPIN(1)
            CALL LOAD(IEXP,3,2,polm,jj)
            CALL PATH(jj)
            CALL LOAD(IEXP,2,2,polm,jj)
            CALL APRAM(IEXP,1,1,jj,ACCA)
            IF ( Ncall.NE.0 ) THEN
               IF ( Icll.NE.3 ) THEN
                  DO indx = 1 , MEMX6
                     CALL APRAM(IEXP,0,indx,jj,ACCA)
                     kx = 0
                     DO i11 = 1 , NMAX
                        IF ( NSTART(i11).NE.0 ) THEN
                           loc = LP3*(indx-1) + i11 + LP11
                           jpp = INT(2.*SPIN(i11)+1.)
                           lpx = MIN(lp,jpp)
                           IF ( ISO.NE.0 ) lpx = NSTOP(i11)
     &                          - NSTART(i11) + 1
                           DO lpxd = 1 , lpx
                              kx = kx + 1
                              ZETA(loc) = ZETA(loc) + fc*DBLE(ARM(kx,5))
     &                           *DBLE(ARM(kx,6))
     &                           /fx + fc*IMAG(ARM(kx,5))
     &                           *IMAG(ARM(kx,6))/fx
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
            CALL TENB(k,Bten,LMAX)
         ENDDO
         IF ( loc.NE.0 ) THEN
            REWIND 14
            WRITE (14,*) (ZETA(i11),i11=LP8,loch)
         ENDIF
         CALL TENS(Bten)
         IF ( Ncall.EQ.0 ) GOTO 200
         IF ( Icll.GE.2 ) GOTO 200
         llx = 28*NMAX
         DO lx = 1 , llx
            ZETA(LP9+lx) = ZETA(lx)
         ENDDO
         IF ( Icll.NE.1 ) GOTO 200
 100     iapx = 0
         issp = 1
         CALL LOAD(IEXP,1,1,0.D0,jj)
         CALL ALLOC(ACCUR)
         CALL SNAKE(IEXP,ZPOL)
         CALL SETIN
         DO k = 1 , LMAX
            polm = DBLE(k-1) - SPIN(1)
            CALL LOAD(IEXP,2,1,polm,kk)
            IF ( IPRM(7).EQ.-1 ) WRITE (22,99001) polm , IEXP
99001       FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',1F4.1,5X,
     &              'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,'SPIN',2X,'M',5X,
     &              'REAL AMPLITUDE',2X,'IMAGINARY AMPLITUDE'//)
            CALL STING(kk)
            CALL PATH(kk)
            CALL INTG(IEXP)
            CALL TENB(k,Bten,LMAX)
            IF ( IPRM(7).EQ.-1 ) THEN
               DO j = 1 , ISMAX
                  WRITE (22,99002) INT(CAT(j,1)) , CAT(j,2) , CAT(j,3) , 
     &                             DBLE(ARM(j,5)) , IMAG(ARM(j,5))
99002             FORMAT (7X,1I2,3X,1F4.1,2X,1F4.1,2X,1E14.6,2X,1E14.6)
               ENDDO
            ENDIF
         ENDDO
         CALL TENS(Bten)
         IF ( IPRM(7).EQ.-1 ) THEN
            DO jjgg = 2 , NMAX
               loct = (jjgg-1)*28 + 1
               WRITE (22,99003) jjgg , ZETA(loct)
99003          FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
            ENDDO
         ENDIF
         GOTO 200
 150     IF ( iflg.EQ.1 ) THEN
            itemp = 1
            iflg = 2
            GOTO 50
         ELSE
            IF ( iflg.EQ.2 ) GOTO 300
            itemp = 2
            iflg = 1
            GOTO 100
         ENDIF
 200     CALL CEGRY(Chisq,itemp,Chilo,Idr,nwyr,Icll,issp,0)
         issp = 0
         IF ( Ncall.EQ.0 .AND. JSKIP(jkl).NE.0 ) THEN
            IF ( Ncall.EQ.0 ) GOTO 150
            GOTO 200
         ELSE
            NWR = NWR + nwyr
            IF ( Icll.LE.2 .AND. JSKIP(jkl).NE.0 ) THEN
               IF ( IEXP.EQ.1 ) chish = CHIS11
               IF ( Icll.EQ.1 ) chis1 = CHIS11
               IF ( Icll.EQ.0 ) chis1 = Chisq
               LFL2 = 0
               IGRD = 1
               IF ( ITAK2.EQ.-1 ) LFL = 1
               REWIND 14
               READ (14,*) (ZETA(i11),i11=LP8,loch)
               DO larm = 1 , 4
                  DO karm = 1 , LP10
                     ARM(karm,larm) = (0.,0.)
                  ENDDO
               ENDDO
               chisx = 0.
               llx = 28*NMAX
               DO lix = 1 , llx
                  ZETA(LP9+lix) = ZETA(lix)
               ENDDO
               CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,1)
               DO knm = 1 , MEMAX
                  INM = knm
                  chisx = 0.
                  EMH = ELM(INM)
                  ELM(INM) = 1.05*EMH
                  lcc = LP3*(INM-1) + LP11
                  DO lst = 2 , NMAX
                     wz = ZETA(lst+lcc)
                     inpx = (lst-1)*28
                     DO jy = 1 , 4
                        inp = inpx + (jy-1)*7
                        IF ( jy.EQ.1 ) pr = ZETA(LP13+inp) + 1.E-12
                        jmf = 2*jy - 1
                        IF ( IAXS(IEXP).EQ.0 ) jmf = 1
                        DO jm = 1 , jmf
                           inp = inp + 1
                           ZETA(inp) = ZETA(inp+LP9)*(1.+.1*EMH*wz/pr)
                        ENDDO
                     ENDDO
                  ENDDO
                  CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,0)
                  ELM(INM) = EMH
               ENDDO
               IF ( ITAK2.EQ.-1 .AND. LFL1.NE.0 ) THEN
                  IF ( IPRM(17).NE.0 ) THEN
                     kmt = ABS(IPRM(17))
                     WRITE (22,99004) IEXP
99004                FORMAT (1X///20X,'EXPERIMENT',11X,1I2,5X,
     &                       'D(LOG(P))/D(LOG(ME)) MAP'/20X,52('-')///)
                     nlin = (NMAX-2)/6 + 1
                     nrest = NMAX - 1 - 6*(nlin-1)
                     DO ilin = 1 , nlin
                        npoz = 6
                        IF ( ilin.EQ.nlin ) npoz = nrest
                        inpo = (ilin-1)*6 + 2
                        inko = inpo + npoz - 1
                        lpit = 0
                        DO lm = inpo , inko
                           lpit = lpit + 1
                           jmte(lpit) = lm
                        ENDDO
                        WRITE (22,99005) (jmte(lm),lm=1,lpit)
99005                   FORMAT (5X,'LEVEL',6(8X,1I2,9X))
                        WRITE (22,99006)
     &                         (ZETA(LP13+(jpz-1)*28),jpz=inpo,inko)
99006                   FORMAT (1X,'EXC.PROB.',6(5X,1E10.4,4X))
                        DO jmt = 1 , kmt
                           lput = 0
                           DO ls = inpo , inko
                              lput = lput + 1
                              prop(lput) = 0.
                              DO lm = 1 , MEMX6
                                 inzz = ls + LP3*(lm-1) + LP11
                                 inzzz = LP13 + (ls-1)*28
                                 IF ( ABS(ZETA(inzzz)).LT.1.E-20 )
     &                                ZETA(inzzz) = 1.E-20
                                 val = 2.*ELM(lm)*ZETA(inzz)/ZETA(inzzz)
                                 aval = ABS(val)
                                 IF ( aval.GT.ABS(prop(lput)) ) THEN
                                    prop(lput) = val
                                    lmh = lm
                                    jmte(lput) = lm
                                 ENDIF
                              ENDDO
                              izzz = (lmh-1)*LP3 + LP11 + ls
                              ZETA(izzz) = 0.
                           ENDDO
                           WRITE (22,99007)
     &                            (jmte(lcou),prop(lcou),lcou=1,npoz)
99007                      FORMAT (10X,6(2X,'(',1X,1I3,1X,1E8.2,')',2X))
                        ENDDO
                     ENDDO
                     REWIND 14
                     READ (14,*) (ZETA(i11),i11=LP8,loch)
                     IF ( IPRM(17).LT.0 ) GOTO 300
                  ENDIF
                  LFL = 0
                  WRITE (22,99008) IEXP
99008             FORMAT (10X,'EXPERIMENT',1X,1I2/10X,
     &                    'D(LOG(Y)/D(LOG(ME))',//)
                  ile1 = ILE(1) + NYLDE(IEXP,1) - 1
                  ile3 = ILE(1)
                  licz = 0
                  DO ile2 = ile3 , ile1
                     licz = licz + 1
                     idec = IY(ile2,1)
                     IF ( idec.GT.1000 ) idec = idec/1000
                     luu = 6*licz - 5
                     jk = (luu-1)/LP10 + 1
                     kk = luu - LP10*(jk-1)
                     kk6 = kk + 5
                     WRITE (22,99009) KSEQ(idec,3) , KSEQ(idec,4) , 
     &                                (INT(DBLE(ARM(kkx,jk))),
     &                                IMAG(ARM(kkx,jk)),kkx=kk,kk6)
99009                FORMAT (2X,1I2,'--',1I2,5X,
     &                       6('(',1I3,2X,1E8.2,')',3X))
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
 300  ENDDO
      IF ( ITAK2.EQ.-1 .AND. Icll.LT.2 ) ITAK2 = 0
      IF ( Ncall.NE.0 ) THEN
         IF ( Icll.LE.2 ) THEN
            IF ( Icll.EQ.1 ) CALL CEGRY(Chisq,itemp,Chilo,Idr,nowr,7,
     &                                  issp,0)
         ENDIF
         CALL BRANR(Chisq,NWR,Chilo)
         CALL MIXR(NWR,0,Chisq,Chilo)
         CALL CHMEM(NWR,Chisq,Chilo)
         NWR = NWR + NLIFT
         Chisq = Chisq/NWR
         IF ( INTR.NE.0 ) THEN
            chx = Chisq
            Chisq = Chilo
            Chilo = chx
         ENDIF
      ENDIF
      RETURN
      END
      SUBROUTINE MINI(Chisq,Chiok,Nptl,Conv,Imode,Idr,Xtest,Ips,Is,Jjh,
     &                Bten)
      IMPLICIT NONE
      REAL*8 a , a0 , a1 , b , Bten , c , ccd , chd , chil , chilo , 
     &       Chiok , chirf , CHIS11 , chis12 , chis13 , chisf , chisp , 
     &       Chisq , chiss , chl
      REAL*8 chx , cmax , Conv , CORF , crit , DEVD , DEVU , dl , 
     &       DLOCK , dm , DYEX , ELM , ELMH , ELML , ELMU , EMH , f1 , 
     &       f2 , flt , GRAD
      REAL*8 gradp , HLMLM , ht , p , q , rfk , SA , sel , shl , sumg1 , 
     &       sumg2 , sumht , UPL , uxa , xkat , Xtest , YEXP , YNRM
      INTEGER*4 i , icl1 , icl2 , icount , ICS , Idr , IDRN , IFBFL , 
     &          iht , iin , ILE , Imode , indx1 , INM , inmx , ino , 
     &          INTR , ipas , ipm , IPRM
      INTEGER*4 Ips , IPS1 , Is , istec , ITAK2 , itf , IVAR , IY , j , 
     &          jcoup , jcp , JENTR , jin , Jjh , jjj , jlin , jnm , 
     &          jpr , jsa , jst
      INTEGER*4 KFERR , kh2 , kkk , KVAR , l , LFL , LFL1 , LFL2 , 
     &          LMAXE , lnm , LNY , LOCKF , LOCKS , LP1 , LP10 , LP11 , 
     &          LP12 , LP13 , LP14 , LP2
      INTEGER*4 LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , MAGEXC , MEMAX , 
     &          MEMX6 , metf , mvfl , ncall , nlinn , NLOCK , noflg , 
     &          Nptl , NWR , NYLDE
      DIMENSION ipm(10) , Bten(1200) , gradp(500)
      COMMON /DUMM  / GRAD(500) , HLMLM(500) , ELMH(500)
      COMMON /ILEWY / NWR
      COMMON /CH1T  / CHIS11
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /UWAGA / ITAK2
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /DFTB  / DEVD(500) , DEVU(500)
      COMMON /PRT   / IPRM(20)
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /SEL   / KVAR(500)
      COMMON /FIT   / LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
      COMMON /ERRAN / KFERR
      COMMON /LOGY  / LNY , INTR , IPS1
      COMMON /ERCAL / JENTR , ICS
      DO i = 1 , MEMAX
         gradp(i) = 0.
      ENDDO
      icount = 0
      lnm = 0
      LNY = 0
      INTR = 0
      metf = 0
      LFL1 = 0
      ncall = 0
      ITAK2 = 0
      IF ( Imode.LT.2000 ) THEN
         icl1 = 0
         icl2 = 3
         IF ( Imode.GE.1100 ) metf = 1
         IF ( (Imode-1000-100*metf).GE.10 ) lnm = 1
         IF ( (Imode-1000-100*metf-10*lnm).EQ.1 ) LNY = 1
         IF ( JENTR.EQ.1 ) GOTO 200
         IF ( ICS.NE.0 ) THEN
            REWIND 11
            DO jnm = 1 , LP4
               READ (11) (CORF(jnm,kh2),kh2=1,LP6)
            ENDDO
            ICS = 0
            GOTO 200
         ENDIF
      ELSE
         icl1 = 1
         IF ( Imode.GE.2100 ) metf = 1
         IF ( (Imode-2000-100*metf).GE.10 ) lnm = 1
         IF ( (Imode-2000-100*metf-10*lnm).EQ.1 ) LNY = 1
         icl2 = 4
         IF ( Ips.NE.0 ) THEN
            IF ( Ips.EQ.1 ) THEN
               IF ( IPRM(4).EQ.-1 ) ITAK2 = -2
            ELSE
               IF ( IPRM(4).LT.0 ) ITAK2 = -2
            ENDIF
            icl1 = 4
            IF ( ITAK2.EQ.-2 ) icl1 = 1
            IF ( icl1.EQ.4 ) GOTO 200
         ENDIF
      ENDIF
 100  CALL FTBM(0,chiss,Idr,0,chl,Bten)
      REWIND 11
      DO jnm = 1 , LP4
         WRITE (11) (CORF(jnm,kh2),kh2=1,LP6)
      ENDDO
      IF ( IPS1.EQ.0 ) RETURN
 200  noflg = 0
      ncall = 1
 300  sumht = 0.
      IF ( LNY.EQ.1 ) INTR = 1
      LFL1 = 1
      ITAK2 = ITAK2 + 1
      icount = icount + 1
      IF ( icount.GT.Nptl ) THEN
         IF ( KFERR.EQ.1 ) RETURN
         IF ( Ips.EQ.0 ) WRITE (22,99001) Nptl
99001    FORMAT (5X,'MINIMIZATION STOPPED-NUMBER OF STEPS NPTL=',1I5,1X,
     &           'EXCEEDED')
         IF ( Ips.EQ.0 ) WRITE (22,99010) chil
         INTR = 0
         RETURN
      ELSE
         IF ( ITAK2.EQ.IPRM(4) ) ITAK2 = -1
         IF ( ITAK2.EQ.-1 ) THEN
            IF ( KFERR.NE.1 ) THEN
               CALL FTBM(3,chd,Idr,1,chl,Bten)
               CHIS11 = chd*NWR
               CALL FTBM(icl1,Chisq,Idr,ncall,chilo,Bten)
            ENDIF
         ENDIF
         IF ( Ips.EQ.1 ) RETURN
         IF ( icl1.EQ.1 ) CALL FTBM(4,Chisq,Idr,ncall,chilo,Bten)
         IF ( IPRM(8).EQ.-1 .OR. IPRM(13).EQ.-1 ) THEN
            IF ( IPRM(8).EQ.-1 ) IPRM(8) = -2
            IF ( IPRM(13).EQ.-1 ) IPRM(13) = -2
            CALL FTBM(4,ccd,Idr,ncall,chl,Bten)
            IF ( Ips.EQ.2 ) RETURN
         ENDIF
         CALL FTBM(3,chis12,Idr,ncall,chilo,Bten)
         IF ( icl1.EQ.0 ) Chisq = chis12
         uxa = Chisq
         IF ( INTR.EQ.1 ) uxa = chilo
         ipas = 0
         IF ( uxa.LT.Chiok ) Chisq = uxa
         IF ( uxa.LT.Chiok ) GOTO 700
      ENDIF
 400  ino = 1
      IF ( metf.EQ.1 ) ipas = ipas + 1
      IF ( IFBFL.EQ.1 ) ino = 2
      DO jjj = 1 , ino
         DO jnm = 1 , MEMAX
            GRAD(jnm) = 0.
            IF ( IVAR(jnm).EQ.1 .OR. IVAR(jnm).EQ.2 ) THEN
               DO jcoup = 1 , MEMAX
                  ELMH(jcoup) = ELM(jcoup)
               ENDDO
               DO jcoup = 1 , MEMAX
                  IF ( jnm.NE.jcoup ) THEN
                     IF ( IVAR(jcoup).LT.1000 ) GOTO 410
                     jcp = IVAR(jcoup) - 1000
                     IF ( jcp.NE.jnm ) GOTO 410
                     IF ( IVAR(jnm).EQ.0 ) GOTO 410
                  ENDIF
                  flt = 1.01
                  IF ( jjj.EQ.2 ) flt = .99
                  ELM(jcoup) = ELMH(jcoup)*flt
 410           ENDDO
               CALL FTBM(3,chis13,Idr,ncall,chx,Bten)
               IF ( jjj.EQ.1 ) HLMLM(jnm) = chis13
               IF ( IFBFL.NE.1 .OR. jjj.NE.1 ) THEN
                  IF ( jjj.EQ.2 ) chis12 = chis13
                  GRAD(jnm) = 100.*(HLMLM(jnm)-chis12)/ELMH(jnm)
                  IF ( IFBFL.EQ.1 ) GRAD(jnm) = GRAD(jnm)/2.
                  IF ( lnm.EQ.1 ) GRAD(jnm) = GRAD(jnm)*ABS(ELMH(jnm))
               ENDIF
               DO jcoup = 1 , MEMAX
                  ELM(jcoup) = ELMH(jcoup)
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      IF ( KFERR.EQ.1 ) THEN
         GRAD(Jjh) = 0.
         IF ( Is.EQ.1 .AND. icount.EQ.1 ) WRITE (33,*)
     &        (NWR*GRAD(jnm),jnm=1,MEMAX)
      ENDIF
      IF ( metf.EQ.1 .AND. ipas.EQ.2 ) THEN
         DO jnm = 1 , MEMAX
            ELM(jnm) = DEVU(jnm)
         ENDDO
         shl = dm/20./sumg2
         sumg1 = 0.
         DO jnm = 1 , MEMAX
            GRAD(jnm) = (DEVD(jnm)*sumg2-GRAD(jnm))/shl
            sumg1 = sumg1 + GRAD(jnm)*GRAD(jnm)
         ENDDO
         sumg1 = SQRT(sumg1)
         p = 0.
         DO jnm = 1 , MEMAX
            GRAD(jnm) = GRAD(jnm)/sumg1
            DEVU(jnm) = ELM(jnm)
            sel = dm*GRAD(jnm)/100.
            IF ( lnm.EQ.1 ) sel = sel*ABS(DEVU(jnm))
            p = p + DEVD(jnm)*GRAD(jnm)
            ELM(jnm) = ELM(jnm) + sel
         ENDDO
         CALL FTBM(3,chis13,Idr,ncall,chx,Bten)
         shl = dm/100.
         DO jnm = 1 , MEMAX
            sel = dm*GRAD(jnm)/50.
            IF ( lnm.EQ.1 ) sel = sel*ABS(DEVU(jnm))
            ELM(jnm) = ELM(jnm) - sel
         ENDDO
         CALL FTBM(3,chis12,Idr,ncall,chx,Bten)
         q = (chis12+chis13-2.*Chisq)/shl/shl
         a0 = q*sumg2/sumg1 - p
         a1 = p*p - 1.
         sumg1 = SQRT(a0*a0+a1*a1+2.*a0*a1*p)
         DO jnm = 1 , MEMAX
            ELM(jnm) = DEVU(jnm)
            GRAD(jnm) = (GRAD(jnm)*a1+DEVD(jnm)*a0)/sumg1
         ENDDO
      ELSE
         sumg2 = 0.
         DO jnm = 1 , MEMAX
            IF ( IVAR(jnm).EQ.1 .OR. IVAR(jnm).EQ.2 ) sumg2 = sumg2 + 
     &           GRAD(jnm)*GRAD(jnm)
         ENDDO
         IF ( sumg2.LT.1.E-10 ) GOTO 800
         sumg2 = SQRT(sumg2)
         DO jnm = 1 , MEMAX
            GRAD(jnm) = GRAD(jnm)/sumg2
         ENDDO
         IF ( metf.NE.0 ) THEN
            dm = 0.
            DO jnm = 1 , MEMAX
               IF ( IVAR(jnm).EQ.2 .OR. IVAR(jnm).EQ.1 ) dm = dm + 
     &              ELM(jnm)*ELM(jnm)*GRAD(jnm)*GRAD(jnm)
            ENDDO
            dm = SQRT(dm)
            DO jnm = 1 , MEMAX
               DEVD(jnm) = GRAD(jnm)
               DEVU(jnm) = ELM(jnm)
               sel = dm*GRAD(jnm)/20.
               IF ( lnm.EQ.1 ) sel = sel*ABS(ELM(jnm))
               ELM(jnm) = ELM(jnm) - sel
            ENDDO
            IF ( IFBFL.EQ.0 ) CALL FTBM(3,chis12,Idr,ncall,chx,Bten)
            GOTO 400
         ENDIF
      ENDIF
      LFL1 = 0
      IF ( lnm.NE.0 ) THEN
         DO jnm = 1 , MEMAX
            GRAD(jnm) = GRAD(jnm)*ABS(ELM(jnm))
         ENDDO
      ENDIF
      sumg1 = 0.
      DO jnm = 1 , MEMAX
         sumg1 = sumg1 + GRAD(jnm)*GRAD(jnm)
      ENDDO
      sumg1 = SQRT(sumg1)
      DO jnm = 1 , MEMAX
         GRAD(jnm) = GRAD(jnm)/sumg1
      ENDDO
      IF ( LNY.EQ.1 ) Chisq = chilo
      IF ( noflg.EQ.0 ) chirf = Chisq
      noflg = 1
      chil = Chisq
      IF ( KFERR.NE.1 ) THEN
         IF ( MOD(icount,IPRM(5)).EQ.0 .OR. icount.EQ.1 )
     &        WRITE (22,99010) Chisq
c      WRITE(*,54)CHISQ
         IF ( MOD(icount,IPRM(6)).EQ.0 ) THEN
            WRITE (22,99002)
99002       FORMAT (20X,'GRADIENT'//)
            nlinn = MEMAX/10 + 1
            DO jlin = 1 , nlinn
               jsa = (jlin-1)*10 + 1
               DO jin = 1 , 10
                  ipm(jin) = jsa + jin - 1
               ENDDO
               jst = MIN(jsa+9,MEMAX)
               jpr = MIN(10,MEMAX-jsa+1)
               WRITE (22,99003) (ipm(jin),jin=1,jpr)
99003          FORMAT (5X,10(5X,1I3,4X))
               WRITE (22,99004) (GRAD(jin),jin=jsa,jst)
99004          FORMAT (5X,10(1X,1E10.4,1X)/)
            ENDDO
         ENDIF
      ENDIF
      IF ( chil.LT.Chiok ) GOTO 700
      DO l = 1 , MEMAX
         HLMLM(l) = ELM(l)
      ENDDO
      DO l = 1 , MEMAX
         IF ( ABS(GRAD(l)).LE.DLOCK .AND. LOCKS.EQ.1 .AND. 
     &        icount.EQ.1 .AND. IVAR(l).LE.999 .AND. IVAR(l).NE.0 ) THEN
            IF ( KFERR.NE.1 ) KVAR(l) = 0
            IF ( KFERR.NE.1 ) WRITE (22,99005) l , GRAD(l)
99005       FORMAT (1X,'MATRIX ELEMENT',1X,1I3,1X,'LOCKED',3X,
     &              'DERIVATIVE=',1E14.6)
            IVAR(l) = 0
         ENDIF
      ENDDO
      istec = 0
 500  DO j = 1 , MEMAX
         ELMH(j) = ELM(j)
      ENDDO
      istec = istec + 1
      cmax = 0.
      INTR = 0
      inmx = 1
      DO iht = 1 , MEMAX
         IF ( ABS(GRAD(iht)).GT.cmax ) THEN
            cmax = ABS(GRAD(iht))
            inmx = iht
         ENDIF
      ENDDO
      ht = .01*ABS(ELM(inmx))/cmax
      mvfl = 0
      IF ( icount.NE.1 .AND. istec.EQ.1 ) THEN
         xkat = 0.
         DO j = 1 , MEMAX
            xkat = xkat + GRAD(j)*gradp(j)
         ENDDO
         DO j = 1 , MEMAX
            gradp(j) = GRAD(j)
         ENDDO
         IF ( xkat.GE..8 ) THEN
            a = 0.
            DO j = 1 , MEMAX
               IF ( IVAR(j).NE.0 .AND. IVAR(j).LE.999 ) THEN
                  a = MAX(a,ABS(GRAD(j)))
                  IF ( ABS(a-ABS(GRAD(j))).LT.1.E-9 ) iin = j
               ENDIF
            ENDDO
            WRITE (22,99011) iin
            IVAR(iin) = 0
            GRAD(iin) = 0.
            gradp(iin) = 0.
         ENDIF
      ENDIF
 600  DO j = 1 , MEMAX
         ELM(j) = ELMH(j) - ht*GRAD(j)
      ENDDO
      DO j = 1 , MEMAX
         IF ( IVAR(j).GE.1000 ) THEN
            indx1 = IVAR(j) - 1000
            ELM(j) = ELM(indx1)*SA(j)
         ENDIF
      ENDDO
      IF ( mvfl.EQ.0 ) THEN
         CALL FTBM(icl2,chisp,Idr,ncall,chilo,Bten)
         DO j = 1 , MEMAX
            ELM(j) = 2.*ELMH(j) - ELM(j)
         ENDDO
         CALL FTBM(icl2,chisf,Idr,ncall,chilo,Bten)
         c = (chisp+chisf-2.*chil)/ht/ht
         b = (chisp-chisf)/ht/2.
         dl = b*b - 2.*c*chil
         IF ( dl.GT.0. ) THEN
            f1 = chil
            f2 = b
         ELSE
            f1 = b
            f2 = c
         ENDIF
         mvfl = 1
         IF ( ABS(f2).LT.1.E-10 ) THEN
            ht = 1.
         ELSE
            ht = -f1/f2
         ENDIF
         GOTO 600
      ELSE
         CALL LIMITS
         CALL FTBM(icl2,Chisq,Idr,ncall,chilo,Bten)
         IF ( Chisq.GE.chil ) THEN
            ht = ht/2.
            IF ( ABS(ht).GE.Conv ) GOTO 600
         ELSE
            chil = Chisq
            sumht = sumht + ht
            IF ( ABS(ht/sumht).GE..01 ) GOTO 500
         ENDIF
         crit = 0.
         DO jjj = 1 , MEMAX
            crit = crit + (ELM(jjj)-HLMLM(jjj))**2
         ENDDO
         crit = SQRT(crit)
         IF ( crit.LT.Conv ) GOTO 800
         IF ( Chisq.GE.Chiok ) THEN
            rfk = chirf/Chisq
            IF ( rfk.LE.Xtest .OR. icount.GE.Nptl ) GOTO 300
            GOTO 100
         ENDIF
      ENDIF
 700  chil = Chisq
      IF ( Ips.EQ.0 ) WRITE (22,99006) icount
99006 FORMAT (5X,'AT STEP',1X,1I5,1X,'CHISQ CRITERION FULFILLED')
      IF ( Ips.EQ.0 ) WRITE (22,99010) chil
      RETURN
 800  IF ( LOCKF.EQ.0 ) THEN
         IF ( Chisq.GE.chil ) THEN
            DO jjj = 1 , MEMAX
               ELM(jjj) = ELMH(jjj)
            ENDDO
         ENDIF
         IF ( KFERR.EQ.1 ) RETURN
         IF ( Ips.EQ.0 ) WRITE (22,99007) icount , crit
99007    FORMAT (5X,'AT STEP',1X,1I5,'CONVERGENCE ACHIEVED(',1E14.6,')')
         IF ( Ips.EQ.0 ) WRITE (22,99010) MIN(chil,Chisq)
         INTR = 0
         RETURN
      ELSE
         DO kkk = 1 , NLOCK
            a = 0.
            iin = 1
            DO jjj = 1 , MEMAX
               IF ( IVAR(jjj).NE.0 .AND. IVAR(jjj).LE.999 ) THEN
                  a = MAX(a,ABS(GRAD(jjj)))
                  IF ( ABS(a-ABS(GRAD(jjj))).LT.1.E-9 ) iin = jjj
               ENDIF
            ENDDO
            IVAR(iin) = 0
            WRITE (22,99011) iin
         ENDDO
         itf = 0
         DO jjj = 1 , MEMAX
            IF ( IVAR(jjj).LE.999 ) THEN
               IF ( IVAR(jjj).NE.0 ) itf = itf + 1
            ENDIF
         ENDDO
         IF ( itf.EQ.1 ) THEN
            metf = 0
            WRITE (22,99008)
99008       FORMAT (2x,'Warning - only one matrix element free',//2x,
     &              'Mode reset to single gradient, execution continues'
     &              ,/)
         ENDIF
         IF ( itf.NE.0 ) GOTO 300
         WRITE (22,99009)
99009    FORMAT (1X/////5X,'*****',2X,'ALL MATRIX ELEMENTS LOCKED!',2X,
     &           '*****'/////)
         INTR = 0
         RETURN
      ENDIF
99010 FORMAT (5X,'*** CHISQ=',1E14.6,1X,'***')
99011 FORMAT (1X/5X,'MATRIX ELEMENT',1X,1I3,1X,'LOCKED!')
      END
      SUBROUTINE CEGRY(Chisq,Itemp,Chilo,Idr,Nwyr,Icall,Issp,Iredv)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AGELI , AKS , BETAR , CC , ccc , ccd , 
     &       Chilo , Chisq , CNOR , CNOR1 , CNOR2 , cnr , cocos , CORF , 
     &       d , decen , DELTA , DEV
      REAL*8 DIPOL , DIX , dl , DQ , DSIGS , DYEX , effi , EG , EMH , 
     &       EN , ENDEC , ENZ , EP , EPS , EROOT , fi0 , fi1 , fic , 
     &       FIEX , figl
      REAL*8 fm , g , gth , ODL , PART , PARTL , Q , QCEN , rik , rl , 
     &       rx , ry , rys , rz , sf , sgm , SGW , SPIN , SUBCH1 , 
     &       SUBCH2
      REAL*8 sum3 , SUMCL , sumpr , TACOS , TAU , TETACM , tetrc , 
     &       tfac , thc , TLBDG , TREP , UPL , VACDP , VINF , wf , XA , 
     &       XA1 , XNOR , YEXP , YGN
      REAL*8 YGP , YNRM , ZPOL
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
     &          MCFIX , na , NANG
      INTEGER*4 NDIM , NDST , NEXPT , nf , nf1 , ni , ni1 , NICC , 
     &          NLIFT , NMAX , NMAX1 , Nwyr , NYLDE
      CHARACTER*4 wupl , war , op2
      DIMENSION lic(32) , lth(500) , cnr(32,50)
      COMMON /CLUST / ICLUST(50,200) , LASTCL(50,20) , SUMCL(20,500) , 
     &                IRAWEX(50)
      COMMON /ODCH  / DEV(500)
      COMMON /TUTAJ / PART(32,50,2) , PARTL(32,50,2)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /BREC  / BETAR(50)
      COMMON /DIMX  / DIX(4) , ODL(200)
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      COMMON /CINIT / CNOR(32,75) , INNR , MCFIX , CNOR1(32,75) , 
     &                CNOR2(32,75)
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
                     PARTL(jj,IEXP,jk) = 0.
                     PART(jj,IEXP,jk) = 0.
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
                     ixm = KSEQ(ixl,3)
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
                  ni = KSEQ(l,3)
                  nf = KSEQ(l,4)
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
                           ni1 = KSEQ(l1,3)
                           nf1 = KSEQ(l1,4)
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
                           PART(k9,IEXP,1) = PART(k9,IEXP,1) + YGN(l)
     &                        *YGN(l)*wf*wf/dl
                           PART(k9,IEXP,2) = PART(k9,IEXP,2) - 2.*YGN(l)
     &                        *wf*YEXP(k9,lu)/dl
                           sumpr = sumpr + YEXP(k9,lu)*YEXP(k9,lu)/dl
                           PARTL(k9,IEXP,1) = PARTL(k9,IEXP,1)
     &                        + YEXP(k9,lu)*YEXP(k9,lu)/dl
                           PARTL(k9,IEXP,2) = PARTL(k9,IEXP,2)
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
                     ixm = KSEQ(ixl,3)
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
               cnr(jk,jj) = -.5*PART(jk,jj,2)/PART(jk,jj,1)
               IF ( INNR.NE.0 ) CNOR(jk,jj) = cnr(jk,jj)
            ENDDO
            IF ( INNR.NE.1 ) THEN
               d = 0.
               g = 0.
               DO jj1 = jj , NEXPT
                  IF ( LNORM(jj1).EQ.jj ) THEN
                     k = NDST(jj1)
                     DO jk = 1 , k
                        d = d + YNRM(jk,jj1)*PART(jk,jj1,1)*YNRM(jk,jj1)
                        g = g - .5*YNRM(jk,jj1)*PART(jk,jj1,2)
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
            IF ( MCFIX.EQ.0 ) CNOR(jk,jj) = CNOR1(jk,jj)
            Chilo = Chilo + PARTL(jk,jj,1)*LOG(CNOR(jk,jj))
     &              **2 + PARTL(jk,jj,2)*2.*LOG(CNOR(jk,jj))
            Chisq = Chisq + CNOR(jk,jj)*CNOR(jk,jj)*PART(jk,jj,1)
     &              + CNOR(jk,jj)*PART(jk,jj,2)
         ENDDO
      ENDDO
      Chisq = Chisq + sumpr
      Chilo = Chilo + sum3
      RETURN
99012 FORMAT (1X,1I2,2X,32(1E10.4,1X))
99013 FORMAT (6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,1F6.4,6X,1E9.4,6X,
     &        1E9.4,3X,1F6.1,5X,1F4.1,10X,1A4)
      END
      SUBROUTINE FAKP
      IMPLICIT NONE
      INTEGER*4 i , IP , IPI , k , KF , l
      REAL*8 PILOG , x
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      DO i = 1 , 26
         x = DBLE(IP(i))
         PILOG(i) = LOG(x)
      ENDDO
      DO l = 1 , 26
         KF(1,l) = 0
         KF(2,l) = 0
      ENDDO
      DO k = 3 , 101
         CALL PRIM(k-1)
         DO i = 1 , 26
            KF(k,i) = KF(k-1,i) + IPI(i)
         ENDDO
      ENDDO
      END
      SUBROUTINE PRIM(N)
      IMPLICIT NONE
      INTEGER*4 i , IP , IPI , KF , N , nni , nnk
      REAL*8 PILOG
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      nnk = N
      DO i = 1 , 26
         nni = nnk
         IPI(i) = 0
 50      nni = nni/IP(i)
         IF ( IP(i)*nni.EQ.nnk ) THEN
            IPI(i) = IPI(i) + 1
            nnk = nni
            GOTO 50
         ENDIF
      ENDDO
      END
      SUBROUTINE SEQ(Idr)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , CONV , DELTA , DIPOL , ega , egs , emax , 
     &       EN , ENDEC , ENZ , F , FP , GF , GKP , SPIN , spinf , 
     &       spini , TAU , twoi
      REAL*8 ZPOL
      INTEGER*4 idecay , Idr , indx , inx , inx1 , ir , is , ISO , 
     &          istr1 , istr2 , ITMA , j , js , jsave , k , KLEC , kpa , 
     &          KSEQ , l , la
      INTEGER*4 la1 , LAMDA , LAMMAX , ld , LDNUM , LEAD , LEADF , LP1 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , 
     &          LP6 , LP7 , LP8 , LP9
      INTEGER*4 m , m1 , m6 , MEM , mk , mule , mulm , MULTI , n , n1 , 
     &          NDIM , NMAX , NMAX1 , nob
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /CATLF / FP(4,500,3) , GKP(4,500,2) , KLEC(75)
      m6 = 0
      DO l = 1 , 6
         m6 = m6 + MULTI(l)
      ENDDO
      idecay = 0
      Idr = 0
      DO l = 1 , LP3
         KLEC(l) = 0
      ENDDO
      DO k = 1 , LP2
         DO j = 1 , 3
            DO l = 1 , 4
               FP(l,k,j) = 0.
               IF ( j.NE.3 ) GKP(l,k,j) = 0.
            ENDDO
            DELTA(k,j) = 0.
         ENDDO
      ENDDO
      DO n = 1 , NMAX
         TAU(n) = EN(n)
      ENDDO
      DO n = 1 , NMAX
         emax = 0.
         DO j = 1 , NMAX
            IF ( TAU(j).GE.emax ) THEN
               emax = TAU(j)
               jsave = j
            ENDIF
         ENDDO
         DO is = 1 , NMAX
            DO la = 1 , 8
               IF ( la.LE.3 .OR. la.EQ.7 .OR. la.EQ.8 ) THEN
                  ld = LDNUM(la,is)
                  IF ( ld.NE.0 ) THEN
                     DO ir = 1 , ld
                        m = LEADF(is,ir,la)
                        IF ( m.EQ.jsave .OR. is.EQ.jsave ) THEN
                           IF ( is.NE.jsave .OR. EN(m).LT.EN(is) ) THEN
                              IF ( m.NE.jsave .OR. EN(is).LT.EN(m) )
     &                             THEN
                                 indx = MEM(is,m,la)
                                 idecay = idecay + 1
                                 KSEQ(idecay,1) = m
                                 KSEQ(idecay,2) = is
                                 KSEQ(idecay,3) = indx
                                 KSEQ(idecay,4) = la + 10
                                 IF ( EN(m).LE.EN(is) ) THEN
                                    KSEQ(idecay,1) = is
                                    KSEQ(idecay,2) = m
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         TAU(jsave) = -1.
      ENDDO
      DO l = 1 , idecay
         istr1 = 0
         IF ( KSEQ(l,4).LT.10 ) GOTO 200
         istr2 = 0
         n = KSEQ(l,1)
         m = KSEQ(l,2)
         inx = KSEQ(l,3)
         la = KSEQ(l,4) - 10
         ega = EN(n) - EN(m)
         twoi = 1./SQRT(2.*SPIN(n)+1.)
         spini = SPIN(n) + .001
         spinf = SPIN(m) + .001
         egs = SQRT(ega)*twoi
         js = l + 1
         la1 = 0
         inx1 = 0
         DO j = js , idecay
            IF ( KSEQ(j,4).GE.10 ) THEN
               n1 = KSEQ(j,1)
               m1 = KSEQ(j,2)
               IF ( n1.EQ.n .AND. m1.EQ.m ) THEN
                  inx1 = KSEQ(j,3)
                  la1 = KSEQ(j,4) - 10
                  KSEQ(j,4) = KSEQ(j,4) - 10
               ENDIF
            ENDIF
         ENDDO
         KSEQ(l,4) = KSEQ(l,4) - 10
         Idr = Idr + 1
         mule = 0
         mulm = 0
         nob = 1
 50      IF ( la.LE.3 ) THEN
            IF ( la.EQ.1 ) THEN
               DELTA(Idr,1) = 399.05*ega*egs
               mule = 1
               istr1 = 1
            ELSEIF ( la.EQ.2 ) THEN
               DELTA(Idr,1) = 3.4928*egs*ega*ega
               mule = 2
               istr1 = 2
            ELSEIF ( la.EQ.3 ) THEN
               DELTA(Idr,1) = .02391*ega*ega*ega*egs
               mule = 3
               istr1 = 3
            ELSE
               GOTO 100
            ENDIF
            GOTO 150
         ENDIF
 100     la = la - 6
         IF ( la.EQ.2 ) THEN
            DELTA(Idr,2) = .0368*ega*ega*egs
            mulm = 2
            istr2 = 5
         ELSE
            DELTA(Idr,2) = 4.1952*ega*egs
            mulm = 1
            istr2 = 4
         ENDIF
 150     IF ( nob.NE.2 ) THEN
            IF ( mule.NE.1 ) THEN
               nob = nob + 1
               IF ( la.GT.3 ) inx1 = inx
               IF ( la1.NE.0 ) THEN
                  la = la1
                  GOTO 50
               ENDIF
            ENDIF
            inx1 = 0
         ENDIF
         DELTA(Idr,3) = DELTA(Idr,1)*DELTA(Idr,2)
         DELTA(Idr,1) = DELTA(Idr,1)*DELTA(Idr,1)
         DELTA(Idr,2) = DELTA(Idr,2)*DELTA(Idr,2)
         KSEQ(Idr,1) = inx
         KSEQ(Idr,2) = inx1
         KSEQ(Idr,3) = n
         KSEQ(Idr,4) = m
         IF ( inx.GT.m6 ) THEN
            KSEQ(Idr,2) = inx
            KSEQ(Idr,1) = 0
         ENDIF
         ENDEC(Idr) = EN(n) - EN(m)
         DO mk = 1 , 7 , 2
            kpa = mk/2 + 1
            k = mk - 1
            IF ( mule.GE.3 .OR. k.NE.6 ) THEN
               GKP(kpa,Idr,1) = GF(k,spini,spinf,mule)*DELTA(Idr,1)
     &                          *(1.+CONV(ega,istr1))
               GKP(kpa,Idr,2) = GF(k,spini,spinf,mulm)*DELTA(Idr,2)
     &                          *(1.+CONV(ega,istr2))
               FP(kpa,Idr,1) = F(k,spini,spinf,mule,mule)*DELTA(Idr,1)
               FP(kpa,Idr,3) = F(k,spini,spinf,mulm,mule)*DELTA(Idr,3)
               FP(kpa,Idr,2) = F(k,spini,spinf,mulm,mulm)*DELTA(Idr,2)
            ENDIF
         ENDDO
         DELTA(Idr,1) = DELTA(Idr,1)*(1.+CONV(ega,istr1))
         DELTA(Idr,2) = DELTA(Idr,2)*(CONV(ega,istr2)+1.)
         KLEC(n) = KLEC(n) + 1
 200  ENDDO
      NMAX1 = 0
      DO n = 1 , NMAX
         IF ( KLEC(n).NE.0 ) NMAX1 = NMAX1 + 1
      ENDDO
      END
      FUNCTION GF(K,Sji,Sjf,L)
      IMPLICIT NONE
      REAL*8 GF , phase , Sjf , Sji , WSIXJ
      INTEGER*4 i , ix , jfz , jiz , K , kz , L , lz
      GF = 0.
      IF ( L.EQ.0 ) RETURN
      ix = INT(Sji+Sjf+.0001)
      i = ix + L + K
      phase = 1.
      IF ( i/2*2.NE.i ) phase = -1.
      kz = K*2
      jiz = Sji*2
      jfz = Sjf*2
      lz = L*2
      GF = phase*SQRT((jiz+1.)*(jfz+1.))*WSIXJ(jiz,jiz,kz,jfz,jfz,lz)
      END
      FUNCTION F(K,Sji,Sjf,L1,L2)
      IMPLICIT NONE
      REAL*8 F , phase , Sjf , Sji , WSIXJ , WTHREJ
      INTEGER*4 ix , jfz , jiz , K , kz , l , L1 , l1z , L2 , l2z
      F = 0.
      IF ( (L1*L2).EQ.0 ) RETURN
      ix = INT(Sji+Sjf+.0001)
      l = ix - 1
      phase = 1.
      IF ( l/2*2.NE.l ) phase = -1.
      kz = K*2
      jiz = Sji*2
      jfz = Sjf*2
      l1z = L1*2
      l2z = L2*2
      F = phase*SQRT((l1z+1.)*(l2z+1.)*(jiz+1.)*(kz+1.))
     &    *WTHREJ(l1z,l2z,kz,2,-2,0)*WSIXJ(jiz,jiz,kz,l2z,l1z,jfz)
      END
      FUNCTION CONV(Ega,N)
      IMPLICIT NONE
      REAL*8 AGELI , CC , CONV , cpo , cpo1 , cv , EG , Ega , Q
      INTEGER*4 j , N , n1 , NANG , nen , NICC
      DIMENSION cpo(51) , cpo1(51)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      IF ( N.EQ.0 ) THEN
         CONV = 0.0
      ELSEIF ( ABS(CC(1,N)).LT.1.E-9 ) THEN
         CONV = 0.0
      ELSE
         nen = 4
         DO j = 1 , NICC
            IF ( Ega.LE.EG(j) ) GOTO 50
         ENDDO
 50      n1 = j - 2
         IF ( n1.LT.1 ) n1 = 1
         IF ( (j+1).GT.NICC ) n1 = n1 - 1
         IF ( NICC.LE.4 ) THEN
            n1 = 1
            nen = NICC
         ENDIF
         DO j = 1 , nen
            cpo(j) = CC(n1+j-1,N)
            cpo1(j) = EG(n1+j-1)
         ENDDO
         CALL LAGRAN(cpo1,cpo,4,1,Ega,cv,2,1)
         CONV = cv
         RETURN
      ENDIF
      END
      FUNCTION WTHREJ(J1,J2,J3,M1,M2,M3)
      IMPLICIT NONE
      INTEGER*4 IP , IPI , iz , iza , izb , izc , izd , ize , izexp , 
     &          izf , izmax , izmin , J1 , J2 , J3 , jabc , jabm , 
     &          jbma , jj1 , jj2
      INTEGER*4 jj3 , jjha , jjhb , jjhc , jjhd , jlp , jma , jmax , 
     &          jmb , jmc , jmd , jme , jmf , jta , jtb , jtc , jvo , 
     &          jvora , KF , M1
      INTEGER*4 M2 , M3 , mm1 , mm2 , mm3 , n , nmax
      REAL*8 PILOG , qsumlo , sumlo , vorz , WTHREJ , wthrep , zuthre
      DIMENSION jvora(26)
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      wthrep = 0.E+00
      jjha = (J1+J2-J3)/2 + 1
      jjhb = (J1-J2+J3)/2 + 1
      jjhc = (-J1+J2+J3)/2 + 1
      IF ( (jjha.LT.1) .OR. (jjhb.LT.1) .OR. (jjhc.LT.1) .OR. 
     &     ((M1+M2+M3).NE.0) ) THEN
         WTHREJ = wthrep
         GOTO 99999
      ELSE
         jjhd = (J1+J2+J3+4)/2
         jmax = MAX(J1,J2,J3)
         IF ( jmax.NE.J1 ) THEN
            IF ( jmax.EQ.J2 ) THEN
               jj1 = J3
               jj2 = J1
               jj3 = J2
               mm1 = M3
               mm2 = M1
               mm3 = M2
               GOTO 100
            ELSEIF ( jmax.EQ.J3 ) THEN
               jj1 = J1
               jj2 = J2
               jj3 = J3
               mm1 = M1
               mm2 = M2
               mm3 = M3
               GOTO 100
            ENDIF
         ENDIF
         jj1 = J2
         jj2 = J3
         jj3 = J1
         mm1 = M2
         mm2 = M3
         mm3 = M1
      ENDIF
 100  jma = (jj1+mm1)/2
      jmb = (jj1-mm1)/2
      jmc = (jj2+mm2)/2
      jmd = (jj2-mm2)/2
      jme = (jj3+mm3)/2
      jmf = (jj3-mm3)/2
      jabc = (jj1+jj2-jj3)/2
      jabm = (jj2-jj3-mm1)/2
      jbma = (jj1+mm2-jj3)/2
      izmin = MAX(jabm,jbma,0)
      izmax = MIN(jabc,jmb,jmc)
      nmax = MAX(jjhd,izmax+1)
      DO n = 1 , 26
         IF ( IP(n).GE.nmax ) GOTO 200
      ENDDO
      WTHREJ = wthrep
      GOTO 99999
 200  DO jlp = 1 , n
         jta = KF(jjha,jlp) + KF(jjhb,jlp) + KF(jjhc,jlp) - KF(jjhd,jlp)
         jtb = KF(jma+1,jlp) + KF(jmb+1,jlp) + KF(jmc+1,jlp)
         jtc = KF(jmd+1,jlp) + KF(jme+1,jlp) + KF(jmf+1,jlp)
         jvora(jlp) = jta + jtb + jtc
      ENDDO
      vorz = -1.E+00
      IF ( 2*(izmin/2).EQ.izmin ) vorz = +1.E+00
      IF ( izmin.LE.izmax ) THEN
         DO iz = izmin , izmax
            qsumlo = 0.E+00
            iza = iz + 1
            izb = jabc + 1 - iz
            izc = jmb + 1 - iz
            izd = jmc + 1 - iz
            ize = iz - jabm + 1
            izf = iz - jbma + 1
            DO jlp = 1 , n
               izexp = jvora(jlp) - 2*KF(iza,jlp) - 2*KF(izb,jlp)
     &                 - 2*KF(izc,jlp) - 2*KF(izd,jlp) - 2*KF(ize,jlp)
     &                 - 2*KF(izf,jlp)
               sumlo = izexp
               qsumlo = qsumlo + sumlo*PILOG(jlp)*(.5E+00)
            ENDDO
            zuthre = vorz*EXP(qsumlo)
            wthrep = wthrep + zuthre
            vorz = -vorz
         ENDDO
         jvo = jj1 - jj2 - mm3
         IF ( 4*(jvo/4).NE.jvo ) wthrep = -wthrep
      ENDIF
      WTHREJ = wthrep
99999 END
      FUNCTION WSIXJ(J1,J2,J3,L1,L2,L3)
      IMPLICIT NONE
      INTEGER*4 IP , IPI , irj , irl , isa , isb , isc , isumfa , iva , 
     &          ivb , ivc , ivd , ivorfa , iz , iza , izb , izc , izd , 
     &          ize , izf
      INTEGER*4 izg , izh , izmax , izmin , J1 , J2 , J3 , KF , kqa , 
     &          kqb , kqc , kqd , kra , krb , krc , krd , ksa , ksb , 
     &          ksc , ksd
      INTEGER*4 kta , ktb , ktc , ktd , kua , kub , kuc , L1 , L2 , L3 , 
     &          n , nmax
      REAL*8 PILOG , qsumfa , qsumlo , sumlo , vorz , WSIXJ , wsixp , 
     &       zusix
      DIMENSION isumfa(26) , ivorfa(26)
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      wsixp = 0.E+00
      IF ( ((J1+J2-J3).GE.0) .AND. ((J1-J2+J3).GE.0) .AND. 
     &     ((-J1+J2+J3).GE.0) ) THEN
         IF ( ((J1+L2-L3).GE.0) .AND. ((J1-L2+L3).GE.0) .AND. 
     &        ((-J1+L2+L3).GE.0) ) THEN
            IF ( ((L1+J2-L3).GE.0) .AND. ((L1-J2+L3).GE.0) .AND. 
     &           ((-L1+J2+L3).GE.0) ) THEN
               IF ( ((L1+L2-J3).GE.0) .AND. ((L1-L2+J3).GE.0) .AND. 
     &              ((-L1+L2+J3).GE.0) ) THEN
                  kqa = (J1+J2-J3)/2
                  kqb = (J1-J2+J3)/2
                  kqc = (J2+J3-J1)/2
                  kqd = (J1+J2+J3)/2
                  kra = (J1+L2-L3)/2
                  krb = (J1-L2+L3)/2
                  krc = (L2+L3-J1)/2
                  krd = (J1+L2+L3)/2
                  ksa = (L1+J2-L3)/2
                  ksb = (L1-J2+L3)/2
                  ksc = (J2+L3-L1)/2
                  ksd = (L1+J2+L3)/2
                  kta = (L1+L2-J3)/2
                  ktb = (L1-L2+J3)/2
                  ktc = (L2+J3-L1)/2
                  ktd = (L1+L2+J3)/2
                  izmin = MAX(kqd,krd,ksd,ktd)
                  kua = kqa + kta + J3
                  kub = ksc + ktc + L1
                  kuc = krb + ktb + L2
                  izmax = MIN(kua,kub,kuc)
                  IF ( izmin.LE.izmax ) THEN
                     nmax = MAX(izmax+2,kqd+2,krd+2,ksd+2,ktd+2)
                     DO n = 1 , 26
                        IF ( IP(n).GE.nmax ) GOTO 5
                     ENDDO
                  ENDIF
                  GOTO 100
 5                vorz = -1.E+00
                  IF ( 2*(izmin/2).EQ.izmin ) vorz = +1.E+00
                  DO irl = 1 , n
                     iva = KF(kqa+1,irl) + KF(kqb+1,irl) + KF(kqc+1,irl)
     &                     - KF(kqd+2,irl)
                     ivb = KF(kra+1,irl) + KF(krb+1,irl) + KF(krc+1,irl)
     &                     - KF(krd+2,irl)
                     ivc = KF(ksa+1,irl) + KF(ksb+1,irl) + KF(ksc+1,irl)
     &                     - KF(ksd+2,irl)
                     ivd = KF(kta+1,irl) + KF(ktb+1,irl) + KF(ktc+1,irl)
     &                     - KF(ktd+2,irl)
                     ivorfa(irl) = iva + ivb + ivc + ivd
                  ENDDO
                  DO iz = izmin , izmax
                     sumlo = 0.E+00
                     iza = iz + 2
                     izb = iz - kqd + 1
                     izc = iz - krd + 1
                     izd = iz - ksd + 1
                     ize = iz - ktd + 1
                     izf = kua - iz + 1
                     izg = kub - iz + 1
                     izh = kuc - iz + 1
                     DO irj = 1 , n
                        isa = 2*KF(iza,irj) - 2*KF(izb,irj)
     &                        - 2*KF(izc,irj)
                        isb = -2*KF(izd,irj) - 2*KF(ize,irj)
     &                        - 2*KF(izf,irj)
                        isc = ivorfa(irj) - 2*KF(izg,irj)
     &                        - 2*KF(izh,irj)
                        isumfa(irj) = isa + isb + isc
                        qsumfa = isumfa(irj)
                        sumlo = sumlo + qsumfa*PILOG(irj)
                     ENDDO
                     qsumlo = (.5E+00)*sumlo
                     zusix = EXP(qsumlo)*vorz
                     wsixp = wsixp + zusix
                     vorz = -vorz
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 100  WSIXJ = wsixp
      END
      SUBROUTINE LAGRAN(X,Y,Ndata,Ipc,Xx,Yy,Iscal,Irc)
      IMPLICIT NONE
      REAL*8 arh , FUNC , FUNC1 , t , w , X , Xx , Y , y1 , Yy
      INTEGER*4 i , Ipc , Irc , Iscal , j , Ndata
      DIMENSION X(51) , Y(51) , w(51) , arh(51,51)
      IF ( Irc.EQ.2 ) THEN
      ELSEIF ( Irc.EQ.3 ) THEN
         DO i = 1 , Ndata
            t = 1.
            DO j = 1 , Ndata
               IF ( i.NE.j ) t = t*(Xx-X(j))/(X(i)-X(j))
            ENDDO
            arh(Ipc,i) = t
         ENDDO
         GOTO 100
      ELSEIF ( Irc.EQ.4 ) THEN
         GOTO 100
      ELSE
         DO i = 1 , Ndata
            w(i) = 1.
            DO j = 1 , Ndata
               IF ( i.NE.j ) w(i) = w(i)*(Xx-X(j))/(X(i)-X(j))
            ENDDO
         ENDDO
      ENDIF
      Yy = 0.
      DO j = 1 , Ndata
         y1 = Y(j)
         Yy = Yy + w(j)*FUNC(y1,Iscal)
      ENDDO
      Yy = FUNC1(Yy,Iscal)
      RETURN
 100  Yy = 0.
      DO j = 1 , Ndata
         y1 = Y(j)
         Yy = Yy + arh(Ipc,j)*FUNC(y1,Iscal)
      ENDDO
      Yy = FUNC1(Yy,Iscal)
      END
      FUNCTION FUNC(Y,I)
      IMPLICIT NONE
      REAL*8 FUNC , Y
      INTEGER*4 I
      IF ( I.EQ.2 ) THEN
         IF ( Y.LT.1.E-12 ) Y = 1.E-12
         FUNC = LOG(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC = SQRT(Y)
         GOTO 99999
      ENDIF
      FUNC = Y
      RETURN
99999 END
      FUNCTION FUNC1(Y,I)
      IMPLICIT NONE
      REAL*8 FUNC1 , Y
      INTEGER*4 I
      IF ( I.EQ.2 ) THEN
         FUNC1 = EXP(Y)
         RETURN
      ELSEIF ( I.EQ.3 ) THEN
         FUNC1 = Y*Y
         GOTO 99999
      ENDIF
      FUNC1 = Y
      RETURN
99999 END
      SUBROUTINE GKVAC(Il)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AKS , AVJI , beta , BETAR , DIPOL , DQ , 
     &       EN , EP , EPS , EROOT , FIEL , FIEX , GAMMA , GFAC , GKI , 
     &       POWER , QCEN , sp
      REAL*8 SPIN , SUM , TAU , time , TIMEC , TLBDG , VACDP , VINF , 
     &       XA , XA1 , XLAMB , XNOR , ZPOL
      INTEGER*4 i , IAXS , IBYP , IEXP , Il , ISO , ITTE , IZ , IZ1 , 
     &          KSEQ , NEXPT
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /BREC  / BETAR(50)
      COMMON /GGG   / AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /GVAC  / GKI(3) , SUM(3)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /THTAR / ITTE(50)
      IF ( ABS(XLAMB).GE.1.E-9 ) THEN
         IF ( ITTE(IEXP).EQ.0 ) THEN
            sp = SPIN(Il)
            beta = BETAR(IEXP)
            time = TAU(Il)
            CALL GKK(IZ,beta,sp,time,Il)
            VACDP(1,Il) = GKI(1)
            VACDP(2,Il) = GKI(2)
            VACDP(3,Il) = GKI(3)
            GOTO 99999
         ENDIF
      ENDIF
      DO i = 1 , 3
         VACDP(i,Il) = 1.
      ENDDO
      RETURN
99999 END
      SUBROUTINE GKK(Iz,Beta,Spin,Time,Il)
      IMPLICIT NONE
      REAL*8 AKS , alp , ATS , AVJI , Beta , ccf , down , DQ , dwc , f , 
     &       FIEL , GAMMA , GFAC , GKI , hmean , POWER , QCEN , rk , 
     &       sm , Spin
      REAL*8 SUM , Time , TIMEC , up , upc , VACDP , valmi , w2 , wrt , 
     &       WSIXJ , wsp , xji , xlam , XLAMB , XNOR
      INTEGER*4 i , IBYP , if2 , ifq , Il , imean , inq , irk2 , 
     &          ispin2 , ixji2 , Iz , j , k , k1 , k2 , l , m , ncoup , 
     &          nz
      COMMON /GVAC  / GKI(3) , SUM(3)
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      COMMON /GGG   / AVJI , GAMMA , XLAMB , TIMEC , GFAC , FIEL , POWER
      IF ( IBYP.NE.1 ) THEN
         imean = 0
         CALL XSTATIC(Iz,inq,ifq,Beta)
         l = 0
         DO i = 1 , 6
            AKS(i,Il) = 0.
         ENDDO
 50      IF ( imean.EQ.1 ) inq = 1
         IF ( imean.EQ.1 ) ifq = 1
         DO j = inq , ifq
            l = l + 1
            nz = Iz - j
            xji = ATS(nz)
            sm = Spin
            IF ( imean.EQ.1 ) xji = AVJI
            IF ( Spin.GT.xji ) sm = xji
            ncoup = INT(2.*sm+.5) + 1
            SUM(1) = 0.
            SUM(2) = 0.
            SUM(3) = 0.
            valmi = Spin - xji
            IF ( valmi.LT.0. ) valmi = -valmi
            DO m = 1 , ncoup
               f = valmi + DBLE(m) - 1.
               DO k = 1 , 3
                  rk = 2.*DBLE(k)
                  if2 = f*2. + 0.0001
                  irk2 = rk*2. + 0.0001
                  ispin2 = Spin*2. + 0.0001
                  ixji2 = xji*2. + 0.0001
                  SUM(k) = SUM(k)
     &                     + ((2.*f+1.)*WSIXJ(if2,if2,irk2,ispin2,
     &                     ispin2,ixji2))**2/(2.*xji+1.)
               ENDDO
            ENDDO
            IF ( imean.NE.1 ) THEN
               DO k = 1 , 3
                  k1 = 2*k - 1
                  AKS(k1,Il) = AKS(k1,Il) + SUM(k)
     &                         *EXP(-((QCEN-DBLE(j))/DQ)**2/2.)/XNOR
               ENDDO
               IF ( imean.EQ.0 ) GOTO 100
            ENDIF
            DO k = 1 , 3
               k1 = 2*k
               AKS(k1,Il) = AKS(k1,Il) + SUM(k)
            ENDDO
 100     ENDDO
         imean = imean + 1
         IF ( imean.EQ.1 ) GOTO 50
      ENDIF
      hmean = FIEL*Iz*(Beta**POWER)
      wsp = 4789.*GFAC*hmean/AVJI
      wsp = wsp*TIMEC
      wsp = wsp*wsp*AVJI*(AVJI+1.)/3.
      DO k = 1 , 3
         k2 = 2*k
         k1 = 2*k - 1
         wrt = wsp*k2*(k2+1)
         w2 = wrt
         wrt = -wrt/(1.-AKS(k2,Il))
         xlam = (1.-AKS(k2,Il))*(1.-EXP(wrt))/TIMEC
         up = (GAMMA*Time*AKS(k1,Il)+1.)/(Time*GAMMA+1.)
         up = up*XLAMB*Time + 1.
         down = Time*(xlam+XLAMB) + 1.
         GKI(k) = up/down
         alp = 9.*xlam*xlam + 8.*xlam*TIMEC*(w2-xlam*xlam)
         alp = SQRT(alp) - 3.*xlam
         alp = alp/4./xlam/TIMEC
         upc = xlam*Time*(down-2.*alp*alp*Time*TIMEC)
         dwc = (down+alp*Time)*(down+2.*alp*Time)
         ccf = 1. + upc/dwc
         GKI(k) = GKI(k)*ccf
      ENDDO
      END
      SUBROUTINE XSTATIC(Iz,Ido,Iup,Beta)
      IMPLICIT NONE
      REAL*8 AKS , Beta , DQ , h , QCEN , VACDP , XNOR
      INTEGER*4 IBYP , Ido , Iup , Iz , lq
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      h = 1./(1.+(Iz**.45*.012008/Beta)**1.666667)
      QCEN = Iz*h**.6
      DQ = SQRT(QCEN*(1.-h))/2.
      Iup = INT(QCEN+3.*DQ+.5)
      Ido = INT(QCEN-3.*DQ-.5)
      IF ( Iup.GT.Iz ) Iup = Iz
      IF ( Ido.LT.1 ) Ido = 1
      XNOR = 0.
      DO lq = Ido , Iup
         XNOR = XNOR + EXP(-((QCEN-DBLE(lq))/DQ)**2/2.)
      ENDDO
      END
      FUNCTION ATS(N)
      IMPLICIT NONE
      REAL*8 ATS , x , xm
      INTEGER*4 m , N
      IF ( N.LE.0 .OR. N.GT.96 ) THEN
         ATS = 0.
         RETURN
      ELSE
         x = N/2. + 1
         m = N/2 + 1
         xm = DBLE(m)
         IF ( ABS(x-xm).GE.1.E-9 ) THEN
            IF ( m.EQ.1 .OR. m.EQ.2 .OR. m.EQ.3 .OR. m.EQ.6 .OR. 
     &           m.EQ.7 .OR. m.EQ.10 .OR. m.EQ.15 .OR. m.EQ.16 .OR. 
     &           m.EQ.19 .OR. m.EQ.24 .OR. m.EQ.25 .OR. m.EQ.28 .OR. 
     &           m.EQ.31 .OR. m.EQ.35 .OR. m.EQ.37 .OR. m.EQ.40 .OR. 
     &           m.EQ.41 .OR. m.EQ.44 ) THEN
               ATS = .5
               RETURN
            ELSEIF ( m.EQ.4 .OR. m.EQ.5 .OR. m.EQ.8 .OR. m.EQ.9 .OR. 
     &               m.EQ.11 .OR. m.EQ.17 .OR. m.EQ.18 .OR. m.EQ.20 .OR. 
     &               m.EQ.26 .OR. m.EQ.27 .OR. m.EQ.36 .OR. m.EQ.42 .OR. 
     &               m.EQ.43 .OR. m.EQ.45 ) THEN
               ATS = 1.5
               RETURN
            ELSEIF ( m.EQ.12 .OR. m.EQ.14 .OR. m.EQ.21 .OR. m.EQ.23 .OR. 
     &               m.EQ.32 .OR. m.EQ.39 ) THEN
               ATS = 2.5
               RETURN
            ELSEIF ( m.EQ.13 .OR. m.EQ.22 .OR. m.EQ.38 ) THEN
               ATS = 4.5
               RETURN
            ELSEIF ( m.EQ.29 .OR. m.EQ.30 .OR. m.EQ.48 ) THEN
               ATS = 3.5
               RETURN
            ELSEIF ( m.EQ.33 ) THEN
               ATS = 7.5
               RETURN
            ELSEIF ( m.EQ.34 ) THEN
               ATS = 6.5
               GOTO 99999
            ELSEIF ( m.EQ.46 .OR. m.EQ.47 ) THEN
               ATS = 5.5
               RETURN
            ENDIF
         ENDIF
         m = m - 1
         IF ( m.EQ.4 .OR. m.EQ.8 .OR. m.EQ.17 .OR. m.EQ.26 .OR. 
     &        m.EQ.28 .OR. m.EQ.30 .OR. m.EQ.32 .OR. m.EQ.42 .OR. 
     &        m.EQ.45 .OR. m.EQ.48 ) THEN
            ATS = 2.
            RETURN
         ELSEIF ( m.EQ.10 .OR. m.EQ.36 ) THEN
         ELSEIF ( m.EQ.12 .OR. m.EQ.21 .OR. m.EQ.37 ) THEN
            ATS = 3.
            RETURN
         ELSEIF ( m.EQ.13 .OR. m.EQ.22 .OR. m.EQ.29 .OR. m.EQ.31 .OR. 
     &            m.EQ.34 .OR. m.EQ.38 .OR. m.EQ.47 ) THEN
            ATS = 4.
            RETURN
         ELSEIF ( m.EQ.33 ) THEN
            ATS = 8.
            RETURN
         ELSEIF ( m.EQ.46 ) THEN
            ATS = 6.
            RETURN
         ELSE
            ATS = 0.
            RETURN
         ENDIF
      ENDIF
      ATS = 1.
      RETURN
99999 END
      SUBROUTINE YLM(Theta,Ylmr)
      IMPLICIT NONE
      REAL*8 ct , ctsq , EPS , EROOT , FIEX , st , Theta , Ylmr
      INTEGER*4 i , IAXS , IEXP , j , l , lf , m
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      DIMENSION Ylmr(9,9) , st(7)
      ct = COS(Theta)
      ctsq = ct*ct
      IF ( IAXS(IEXP).EQ.0 ) THEN
         Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
         Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
         Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
         GOTO 99999
      ENDIF
      st(1) = SIN(Theta)
      DO i = 2 , 7
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      Ylmr(1,3) = .1089659406
      Ylmr(1,2) = -.2179318812*ct
      Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
      Ylmr(2,5) = .1248361677
      Ylmr(2,4) = -.3530900028*ct
      Ylmr(2,3) = .0943672726*(7.*ctsq-1.)
      Ylmr(2,2) = -.1334554768*ct*(7.*ctsq-3.)
      Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
      Ylmr(3,7) = .1362755124
      Ylmr(3,6) = -.4720722226*ct
      Ylmr(3,5) = .100646136*(11.*ctsq-1.)
      Ylmr(3,4) = -.1837538634*ct*(11.*ctsq-3.)
      Ylmr(3,3) = .0918769316*((33.*ctsq-18.)*ctsq+1.)
      Ylmr(3,2) = -.1162161475*ct*((33.*ctsq-30.)*ctsq+5.)
      Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
      DO l = 1 , 3
         lf = 2*l + 1
         DO m = 2 , lf
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)
         ENDDO
      ENDDO
      RETURN
99999 END
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
      SUBROUTINE ANGULA(Ygn,Idr,Iful,Fi0,Fi1,Trec,Gth,Figl,Ngl,Op2)
      IMPLICIT NONE
      REAL*8 AGELI , alab , arg , at , attl , BETAR , bt , CC , DELLA , 
     &       DELTA , dsig , DSIGS , EG , ENDEC , ENZ , EP , EPS , 
     &       EROOT , f , Fi0
      REAL*8 fi01 , Fi1 , fi11 , FIEX , Figl , FP , GKP , Gth , Q , qv , 
     &       sm , TAU , TETACM , TLBDG , Trec , TREP , ttx , VINF , XA , 
     &       XA1
      REAL*8 Ygn , ylmr , ZETA
      INTEGER*4 IAXS , Idr , IEXP , ifn , Iful , ig , il , inat , inx1 , 
     &          ipd , is , ITMA , ITTE , iu , ixs , IZ , IZ1 , j , ji , 
     &          jj
      INTEGER*4 jm , k , KLEC , kq , KSEQ , l , lf , lf1 , LZETA , 
     &          mind , NANG , NEXPT , Ngl , NICC , nlv
      CHARACTER*4 Op2
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
         nlv = KSEQ(l,3)
         il = (nlv-1)*28
         inx1 = KSEQ(l,2)
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
            Ygn(l) = alab(1,1)*.0795774715
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
            Ygn(l) = at(1)*.0795774715
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
               ipd = ITMA(IEXP,Ngl)
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,2*jj)*Q(2,ipd,2*jj)+Q(1,ipd,2*jj)*arg)
     &              /(Q(2,ipd,2*jj)+arg)
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ENDIF
      ENDDO
      IF ( Op2.EQ.'INTG' ) RETURN
      dsig = DSIGS(IEXP)
      ttx = TLBDG(IEXP)/57.2957795
      DO j = 1 , Idr
         Ygn(j) = Ygn(j)*dsig*SIN(ttx)
      ENDDO
      END
      SUBROUTINE READY(Idr,Ntap,Ipri)
      IMPLICIT NONE
      REAL*8 ap , CORF , DYEX , EP , TAU , TLBDG , u , UPL , VINF , w , 
     &       waga , XA , XA1 , xep , YEXP , YNRM , zp
      INTEGER*4 idc , idc1 , idcx , Idr , IDRN , ii , ILE , Ipri , IY , 
     &          iytot , iytt , IZ , IZ1 , j , k , kk , kkl , KSEQ , 
     &          lbg , LP1
      INTEGER*4 LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , 
     &          LP6 , LP7 , LP8 , LP9 , lxp , nanx , nde , nde1 , NDST , 
     &          ne , NEXPT , ns1
      INTEGER*4 ns2 , ns3 , ns4 , nsxh , nsyh , Ntap , nval , NYLDE
      DIMENSION iytot(32)
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CCCDS / NDST(50)
      REWIND Ntap
      DO k = 1 , LP6
         iytot(k) = 0
      ENDDO
      IF ( Ipri.EQ.1 ) WRITE (22,99001)
99001 FORMAT (5X/47X,'REPRINT OF EXPERIMENTAL DATA TO BE FITTED'//)
      DO lxp = 1 , NEXPT
         DO kkl = 1 , LP6
            NYLDE(lxp,kkl) = 0
         ENDDO
         ii = NDST(lxp)
         DO kk = 1 , ii
            READ (Ntap,*) ne , nanx , zp , ap , xep , nval , waga
            IF ( Ipri.EQ.1 ) WRITE (22,99002) ne , zp , ap , xep , 
     &                              NDST(ne) , waga
99002       FORMAT (1X,///5X,'EXPERIMENT',1X,1I2/2X,'PROJECTILE',1X,'(',
     &              1F4.0,',',1F4.0,')',1X,1F7.3,1X,'MEV',1X,'---',1I1,
     &              1X,'GE(LI) DETECTOR(S)',2X,'WEIGHT=',1E8.2/20X,
     &              '** EXPERIME','NTAL YIELDS **')
            IF ( Ipri.EQ.1 ) WRITE (22,99003)
99003       FORMAT (4X,'DECAY',1X,'IS',2X,'IF',1(9X,'YIELD+/-ERROR',9X)
     &              /)
            DO j = 1 , nval
               READ (Ntap,*) ns1 , ns2 , u , w
               nsxh = ns1
               nsyh = ns2
               IF ( ns1.GE.100 ) THEN
                  ns1 = ns1/100
                  ns2 = ns2/100
               ENDIF
               DO nde = 1 , Idr
                  IF ( ns1.EQ.KSEQ(nde,3) .AND. ns2.EQ.KSEQ(nde,4) )
     &                 GOTO 10
               ENDDO
               IF ( Ipri.EQ.1 ) WRITE (22,99005) ns1 , ns2
               GOTO 40
 10            idc = nde
               iytot(kk) = iytot(kk) + 1
               idc1 = 0
               IF ( nsxh.GE.100 ) THEN
                  ns3 = nsxh - 100*ns1
                  ns4 = nsyh - 100*ns2
                  DO nde1 = 1 , Idr
                     IF ( ns3.EQ.KSEQ(nde1,3) .AND. ns4.EQ.KSEQ(nde1,4)
     &                    ) GOTO 20
                  ENDDO
                  IF ( Ipri.EQ.1 ) WRITE (22,99005) ns3 , ns4
               ENDIF
               GOTO 30
 20            idcx = idc*1000 + nde1
               IF ( idc.GT.nde1 ) idcx = nde1*1000 + idc
               idc = idcx
 30            idc1 = idc
               IF ( idc1.GT.1000 ) idc1 = idc/1000
               IF ( Ipri.EQ.1 ) WRITE (22,99004) idc , KSEQ(idc1,3) , 
     &                                 KSEQ(idc1,4) , u , w
99004          FORMAT (2X,1I6,2X,1I2,2X,1I2,1(1E14.6,3X,1E14.6))
               iytt = iytot(kk)
               YEXP(kk,iytt) = u
               DYEX(kk,iytt) = w/(SQRT(waga)+1.E-4)
               IY(iytt,kk) = idc
 40         ENDDO
            iytt = iytot(kk)
            lbg = iytt - nval + 1
            CALL SZEREG(lbg,iytt,kk)
            NYLDE(lxp,kk) = nval
         ENDDO
      ENDDO
99005 FORMAT (1X///5X,'ERROR-NO MATRIX ELEMENT BETWEEN STATES',1X,1I2,
     &        ' AND ',1I2,/10X,'THIS TRANSITION IGNORED',//)
      END
      SUBROUTINE BRANR(Chisq,Nwyr,Chilo)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , BRAT , ch1 , ch2 , Chilo , Chisq , CONV , 
     &       DELTA , DIPOL , ELM , ELML , ELMU , EN , ENDEC , eng1 , 
     &       eng2 , ENZ , SA , SPIN
      REAL*8 TAU , u , ZPOL
      INTEGER*4 i1 , i2 , IBRC , iflg , iout , IPRM , ISO , ITMA , itt , 
     &          j1 , j2 , k , KSEQ , lab1 , lab2 , LAMDA , LAMMAX , 
     &          LDNUM , LEAD , mul2
      INTEGER*4 MULTI , n1 , n2 , NBRA , Nwyr
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /PRT   / IPRM(20)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      IF ( NBRA.EQ.0 ) RETURN
      IF ( IPRM(3).EQ.-1 ) WRITE (22,99001)
99001 FORMAT (1X,///10X,'EXP. AND CALCULATED BRANCHING RATIOS',//5X,
     &        'NS1',5X,'NF1',5X,'NS2',5X,'NF2',5X,'RATIO(1:2)',9X,
     &        'ERROR',7X,'CALC.RATIO',5X,'(EXP-CAL)/ERROR',//)
      Nwyr = Nwyr + NBRA
      mul2 = MULTI(1) + MULTI(2)
      DO k = 1 , NBRA
         ch1 = 0.
         ch2 = 0.
         iflg = 1
         itt = 1
         iout = 0
         n1 = IBRC(1,k)
         n2 = IBRC(2,k)
         i1 = KSEQ(n1,1)
         i2 = KSEQ(n2,1)
         eng1 = EN(KSEQ(n1,3)) - EN(KSEQ(n1,4))
         eng2 = EN(KSEQ(n2,3)) - EN(KSEQ(n2,4))
         IF ( i1.NE.0 ) THEN
            IF ( i1.LE.MULTI(1) ) lab1 = 1
            IF ( i1.GT.MULTI(1) .AND. i1.LE.mul2 ) lab1 = 2
            IF ( i1.GT.mul2 ) lab1 = 3
         ENDIF
         IF ( i2.NE.0 ) THEN
            IF ( i2.LE.MULTI(1) ) lab2 = 1
            IF ( i2.GT.MULTI(1) .AND. i2.LE.mul2 ) lab2 = 2
            IF ( i2.GT.mul2 ) lab2 = 3
         ENDIF
         IF ( i1.NE.0 ) ch1 = ELM(i1)*ELM(i1)*DELTA(n1,1)
     &                        /(1.+CONV(eng1,lab1))
         IF ( i2.NE.0 ) ch2 = ELM(i2)*ELM(i2)*DELTA(n2,1)
     &                        /(1.+CONV(eng2,lab2))
         j1 = KSEQ(n1,2)
         IF ( j1.NE.0 ) THEN
            iflg = iflg + 1
            lab1 = lab1 + 2
            ch1 = ch1 + ELM(j1)*ELM(j1)*DELTA(n1,2)/(1.+CONV(eng1,lab1))
         ENDIF
         j2 = KSEQ(n2,2)
         IF ( j2.NE.0 ) THEN
            iflg = iflg + 1
            lab2 = lab2 + 2
            ch2 = ch2 + ELM(j2)*ELM(j2)*DELTA(n2,2)/(1.+CONV(eng2,lab2))
         ENDIF
         u = (ch1/ch2-BRAT(k,1))/BRAT(k,2)
         Chisq = Chisq + u*u
         Chilo = Chilo + (BRAT(k,1)*LOG(ch1/ch2/BRAT(k,1))/BRAT(k,2))**2
         IF ( IPRM(3).EQ.-1 ) WRITE (22,99002) KSEQ(n1,3) , KSEQ(n1,4) , 
     &                               KSEQ(n2,3) , KSEQ(n2,4) , BRAT(k,1)
     &                               , BRAT(k,2) , ch1/ch2 , -u
99002    FORMAT (5X,3(1I2,6X),1I2,5X,3(1F10.5,5X),5X,1F4.1)
      ENDDO
      IF ( IPRM(3).EQ.-1 ) IPRM(3) = 0
      RETURN
      END
      SUBROUTINE LIMITS
      IMPLICIT NONE
      REAL*8 ELM , ELML , ELMU , SA
      INTEGER*4 IVAR , j , LMAXE , MAGEXC , MEMAX , MEMX6
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      DO j = 1 , MEMAX
         IF ( IVAR(j).NE.0 ) THEN
            IF ( ELM(j).GT.ELMU(j) .OR. ELM(j).LT.ELML(j) ) THEN
               IF ( ELM(j).GT.ELMU(j) ) THEN
                  ELM(j) = ELMU(j)
                  WRITE (22,99001) j , ELM(j)
               ELSE
                  ELM(j) = ELML(j)
                  WRITE (22,99001) j , ELM(j)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
99001 FORMAT (2X,'Warning - matrix element ',1I3,' reset to ',1F10.6)
      END
      SUBROUTINE SZEREG(Lst,Ls,L)
      IMPLICIT NONE
      REAL*8 CORF , DYEX , dyh , UPL , YEXP , yh , YNRM
      INTEGER*4 ia , ib , IDRN , ih , ILE , inx , IY , k , L , Ls , 
     &          lsp , Lst , lst1 , NYLDE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      IF ( Lst.EQ.Ls ) RETURN
      lst1 = Lst
      lsp = Ls - 1
 100  ia = IY(lst1,L)
      IF ( ia.GT.1000 ) ia = ia/1000
      inx = lst1
      DO k = lst1 , lsp
         ib = IY(k+1,L)
         IF ( ib.GT.1000 ) ib = ib/1000
         ia = MIN(ia,ib)
         IF ( ia.EQ.ib ) inx = k + 1
      ENDDO
      IF ( inx.NE.lst1 ) THEN
         ih = IY(lst1,L)
         IY(lst1,L) = IY(inx,L)
         IY(inx,L) = ih
         yh = YEXP(L,lst1)
         dyh = DYEX(L,lst1)
         YEXP(L,lst1) = YEXP(L,inx)
         DYEX(L,lst1) = DYEX(L,inx)
         YEXP(L,inx) = yh
         DYEX(L,inx) = dyh
      ENDIF
      lst1 = lst1 + 1
      IF ( lst1.GT.lsp ) RETURN
      GOTO 100
      END
      SUBROUTINE SIXEL(Rik,Rv,Em,Jk,Kk,Indx,Lu)
      IMPLICIT NONE
      REAL*8 a1 , al , al1 , c1 , c2 , DEV , Em , EPS , EROOT , FIEX , 
     &       Rik , rn , Rv , rx
      INTEGER*4 IAXS , IEXP , Indx , ITS , j , j1 , Jk , Kk , kk6 , 
     &          KVAR , l , l1 , Lu
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(600,7)
      COMMON /ODCH  / DEV(500)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /TRB   / ITS
      COMMON /SEL   / KVAR(500)
      kk6 = Kk + 5
      rn = DEV(Lu)
      al = (Rv-rn)*20./Rik
      IF ( ITS.EQ.1 .AND. KVAR(Indx).NE.0 ) WRITE (18,*) Lu , Indx , 
     &     IEXP , al/Em
      al1 = ABS(al)
      IF ( ITS.EQ.2 ) WRITE (18,*) Lu , Indx , IEXP , al1
      IF ( al1.LE.ABS(IMAG(ARM(kk6,Jk))) ) RETURN
      DO j = Kk , kk6
         a1 = ABS(IMAG(ARM(j,Jk)))
         IF ( al1.GT.a1 ) THEN
            j1 = j + 1
            DO l = j1 , kk6
               l1 = kk6 + j1 - l
               c1 = DBLE(ARM(l1-1,Jk))
               c2 = IMAG(ARM(l1-1,Jk))
               ARM(l1,Jk) = CMPLX(c1,c2)
            ENDDO
            rx = DBLE(Indx)
            ARM(j,Jk) = CMPLX(rx,al)
            GOTO 99999
         ENDIF
      ENDDO
99999 END
      SUBROUTINE PRELM(Iop)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , b , DIPOL , ELM , ELM25 , ELM26 , ELML , 
     &       ELMU , EN , pv , SA , SPIN , ste , ZPOL
      INTEGER*4 IBPS , inx , Iop , ISO , isp , IVAR , j , JZB , k , kk , 
     &          l , LAMDA , LAMMAX , LDNUM , LEAD , LMAXE , m , MAGEXC , 
     &          MEMAX , MEMX6
      INTEGER*4 MULTI , NDIM , NMAX , NMAX1
      CHARACTER*3 wrn
      COMMON /SWITCH/ JZB , IBPS
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /RESC  / ELM25(500) , ELM26(500)
      inx = 0
      WRITE (22,99001)
99001 FORMAT (2X/40X,'MATRIX ELEMENTS',//)
      DO j = 1 , 8
         m = MULTI(j)
         IF ( m.NE.0 ) THEN
            WRITE (22,99002) j
99002       FORMAT (5X,'MULTIPOLARITY=',1I1)
            IF ( Iop.EQ.1 ) WRITE (22,99003)
99003       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME')
            IF ( Iop.EQ.2 .OR. Iop.EQ.4 ) WRITE (22,99004)
99004       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',15X,'LIMITS')
            IF ( Iop.EQ.3 ) WRITE (22,99005)
99005       FORMAT (4X,'INDEX',3X,'NF',5X,'NS',10X,'ME',10X,'PC CHANGE',
     &              5X,'RED. TRANS. PROB.')
            DO k = 1 , NMAX
               l = LDNUM(j,k)
               IF ( l.NE.0 ) THEN
                  DO kk = 1 , l
                     inx = inx + 1
                     IF ( Iop.EQ.2 ) GOTO 2
                     IF ( Iop.EQ.3 ) THEN
                        isp = LEAD(2,inx)
                        pv = (ELMU(inx)-ELML(inx))/100.
                        wrn = '   '
                        IF ( (ELM(inx)-ELML(inx)).LT.pv ) wrn = '*?*'
                        IF ( (ELMU(inx)-ELM(inx)).LT.pv ) wrn = '*?*'
                        IF ( JZB.EQ.25 ) ste = ELM25(inx)
                        IF ( JZB.EQ.26 ) ste = ELM26(inx)
                        b = ELM(inx)*ELM(inx)/(2.*SPIN(isp)+1.)
                        IF ( LEAD(1,inx).EQ.LEAD(2,inx) ) b = 9999999.
                        WRITE (22,99009) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx) , 100.*(ELM(inx)-ste)/ste , 
     &                         b , wrn
                     ELSEIF ( Iop.EQ.4 ) THEN
                        GOTO 2
                     ELSE
                        WRITE (22,99006) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx)
99006                   FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5)
                     ENDIF
                     GOTO 5
 2                   IF ( IVAR(inx).EQ.0 ) THEN
                        WRITE (22,99007) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx)
99007                   FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,
     &                          'FIXED')
                     ELSEIF ( IVAR(inx).GT.1000 ) THEN
                        WRITE (22,99008) inx , LEAD(1,inx) , LEAD(2,inx)
     &                         , ELM(inx) , (IVAR(inx)-1000)
99008                   FORMAT (5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,
     &                          'COUPLED TO',1X,1I3)
                     ELSE
                        IF ( Iop.EQ.2 ) WRITE (22,99009) inx , 
     &                       LEAD(1,inx) , LEAD(2,inx) , ELM(inx) , 
     &                       ELML(inx) , ELMU(inx)
                        IF ( Iop.EQ.4 .AND. JZB.EQ.25 ) WRITE (22,99009)
     &                       inx , LEAD(1,inx) , LEAD(2,inx) , 
     &                       ELM25(inx) , ELML(inx) , ELMU(inx)
                        IF ( Iop.EQ.4 .AND. JZB.EQ.26 ) WRITE (22,99009)
     &                       inx , LEAD(1,inx) , LEAD(2,inx) , 
     &                       ELM26(inx) , ELML(inx) , ELMU(inx)
                     ENDIF
 5                ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO
99009 FORMAT (5X,1I3,5X,1I2,5X,1I2,3(5X,1F10.5),1A3)
      END
      SUBROUTINE RECOIL(Alab,Attl,Beta,Theta)
      IMPLICIT NONE
      REAL*8 Alab , atemp , Attl , Beta , betasq , dum , hold , test , 
     &       Theta
      INTEGER*4 i , i1 , j , l , m
      DIMENSION Alab(9,9) , Attl(9,9) , atemp(16)
      hold = Alab(1,1)
      IF ( ABS(hold).LT.1.E-9 ) RETURN
      CALL ROTATE(Alab,Attl,-Theta,7,2)
      Attl(2,1) = (2./SQRT(15.))*(SQRT(5.)*Attl(1,1)-Attl(3,1))
      Attl(2,2) = -Attl(3,2)/SQRT(5.)
      Attl(4,1) = (4./SQRT(35.))*(3.*Attl(3,1)-SQRT(5.)*Attl(5,1))
      Attl(4,2) = (8.*SQRT(2.)*Attl(3,2)-5.*SQRT(3.)*Attl(5,2))
     &            /SQRT(35.)
      Attl(4,3) = (2./SQRT(7.))*(2.*Attl(3,3)-SQRT(3.)*Attl(5,3))
      Attl(4,4) = -Attl(5,4)
      Attl(6,1) = (10./SQRT(11.))*(Attl(5,1)-(3./SQRT(13.))*Attl(7,1))
      Attl(6,2) = (1./SQRT(11.))
     &            *(4.*SQRT(6.)*Attl(5,2)-5.*SQRT(35./13.)*Attl(7,2))
      Attl(6,3) = SQRT(4./11.)
     &            *(SQRT(21.)*Attl(5,3)-10.*SQRT(2./13.)*Attl(7,3))
      Attl(6,4) = SQRT(1./11.)*(8.*Attl(5,4)-15.*SQRT(3./13.)*Attl(7,4))
      Attl(6,5) = SQRT(4./11.)*(3.*Attl(5,5)-5.*SQRT(5./13.)*Attl(7,5))
      Attl(6,6) = -Attl(7,6)*SQRT(25./13.)
      Attl(8,1) = (56./SQRT(195.))*Attl(7,1)
      Attl(8,2) = (32./SQRT(65.))*Attl(7,2)
      Attl(8,3) = (8.*SQRT(3./13.))*Attl(7,3)
      Attl(8,4) = (16.*SQRT(2./39.))*Attl(7,4)
      Attl(8,5) = (8.*SQRT(11./65.))*Attl(7,5)
      Attl(8,6) = (16.*SQRT(2./65.))*Attl(7,6)
      Attl(8,7) = (8./SQRT(15.))*Attl(7,7)
      DO l = 2 , 8 , 2
         DO m = 1 , l
            Attl(l,m) = Beta*Attl(l,m)
         ENDDO
      ENDDO
      betasq = Beta*Beta
      IF ( betasq.GE.1.0E-10 ) THEN
         i1 = 0
         DO i = 1 , 7 , 2
            DO j = 1 , i
               i1 = i1 + 1
               atemp(i1) = Attl(i,j)
            ENDDO
         ENDDO
         dum = (2./5.)*SQRT(5.)*atemp(1) - (10./7.)*atemp(2) + (12./35.)
     &         *SQRT(5.)*atemp(5)
         Attl(3,1) = atemp(2) + betasq*dum
         dum = -(17./14.)*atemp(3) + (2./7.)*SQRT(6.)*atemp(6)
         Attl(3,2) = atemp(3) + betasq*dum
         dum = -(4./7.)*atemp(4) + (2./7.)*SQRT(3.)*atemp(7)
         Attl(3,3) = atemp(4) + betasq*dum
         dum = (8./7.)*SQRT(5.)*atemp(2) - (380./77.)*atemp(5)
     &         + (100./11.)*SQRT(1./13.)*atemp(10)
         Attl(5,1) = atemp(5) + betasq*dum
         dum = (20./21.)*SQRT(6.)*atemp(3) - (723./154.)*atemp(6)
     &         + (20./11.)*SQRT(70./39.)*atemp(11)
         Attl(5,2) = atemp(6) + betasq*dum
         dum = (20./21.)*SQRT(3.)*atemp(4) - (306./77.)*atemp(7)
     &         + (40./11.)*SQRT(14./39.)*atemp(12)
         Attl(5,3) = atemp(7) + betasq*dum
         dum = -(61./22.)*atemp(8) + (40./11.)*SQRT(3./13.)*atemp(13)
         Attl(5,4) = atemp(8) + betasq*dum
         dum = -(12./11.)*atemp(9) + (20./11.)*SQRT(5./13.)*atemp(14)
         Attl(5,5) = atemp(9) + betasq*dum
         dum = (210./11.)*SQRT(1./13.)*atemp(5) - (574./55.)*atemp(10)
         Attl(7,1) = atemp(10) + betasq*dum
         dum = (14./11.)*SQRT(210./13.)*atemp(6) - (1121./110.)
     &         *atemp(11)
         Attl(7,2) = atemp(11) + betasq*dum
         dum = (28./11.)*SQRT(42./13.)*atemp(7) - (104./11.)*atemp(12)
         Attl(7,3) = atemp(12) + betasq*dum
         dum = (84./11.)*SQRT(3./13.)*atemp(8) - (181./22.)*atemp(13)
         Attl(7,4) = atemp(13) + betasq*dum
         dum = (42./11.)*SQRT(5./13.)*atemp(9) - (358./55.)*atemp(14)
         Attl(7,5) = atemp(14) + betasq*dum
         Attl(7,6) = atemp(15)*(1.-(43./10.)*betasq)
         Attl(7,7) = atemp(16)*(1.-(8./5.)*betasq)
         Attl(9,1) = (672./5.)*SQRT(1./221.)*atemp(10)*betasq
         Attl(9,2) = (144./5.)*SQRT(21./221.)*atemp(11)*betasq
         Attl(9,3) = 36.*SQRT(12./221.)*atemp(12)*betasq
         Attl(9,4) = 24.*SQRT(22./221.)*atemp(13)*betasq
         Attl(9,5) = (144./5.)*SQRT(11./221.)*atemp(14)*betasq
         Attl(9,6) = (72./5.)*SQRT(2./17.)*atemp(15)*betasq
         Attl(9,7) = (24./5.)*SQRT(7./17.)*atemp(16)*betasq
      ENDIF
      CALL ROTATE(Attl,Alab,Theta,9,1)
      test = ABS(1.0-Alab(1,1)/hold)
      IF ( test.GT.1.0E-07 ) THEN
         WRITE (22,99001) test
99001    FORMAT (' ERROR IN ROTATION',1X,1E10.3/)
      ENDIF
      RETURN
      END
      SUBROUTINE ROTATE(Alab,Attl,Theta,K2,Kd)
      IMPLICIT NONE
      REAL*8 Alab , Attl , djarg , DJMM , dkkk , sum , Theta
      INTEGER*4 idj , idm , idmp , j , k , K2 , ka , kappa , kapri , Kd
      DIMENSION Alab(9,9) , Attl(9,9)
      IF ( ABS(Theta).GT..01 ) THEN
         djarg = Theta
         DO ka = 1 , K2 , Kd
            idj = ka - 1
            DO kappa = 1 , ka
               idmp = kappa - 1
               sum = 0.0
               DO kapri = 1 , ka
                  idm = kapri - 1
                  dkkk = DJMM(djarg,idj,idm,idmp)
                  sum = sum + dkkk*Alab(ka,kapri)
               ENDDO
               IF ( ka.NE.1 ) THEN
                  DO kapri = 2 , ka
                     idm = -kapri + 1
                     dkkk = DJMM(djarg,idj,idm,idmp)
                     sum = sum + dkkk*Alab(ka,kapri)*(-1.0)**(kapri-1)
                  ENDDO
               ENDIF
               Attl(ka,kappa) = sum
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
      DO j = 1 , 9
         DO k = 1 , 9
            Attl(j,k) = Alab(j,k)
         ENDDO
      ENDDO
      RETURN
99999 END
      SUBROUTINE YLM1(Theta,Ylmr)
      IMPLICIT NONE
      REAL*8 ct , ctsq , st , Theta , Ylmr
      INTEGER*4 i , j , l , m
      DIMENSION Ylmr(9,9) , st(9)
      ct = COS(Theta)
      ctsq = ct*ct
      st(1) = SIN(Theta)
      DO i = 2 , 9
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      DO l = 2 , 9
         DO m = 1 , 9
            Ylmr(l,m) = 0.0
         ENDDO
      ENDDO
      Ylmr(2,2) = -SQRT(6.)/2.
      Ylmr(2,1) = SQRT(3.)*ct
      Ylmr(3,3) = SQRT(30.)/4.
      Ylmr(3,2) = -(SQRT(30.)/2.)*ct
      Ylmr(3,1) = (SQRT(5.)/2.)*(3.*ctsq-1.)
      Ylmr(4,4) = -SQRT(35.)/4.
      Ylmr(4,3) = (SQRT(210.)/4.)*ct
      Ylmr(4,2) = -(SQRT(21.)/4.)*(5.*ctsq-1.)
      Ylmr(4,1) = (SQRT(7.)/2.)*ct*(5.*ctsq-3.)
      Ylmr(5,5) = 3.*SQRT(70.)/16.
      Ylmr(5,4) = -(3.*SQRT(35.)/4.)*ct
      Ylmr(5,3) = (3.*SQRT(10.)/8.)*(7.*ctsq-1.)
      Ylmr(5,2) = -(3.*SQRT(5.)/4.)*ct*(7.*ctsq-3.)
      Ylmr(5,1) = (3./8.)*((35.*ctsq-30.)*ctsq+3.)
      Ylmr(6,6) = -3.*SQRT(77.)/16.
      Ylmr(6,5) = (3.*SQRT(770.)/16.)*ct
      Ylmr(6,4) = -(SQRT(385.)/16.)*(9.*ctsq-1.)
      Ylmr(6,3) = (SQRT(2310.)/8.)*ct*(3.*ctsq-1.)
      Ylmr(6,2) = -(SQRT(330.)/16.)*((21.*ctsq-14.)*ctsq+1.)
      Ylmr(6,1) = (SQRT(11.)/8.)*ct*((63.*ctsq-70.)*ctsq+15.)
      Ylmr(7,7) = SQRT(3003.)/32.
      Ylmr(7,6) = -(3.*SQRT(1001.)/16.)*ct
      Ylmr(7,5) = (3.*SQRT(182.)/32.)*(11.*ctsq-1.)
      Ylmr(7,4) = -(SQRT(1365.)/16.)*ct*(11.*ctsq-3.)
      Ylmr(7,3) = (SQRT(1365.)/32.)*((33.*ctsq-18.)*ctsq+1.)
      Ylmr(7,2) = -(SQRT(546.)/16.)*ct*((33.*ctsq-30.)*ctsq+5.)
      Ylmr(7,1) = (SQRT(13.)/16.)*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
      Ylmr(8,8) = -3.*SQRT(1430.)/64.
      Ylmr(8,7) = (3.*SQRT(5005.)/32.)*ct
      Ylmr(8,6) = -(3.*SQRT(770.)/64.)*(13.*ctsq-1.)
      Ylmr(8,5) = (3.*SQRT(770.)/32.)*(13.*ctsq-3.)*ct
      Ylmr(8,4) = -(3.*SQRT(70.)/64.)*((143.*ctsq-66.)*ctsq+3.)
      Ylmr(8,3) = (3.*SQRT(35.)/32.)*((143.*ctsq-110.)*ctsq+15.)*ct
      Ylmr(8,2) = -(SQRT(210.)/64.)
     &            *(((429.*ctsq-495.)*ctsq+135.)*ctsq-5.)
      Ylmr(8,1) = (SQRT(15.)/16.)
     &            *(((429.*ctsq-693.)*ctsq+315.)*ctsq-35.)*ct
      Ylmr(9,9) = 3.*SQRT(24310.)/256.
      Ylmr(9,8) = -(3.*SQRT(24310.)/64.)*ct
      Ylmr(9,7) = (SQRT(7293.)/64.)*(15.*ctsq-1.)
      Ylmr(9,6) = -(3.*SQRT(34034.)/64.)*(5.*ctsq-1.)*ct
      Ylmr(9,5) = (3.*SQRT(2618.)/128.)*((65.*ctsq-26.)*ctsq+1.)
      Ylmr(9,4) = -(SQRT(39270.)/64.)*((39.*ctsq-26.)*ctsq+3.)*ct
      Ylmr(9,3) = (3.*SQRT(595.)/64.)
     &            *(((143.*ctsq-143.)*ctsq+33.)*ctsq-1.)
      Ylmr(9,2) = -(3.*SQRT(34.)/64.)
     &            *(((715.*ctsq-1001.)*ctsq+385.)*ctsq-35.)*ct
      Ylmr(9,1) = (SQRT(17.)/128.)
     &            *((((6435.*ctsq-12012.)*ctsq+6930.)*ctsq-1260.)
     &            *ctsq+35.)
      DO l = 2 , 9
         Ylmr(l,1) = Ylmr(l,1)*.0795774715
         DO m = 2 , l
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)*.0795774715
         ENDDO
      ENDDO
      END
      SUBROUTINE FIINT(Fi0,Fi1,At,Ixs)
      IMPLICIT NONE
      REAL*8 At , Fi0 , Fi1 , wsp
      INTEGER*4 Ixs , j , jf , js , m , mm
      DIMENSION At(28)
      IF ( Ixs.NE.0 ) THEN
         DO m = 2 , 7
            js = m/2
            mm = m - 1
            wsp = (SIN(mm*Fi1)-SIN(mm*Fi0))/mm
            js = js*7 + m
            jf = m + 21
            DO j = js , jf , 7
               At(j) = At(j)*wsp
            ENDDO
         ENDDO
         wsp = Fi1 - Fi0
      ENDIF
      IF ( Ixs.EQ.0 ) wsp = 6.283185308
      DO j = 1 , 4
         js = (j-1)*7 + 1
         At(js) = At(js)*wsp
      ENDDO
      END
      SUBROUTINE FIINT1(Fi0,Fi1,Alab,Ixs)
      IMPLICIT NONE
      REAL*8 Alab , Fi0 , Fi1 , wsp
      INTEGER*4 Ixs , j , m , mm
      DIMENSION Alab(9,9)
      IF ( Ixs.NE.0 ) THEN
         DO m = 2 , 9
            mm = m - 1
            wsp = (SIN(mm*Fi1)-SIN(mm*Fi0))/mm
            DO j = 1 , 9
               Alab(j,m) = Alab(j,m)*wsp
            ENDDO
         ENDDO
         wsp = Fi1 - Fi0
      ENDIF
      IF ( Ixs.EQ.0 ) wsp = 6.283185308
      DO j = 1 , 9
         Alab(j,1) = Alab(j,1)*wsp
      ENDDO
      END
      SUBROUTINE TAPMA(Lx,Iske,Isko,Iskf,Nflr,Idr,Nco,Nft,Enb)
      IMPLICIT NONE
      REAL*8 DS , DSE , DSG , emn , emx , en0 , Enb , tmn , tmx , tta , 
     &       XV , YGN , YGP , YV , ZETA , ZV
      INTEGER*4 Idr , IFMO , Iske , Iskf , Isko , j , jf , jj , js , k , 
     &          Lx , lx1 , LZETA , na , Nco , ne , nfil , nfilt , Nflr , 
     &          Nft
      INTEGER*4 ng , ng1 , ntt
      COMMON /VLIN  / XV(51) , YV(51) , ZV(20) , DSG(20) , DSE(20) , DS
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /YTEOR / YGN(500) , YGP(500) , IFMO
      Nft = 0
      nfilt = 0
      REWIND 14
      IF ( Iske.NE.0 ) THEN
 50      READ (14,*) ne , ntt , emn , emx , tmn , tmx , na , tmx , tmx , 
     &               tmx
         nfil = ne*ntt*na
         nfilt = nfilt + nfil
         DO j = 1 , nfil
            READ (14,*) lx1 , Enb , tta , ng , DS , (YGN(k),k=1,Idr)
         ENDDO
         IF ( nfilt.NE.Iske ) GOTO 50
      ENDIF
      IF ( Nco.EQ.0 ) RETURN
      READ (14,*) ne , ntt , emn , emx , tmn , tmx , na , tmx , tmx , 
     &            tmx
      IF ( Isko.NE.0 ) THEN
         DO j = 1 , Isko
            READ (14,*) lx1 , Enb , tta , ng , DS , (YGN(k),k=1,Idr)
         ENDDO
      ENDIF
      DO j = 1 , Nflr
         js = (j-1)*Idr + 1
         jf = js + Idr - 1
         READ (14,*) lx1 , Enb , tta , ng1 , DS , (ZETA(k),k=js,jf)
         IF ( lx1.NE.Lx ) Nft = 1
         IF ( Nft.EQ.1 ) GOTO 100
         XV(j) = tta/57.2957795
         IF ( Iskf.NE.0 .AND. j.NE.Nflr ) THEN
            DO jj = 1 , Iskf
               READ (14,*) lx1 , en0 , tta , ng , DS , (YGN(k),k=1,Idr)
            ENDDO
         ENDIF
      ENDDO
      RETURN
 100  WRITE (22,99001)
99001 FORMAT (10X///10X,'TAPE READ ERROR'/10X,'JOB ABORTED')
      END
      FUNCTION SIMIN(Np,H,Y)
      IMPLICIT NONE
      REAL*8 ee , H , SIMIN , sm , Y
      INTEGER*4 ik , in , Np
      DIMENSION Y(101)
      IF ( Np.GE.3 ) THEN
         ik = Np - 2
         sm = Y(1) + Y(Np)
         DO in = 1 , ik
            ee = in/2.
            sm = sm + 2.*Y(in+1)/(1.+INT(ee)-ee)
         ENDDO
         SIMIN = sm*H/3.
         RETURN
      ELSEIF ( Np.EQ.1 ) THEN
         SIMIN = Y(1)
         GOTO 99999
      ENDIF
      SIMIN = (Y(1)+Y(2))*H/2.
      RETURN
99999 END
      SUBROUTINE MIXUP
      IMPLICIT NONE
      REAL*8 ELM , ELML , ELMU , RNDM , SA , SE
      INTEGER*4 IVAR , k , k1 , LMAXE , MAGEXC , MEMAX , MEMX6
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /XRA   / SE
      DO k = 1 , MEMAX
         IF ( IVAR(k).NE.0 .AND. IVAR(k).LE.999 ) ELM(k) = ELML(k)
     &        + RNDM(SE)*(ELMU(k)-ELML(k))
      ENDDO
      DO k = 1 , MEMAX
         IF ( IVAR(k).GE.999 ) THEN
            k1 = IVAR(k) - 1000
            IF ( ABS(ELMU(k1)).LT.1.E-9 ) THEN
               ELM(k) = 0.
            ELSE
               ELM(k) = ELM(k1)*SA(k)
            ENDIF
         ENDIF
      ENDDO
      END
      FUNCTION FXIS1(I,N)
      IMPLICIT NONE
      REAL*8 FXIS1 , XI
      INTEGER*4 I , N
      COMMON /CXI   / XI(500)
      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS1 = 1.
         GOTO 99999
      ENDIF
      FXIS1 = -SIGN(1.,XI(I))
      RETURN
99999 END
      FUNCTION FXIS2(I,N)
      IMPLICIT NONE
      REAL*8 FXIS2 , XI
      INTEGER*4 I , N
      COMMON /CXI   / XI(500)
      IF ( N.EQ.2 .OR. N.EQ.3 .OR. N.EQ.5 .OR. N.EQ.6 ) THEN
         FXIS2 = -SIGN(1.,XI(I))
         GOTO 99999
      ENDIF
      FXIS2 = 1.
      RETURN
99999 END
      SUBROUTINE PODZIEL(I,J)
      IMPLICIT NONE
      INTEGER*4 I , IAPR , IDIVE , ISEX , J , k , l , l1 , l2 , LERF , 
     &          LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , 
     &          LP4 , LP6
      INTEGER*4 LP7 , LP8 , LP9
      REAL*8 QAPR
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /APRX  / LERF , IDIVE(50,2)
      IF ( I.NE.3 ) THEN
         IF ( I.EQ.1 ) THEN
            l1 = IDIVE(J,1)
            IDIVE(J,1) = l1 + 1
            GOTO 100
         ELSE
            l1 = IDIVE(J,2)
            IDIVE(J,2) = l1 + 1
         ENDIF
      ENDIF
      l2 = IDIVE(J,2)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2
         DO l = 1 , 7
            QAPR(k,2,l) = QAPR(k,2,l)*l1/l2
         ENDDO
      ENDDO
      IF ( I.EQ.2 ) WRITE (22,99001) J , IDIVE(J,1) , l2
      IF ( I.NE.3 ) RETURN
 100  l2 = IDIVE(J,1)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2
         DO l = 1 , 7
            QAPR(k,1,l) = QAPR(k,1,l)*l1/l2
         ENDDO
      ENDDO
      IF ( I.EQ.1 ) WRITE (22,99001) J , l2 , IDIVE(J,2)
      RETURN
99001 FORMAT (5X,'*****',1X,'EXP(A) EXPANSION FAILURE!',1X,'*****'/5X,
     &        'EXPERIMENT',1X,1I2,3X,'NEW SUBDIVISION',1X,'(',1I1,',',
     &        1I1,')')
      END
      SUBROUTINE KLOPOT(K,Rlr)
      IMPLICIT NONE
      REAL*8 a , al , al1 , b , c , ch , CORF , d , dy , DYEX , e , 
     &       ELM , ELML , ELMU , EP , g , g1 , g2 , rl , Rlr
      REAL*8 SA , sgm , TLBDG , u , umm , ump , UPL , ux , VINF , XA , 
     &       XA1 , YEXP , YNRM , ZETA
      INTEGER*4 i , IDRN , iex , iexh , iexp , ILE , indx , inh , ipf , 
     &          IVAR , IY , IZ , IZ1 , j , jm , jp , K , KVAR , l , lc
      INTEGER*4 ll , LMAXE , lngt , loc , LP1 , LP10 , LP11 , LP12 , 
     &          LP13 , LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &          lu , LZETA , MAGEXC
      INTEGER*4 MEMAX , MEMX6 , NEXPT , nf , ni , nm , np , NYLDE
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /SEL   / KVAR(500)
      REWIND 14
      REWIND 18
      ipf = 1
      lngt = 0
      indx = 1
      REWIND 15
      REWIND 17
      DO i = 1 , MEMAX
         ELM(i) = 0.
         ELMU(i) = 0.
         ELML(i) = 0.
      ENDDO
      iexh = 1
 100  g = 0.
      d = 0.
 200  READ (15,*) iex , a , b , c , e
      IF ( iex.NE.iexh ) THEN
         EP(iexh) = g/d
         TLBDG(iexh) = g
         VINF(iexh) = d
         iexh = iex
         BACKSPACE 15
         IF ( iex.NE.0 ) GOTO 100
         REWIND 15
         iexp = 1
      ELSE
         g = g + e*a/c/c
         d = d + a*a/c/c
         lngt = lngt + 1
         GOTO 200
      ENDIF
 300  g1 = 0.
      g2 = 0.
      inh = indx
      iexh = iexp
 400  READ (18,*) lu , indx , iexp , al
      IF ( indx.NE.0 ) THEN
         READ (17,*) ni , nf , sgm , al1
         READ (15,*) iex , a , b , c , e
         IF ( iexp.NE.1 .AND. ipf.NE.1 ) THEN
 420        READ (15,*) iex , a , b , c , e
            IF ( iexp.NE.iex ) GOTO 420
            ipf = 1
         ENDIF
         IF ( indx.EQ.inh ) THEN
            dy = al*al1/b/c
            g1 = e*dy + g1
            g2 = -2.*dy*a + g2
            WRITE (14,*) indx , iexp , ni , nf , dy , a , e , c
            GOTO 400
         ENDIF
      ENDIF
      loc = (iexh-1)*LP2 + inh
      ipf = 0
      ZETA(loc) = (VINF(iexh)*g1+TLBDG(iexh)*g2)/VINF(iexh)/VINF(iexh)
      inh = indx
      REWIND 15
      BACKSPACE 17
      BACKSPACE 18
      IF ( indx.NE.0 ) GOTO 300
      WRITE (14,*) indx , iexp , ni , nf , dy , a , e , b
      REWIND 14
      REWIND 17
 500  READ (14,*) indx , iexp , ni , nf , dy , a , e , b
      IF ( indx.EQ.0 ) THEN
         WRITE (17,*) indx , iexp , sgm , ni , nf , u
         REWIND 17
         ll = 0
         ch = 0.
 550     READ (17,*) indx , iexp , sgm , ni , nf , u
         IF ( indx.EQ.0 ) THEN
            WRITE (22,99001)
99001       FORMAT (2X////40X,'TROUBLESHOOTING ROUT',
     &              'INE HAS BEEN ACTIVATED...'//5X,
     &              'LOCAL MINIMUM ANALYSIS FOLLOWS:'//)
            WRITE (22,99002) ch/ll
99002       FORMAT (2X//5X,'CHISQ FOR FIRST GE(LI)S ONLY ',
     &              'WITH INDEPENDENT NORMALIZATION=',1E12.4//5X,
     &              'NORM.CONSTANTS:'//)
            DO i = 1 , NEXPT
               WRITE (22,99003) i , EP(i)
99003          FORMAT (5X,'EXP.',1X,1I2,5X,'C=',1E14.6)
            ENDDO
            WRITE (22,99004)
99004       FORMAT (1X//5X,'M.E.',20X,'RL',20X,'STRENGTH',//)
            DO i = 1 , MEMAX
               IF ( KVAR(i).NE.0 ) THEN
                  rl = LOG10(ELMU(i)/ABS(ELM(i)))
                  IF ( rl.GE.Rlr ) ELML(i) = 1.
                  WRITE (22,99005) i , rl , ELMU(i)/lngt
99005             FORMAT (6X,1I3,18X,1F4.1,20X,1E7.2)
               ENDIF
            ENDDO
            WRITE (22,99006)
99006       FORMAT (2X////40X,'ANALYSIS OF SIGNIFICANT DEPENDENCES'//)
            DO i = 1 , MEMAX
               IF ( KVAR(i).NE.0 ) THEN
                  lc = 0
                  IF ( ELML(i).GE..5 ) THEN
                     REWIND 17
 552                 READ (17,*) indx , iexp , sgm , ni , nf , al
                     IF ( indx.EQ.0 ) THEN
                        np = 0
                        nm = 0
                        DO j = 1 , lc
                           u = CORF(j,1)*CORF(j,2)*2.
                           IF ( ABS(u)/ELMU(i).GE..05 ) THEN
                              IF ( u.LT.0. ) nm = nm + 1
                              IF ( u.GT.0. ) np = np + 1
                           ENDIF
                        ENDDO
                        WRITE (22,99007) i , np , nm
99007                   FORMAT (1X/5X,10('*'),5X,'M.E.',1X,1I3,5X,1I3,
     &                          1X,'POSITIVE COMPONENTS',20X,1I3,1X,
     &                          'NEGATIVE COMPONENTS'///30X,'POSITIVE',
     &                          52X,'NEGATIVE'//5X,'EXP',2X,
     &                          'TRANSITION',2X,'SIGMA',3X,'DERIVATIVE',
     &                          3X,'D(SIGMA**2)/D(ME)',4X,'I',1X,'EXP',
     &                          2X,'TRANSITION',2X,'SIGMA',3X,
     &                          'DERIVATIVE',3X,'D(SIGMA**2)/D(ME)')
                        DO l = 1 , K
                           ump = 0.
                           umm = 0.
                           DO j = 1 , lc
                              u = 2.*CORF(j,1)*CORF(j,2)
                              IF ( u.LT.0. ) THEN
                                 IF ( u.LE.umm ) THEN
                                    umm = u
                                    jm = j
                                 ENDIF
                              ELSEIF ( u.GE.ump ) THEN
                                 ump = u
                                 jp = j
                              ENDIF
                           ENDDO
                           WRITE (22,99008) IY(jp,1) , IY(jp,2) , 
     &                            IY(jp,3) , CORF(jp,1) , CORF(jp,2) , 
     &                            ump , IY(jm,1) , IY(jm,2) , IY(jm,3) , 
     &                            CORF(jm,1) , CORF(jm,2) , umm
99008                      FORMAT (6X,1I2,3X,1I2,'--',1I2,5X,1F4.1,4X,
     &                             1E9.2,7X,1E9.2,9X,'I',2X,1I2,3X,1I2,
     &                             '--',1I2,5X,1F4.1,4X,1E9.2,7X,1E9.2)
                           CORF(jp,1) = 0.
                           CORF(jm,1) = 0.
                        ENDDO
                     ELSE
                        IF ( indx.EQ.i ) THEN
                           lc = lc + 1
                           IY(lc,1) = iexp
                           IY(lc,2) = ni
                           IY(lc,3) = nf
                           CORF(lc,1) = sgm
                           CORF(lc,2) = al
                        ENDIF
                        GOTO 552
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
            RETURN
         ELSE
            ll = ll + 1
            ch = ch + sgm*sgm
            ux = 2.*sgm*u
            ELM(indx) = ELM(indx) + ux
            ELMU(indx) = ELMU(indx) + ABS(ux)
            GOTO 550
         ENDIF
      ELSE
         loc = (iexp-1)*LP2 + indx
         sgm = (e-a*EP(iexp))/b
         u = dy*EP(iexp)*b + a*ZETA(loc)/b
         WRITE (17,*) indx , iexp , sgm , ni , nf , u
         GOTO 500
      ENDIF
      END
      SUBROUTINE MIXR(Nw,Ipsw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , dl , DMIX , DMIXE , ELM , ELML , ELMU , SA , 
     &       TAU
      INTEGER*4 i , IMIX , INTR , inx , inx1 , IPS1 , Ipsw , it , KSEQ , 
     &          LNY , NDL , Nw
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /MIXD  / DMIXE(20,2) , DMIX(20) , IMIX(20) , NDL
      COMMON /LOGY  / LNY , INTR , IPS1
      IF ( NDL.EQ.0 ) RETURN
      Nw = Nw + NDL
      DO i = 1 , NDL
         it = IMIX(i)
         inx = KSEQ(it,1)
         inx1 = KSEQ(it,2)
         IF ( ABS(ELM(inx1)).LT.1.E-5 ) ELM(inx1) = 1.E-5
         dl = DMIX(i)*ELM(inx)/ELM(inx1)
         IF ( Ipsw.EQ.1 ) DMIX(i) = dl
         Chi = Chi + (dl-DMIXE(i,1))**2/DMIXE(i,2)/DMIXE(i,2)
         IF ( LNY.EQ.1 ) Chilo = Chilo + 
     &                           (DMIXE(i,1)*LOG(ABS(dl/DMIXE(i,1)))
     &                           /DMIXE(i,2))**2
      ENDDO
      IF ( Ipsw.EQ.0 ) RETURN
      WRITE (22,99001)
99001 FORMAT (1X//10X,'E2/M1 MIXING RATIOS'/10X,'TRANSITION',10X,
     &        'EXP.DELTA',10X,'CALC.DELTA',10X,'SIGMA'/)
      DO i = 1 , NDL
         dl = (DMIX(i)-DMIXE(i,1))/DMIXE(i,2)
         it = IMIX(i)
         WRITE (22,99002) KSEQ(it,3) , KSEQ(it,4) , DMIXE(i,1) , DMIX(i)
     &                    , dl
99002    FORMAT (10X,1I2,'---',1I2,14X,1F7.2,12X,1F7.2,13X,1F5.2)
      ENDDO
      RETURN
      END
      SUBROUTINE COORD(Wth,Wph,Wthh,Naa,Ifw,Pfi,Wpi,Wtlb,Lz,Tyy,Tzz)
      IMPLICIT NONE
      REAL*8 DS , DSE , DSG , EP , EPS , EROOT , FIEX , ga , gi , Pfi , 
     &       rade , rmass , TACOS , TASIN , thetb , TLBDG , ttcm , Tyy , 
     &       Tzz , VINF
      REAL*8 wpa , Wph , Wpi , ws , Wth , Wthh , Wtlb , XA , XA1 , xaa , 
     &       xph , xth , xthh , XV , YV , za , za1 , zb , zl , ZV
      INTEGER*4 i , IAXS , IEXP , Ifw , ISKIN , IZ , IZ1 , Lz , Naa , 
     &          NEXPT
      DIMENSION Pfi(101) , Wpi(11,2)
      COMMON /VLIN  / XV(51) , YV(51) , ZV(20) , DSG(20) , DSE(20) , DS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /SECK  / ISKIN(50)
      DATA rade/57.2957795/
      IF ( Ifw.EQ.0 ) THEN
         Tyy = Wth - Wthh
         Tzz = Wth + Wthh
      ENDIF
      xth = Wth/rade
      xph = Wph/rade
      xthh = Wthh/rade
      zl = TAN(xthh)
      za = COS(xth)
      za1 = SIN(xth)
      zb = COS(xthh)
      rmass = XA1(Lz)/XA
      IF ( IZ1(Lz).LT.0 ) rmass = 1./rmass
      IF ( Ifw.NE.2 ) THEN
         ws = (Tzz-Tyy)/(Naa+1)
         IF ( Ifw.EQ.1 ) ws = (Tzz-Tyy)/(Naa-1)
      ENDIF
      DO i = 1 , Naa
         IF ( Ifw.NE.2 ) THEN
            IF ( Ifw.EQ.0 ) YV(i) = Tyy + i*ws
            xaa = (Tyy+ws*(i-1))/rade
            IF ( Ifw.EQ.1 .AND. (i.EQ.1 .OR. i.EQ.Naa) ) THEN
               Pfi(i) = 0.
               GOTO 100
            ELSE
               IF ( Ifw.EQ.0 ) xaa = YV(i)/rade
            ENDIF
         ELSE
            xaa = ABS(Wtlb)/rade
            IF ( Wtlb.GT.0. ) GOTO 50
            IF ( IZ1(Lz).LT.0 ) THEN
               IF ( XA.LE.XA1(Lz) ) GOTO 20
            ELSEIF ( XA1(Lz).LE.XA ) THEN
               GOTO 20
            ENDIF
            IF ( ISKIN(Lz).EQ.0 ) THEN
               ttcm = xaa - TASIN(rmass*SIN(xaa))
               xaa = ABS(ttcm)/2.
               GOTO 50
            ENDIF
 20         ttcm = xaa + TASIN(rmass*SIN(xaa))
            xaa = (3.14159265-ttcm)/2.
         ENDIF
 50      gi = (za-COS(xaa)/zb)/(zl*za1)
         ga = TACOS(gi)
         wpa = ATAN(zl*SIN(ga)/(za1+zl*COS(ga)*za))
         wpa = ABS(wpa)
         IF ( Ifw.EQ.2 ) THEN
            FIEX(Lz,1) = (xph-wpa)
            FIEX(Lz,2) = (xph+wpa)
         ELSEIF ( Ifw.EQ.1 ) THEN
            Pfi(i) = 2.*wpa*rade
         ELSE
            Wpi(i,1) = (xph-wpa)*rade
            Wpi(i,2) = (xph+wpa)*rade
         ENDIF
 100  ENDDO
      IF ( Wtlb.LT.0. .AND. Ifw.EQ.0 ) THEN
         DO i = 1 , Naa
            xaa = YV(i)/rade
            thetb = ATAN(SIN(2.*xaa)/(rmass-COS(2.*xaa)))*rade
            IF ( thetb.LT.0. ) thetb = 180. + thetb
            YV(i) = -1.*thetb
            Wpi(i,1) = Wpi(i,1) + 180.
            Wpi(i,2) = Wpi(i,2) + 180.
         ENDDO
      ENDIF
      END
      SUBROUTINE CHMEM(Nw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , di , EAMX , ELM , ELML , ELMU , SA
      INTEGER*4 ia , IAMX , IAMY , ib , NAMX , Nw
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      IF ( NAMX.EQ.0 ) RETURN
      Nw = Nw + NAMX
      DO ia = 1 , NAMX
         ib = IAMX(ia)
         IF ( IAMY(ia,1).NE.IAMY(ia,2) ) THEN
            di = (ELM(ib)-EAMX(ia,1))/EAMX(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &              /EAMX(ia,2))**2
            Chi = Chi + di*di
         ELSE
            di = (ELM(ib)-EAMX(ia,1))/EAMX(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &              /EAMX(ia,2))**2
            Chi = Chi + di*di
         ENDIF
      ENDDO
      END
      SUBROUTINE PTICC(Idr)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , cone1 , cone2 , conm1 , CONV , DIPOL , EN , 
     &       enet , SPIN , TAU , ZPOL
      INTEGER*4 Idr , iinx , ISO , KSEQ , l , LAMDA , LAMMAX , LDNUM , 
     &          LEAD , MULTI , nf , ni
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      WRITE (22,99001)
99001 FORMAT (1X//20X,'CALCULATED INTERNAL CONVERSION ',
     &        'COEFFICIENTS FOR E1,E2 AND M1'//5X,'NI',5X,'NF',7X,'II',
     &        8X,'IF',9X,'ENERGY(MEV)',6X,'ICC(E1)',8X,'ICC(E2)',8X,
     &        'ICC(M1)')
      DO l = 1 , Idr
         iinx = KSEQ(l,1)
         ni = KSEQ(l,3)
         nf = KSEQ(l,4)
         enet = EN(ni) - EN(nf)
         cone2 = CONV(enet,2)
         IF ( ABS(SPIN(ni)-SPIN(nf)).GT.2. ) cone2 = 0.
         conm1 = 0.
         cone1 = 0.
         IF ( iinx.LE.MULTI(1) ) cone1 = CONV(enet,1)
         IF ( ABS(SPIN(ni)-SPIN(nf)).LT.2. ) conm1 = CONV(enet,4)
         WRITE (22,99002) ni , nf , SPIN(ni) , SPIN(nf) , enet , cone1 , 
     &                    cone2 , conm1
99002    FORMAT (5X,I2,5X,I2,7X,F4.1,6X,F4.1,9X,F6.4,8X,E9.4,6X,E9.4,6X,
     &           E9.4)
      ENDDO
      END
      FUNCTION RNDM(Se)
      IMPLICIT NONE
      REAL*8 ai , p , r , RNDM , rxdm , Se , t , u
      INTEGER*4 i
      IF ( Se.GT.32000. ) Se = 100.*t + .511
      Se = Se*Se
      u = LOG10(Se)
      i = INT(u) + 1
      t = Se/(10.**i)
      r = SQRT(SQRT(SQRT(t)))
      p = SQRT(SQRT(SQRT(.1)))
      rxdm = (r-p)/(1.-p)
      rxdm = 10.*rxdm
      ai = DBLE(INT(rxdm))
      RNDM = rxdm - ai
      RETURN
      END
      SUBROUTINE KONTUR(Idr,Chis0,Chil,Ifbf,Inpo,Jj,Sh,Bten,Rem)
      IMPLICIT NONE
      REAL*8 ac , Bten , c , Chil , chilo , Chis0 , chis1 , chis2 , d1 , 
     &       d2 , DEVD , DEVU , DS , DSE , DSG , ELM , ELML , ELMU , f , 
     &       h
      REAL*8 HLM , Rem , RK4 , SA , sajj , Sh , t , v , ww , x , XV , 
     &       y , YV , ZV
      INTEGER*4 i , Idr , Ifbf , Inpo , INTR , IPS1 , itl , IVAR , ix , 
     &          j , Jj , l , LMAXE , LNY , m , MAGEXC , MEMAX , MEMX6 , 
     &          NWR
      DIMENSION f(3) , Bten(1200)
      COMMON /VLIN  / XV(51) , YV(51) , ZV(20) , DSG(20) , DSE(20) , DS
      COMMON /DFTB  / DEVD(500) , DEVU(500)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /HHH   / HLM(500)
      COMMON /ILEWY / NWR
      COMMON /LOGY  / LNY , INTR , IPS1
      LNY = 0
      h = .05*ABS(HLM(Jj))
      IF ( Inpo.NE.-1 ) h = ABS(Sh)
 100  INTR = 0
      sajj = ABS(SA(Jj))
      DO l = 1 , MEMAX
         ELM(l) = HLM(l)
         SA(l) = SA(l)/sajj
      ENDDO
      YV(1) = 0.
      XV(1) = HLM(Jj)
      f(3) = 1.
      i = 1
 200  itl = 0
      v = ELMU(Jj) - ELM(Jj)
      IF ( SA(Jj).LT.0. ) v = ELM(Jj) - ELML(Jj)
      IF ( h.GT.v ) itl = 1
      IF ( h.GT.v ) h = v
      i = i + 1
      f(1) = f(3)
      DO j = 1 , MEMAX
         ELM(j) = .5*h*SA(j) + ELM(j)
      ENDDO
      CALL LIMITS
      CALL FTBM(3,chis1,Idr,1,chilo,Bten)
      IF ( chis1.LE.Chis0 ) THEN
         IF ( Inpo.EQ.-1 ) WRITE (22,99003) Jj , ELM(Jj) , chis1
         IF ( chis1.LE.Chil .AND. Inpo.NE.-1 ) THEN
            Ifbf = 1
            ix = 1
            Chil = chis1
            WRITE (22,99004) Chil
            GOTO 500
         ENDIF
      ENDIF
 300  ww = .5*(Chis0-chis1)*NWR
      IF ( ww.GE.Rem ) GOTO 700
      f(2) = EXP(ww)
      IF ( i.EQ.2 .AND. f(2).LT..1 .AND. ABS(XV(1)-HLM(Jj)).LT.1E-9 )
     &     THEN
         h = h/2.
         GOTO 100
      ELSE
         DO j = 1 , MEMAX
            ELM(j) = ELM(j) + .5*SA(j)*h
         ENDDO
         v = ELM(Jj)
         CALL LIMITS
         IF ( ABS(v-ELM(Jj)).GT.1.E-6 ) itl = 1
         CALL FTBM(3,chis2,Idr,1,chilo,Bten)
         IF ( chis2.LE.Chis0 ) THEN
            IF ( Inpo.EQ.-1 ) WRITE (22,99003) Jj , ELM(Jj) , chis2
            IF ( chis2.LE.Chil .AND. Inpo.NE.-1 ) THEN
               Ifbf = 1
               ix = 2
               Chil = chis2
               WRITE (22,99004) Chil
               GOTO 500
            ENDIF
         ENDIF
      ENDIF
 400  ww = .5*(Chis0-chis2)*NWR
      IF ( ww.GT.Rem ) GOTO 700
      f(3) = EXP(ww)
      IF ( itl.EQ.1 ) WRITE (22,99001) Jj
99001 FORMAT (5X,'WARNING-ME(',1I3,')',5X,
     &        'INTEGRATION STOPPED AT THE LIMIT')
      IF ( i.EQ.2 ) THEN
         IF ( itl.NE.1 ) THEN
            IF ( f(3).LT..1 .AND. ABS(XV(1)-HLM(Jj)).LT.1.E-9 ) THEN
               h = h/2.
               GOTO 100
            ELSEIF ( f(1).LE.f(2) .OR. f(2).LE.f(3) ) THEN
               IF ( f(1).LT.f(2) .AND. f(2).GT.f(3) ) THEN
                  d1 = f(2) - f(1)
                  d2 = f(3) - f(1)
                  ac = (d2-4.*d1)*h/(d2-2.*d1)/4.
                  DO l = 1 , MEMAX
                     ELM(l) = (ELM(l)-h*SA(l)) + ac*SA(l)
                  ENDDO
                  CALL LIMITS
                  XV(1) = ELM(Jj)
                  i = 1
                  CALL FTBM(3,chis1,Idr,1,chilo,Bten)
                  ww = .5*(Chis0-chis1)*NWR
                  IF ( ww.GE.Rem ) GOTO 700
                  f(3) = EXP(ww)
                  GOTO 200
               ELSE
                  i = 1
                  XV(1) = ELM(Jj)
                  IF ( Inpo.EQ.-1 ) h = 2.*h
                  GOTO 200
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      y = YV(i-1)
      YV(i) = RK4(y,h,f)
      XV(i) = ELM(Jj)
      IF ( NWR*(chis2-Chis0).LT.2. .AND. Inpo.EQ.-1 ) h = 2.*h
      IF ( itl.EQ.1 ) GOTO 600
      IF ( f(3).GE.1.E-3 ) GOTO 200
      GOTO 600
 500  REWIND 17
      DO l = 1 , MEMAX
         WRITE (17,*) ELM(l)
      ENDDO
      IF ( ix.EQ.1 ) GOTO 300
      IF ( ix.NE.2 ) GOTO 200
      GOTO 400
 600  c = YV(i)
      m = 0
      DO l = 1 , i
         YV(l) = 1.00001 - YV(l)/c
         IF ( m.EQ.0 .AND. YV(l).LT..317 ) m = l
      ENDDO
      x = (XV(m)-XV(m-1))*(.317-YV(m))/(YV(m-1)-YV(m))
      t = XV(m) - x - HLM(Jj)
      IF ( t.GE.0. ) DEVU(Jj) = t
      IF ( t.LT.0. ) DEVD(Jj) = t
      RETURN
 700  WRITE (22,99002) Jj
99002 FORMAT (5X,'** WARNING **',/,2X,'ME=',1I3,2X,
     &     'TOO FAR FROM THE MINIMUM TO CARRY OUT THE ERROR ESTIMATION!'
     &     ,/)
99003 FORMAT (5X,'ELM(',1I3,')=',1F10.6,5X,'CHISQ=',1E12.4)
99004 FORMAT (10X,'BETTER POINT FOUND...MATRIX ELEMENTS WRITTEN ON 17',
     &        3X,'CHISQ=',1E12.4)
      END
      FUNCTION RK4(Y,H,F)
      IMPLICIT NONE
      REAL*8 F , H , RK4 , Y
      DIMENSION F(3)
      RK4 = Y + H*(F(1)+4.*F(2)+F(3))/6.
      END
      SUBROUTINE QFIT(Qui,Tau1,Tau2,Eng,Xl1,Cf,Nl,Ind)
      IMPLICIT NONE
      REAL*8 ca , cb , Cf , cm , cn , co , d , d1 , d2 , Eng , Qui , 
     &       Tau1 , Tau2 , Xl1
      INTEGER*4 Ind , ind1 , k , Nl
      DIMENSION Tau1(10) , Eng(10) , Tau2(10,7) , Xl1(7) , Qui(8,10) , 
     &          Cf(8,2)
      CALL GAMATT(Qui,Tau1,Tau2,Xl1,Nl)
      ind1 = 5
      IF ( Ind.EQ.4 ) ind1 = 6
      IF ( Ind.EQ.5 ) ind1 = 7
      DO k = 1 , 8
         co = Qui(k,Ind)
         cn = Qui(k,10)
         cm = Qui(k,ind1)
         ca = (Eng(ind1)-Eng(Ind))**2
         cb = (Eng(10)-Eng(Ind))**2
         d = ca*(co-cn) - cb*(co-cm)
         d1 = ca*cm*(co-cn) - cb*cn*(co-cm)
         d2 = ca*cb*(cn-cm)
         Cf(k,1) = d1/d
         Cf(k,2) = d2/d
      ENDDO
      END
      SUBROUTINE GAMATT(Qui,Tau1,Tau2,Xl1,Nl)
      IMPLICIT NONE
      INTEGER*4 i , i1 , k , Nl
      REAL*8 q , Qui , tau , Tau1 , Tau2 , thing , thing1 , thing3 , Xl1
      DIMENSION Tau1(10) , Tau2(10,7) , Xl1(7) , thing3(10) , q(9) , 
     &          Qui(8,10)
      DO i = 1 , 10
         i1 = 1
         thing3(i) = 0.
 50      thing1 = -Tau2(i,i1)*Xl1(i1) + thing3(i)
         i1 = i1 + 1
         thing3(i) = thing1
         IF ( i1.LE.Nl ) GOTO 50
      ENDDO
      DO i = 1 , 10
         tau = Tau1(i)
         thing = thing3(i)
         CALL GCF(tau,thing,q)
         DO k = 2 , 9
            Qui(k-1,i) = q(k)
         ENDDO
      ENDDO
      END
      SUBROUTINE GCF(Tau,Thing,Q)
      IMPLICIT NONE
      REAL*8 A , b , D , dl , ev , ex , f , fint , od , ODL , Q , R , 
     &       Tau , Thing , XL , xm , yl , yu
      INTEGER*4 i , j , k , m
      COMMON /DIMX  / A , R , XL , D , ODL(200)
      DIMENSION f(101) , b(4) , Q(9)
      b(1) = ATAN2(A,D+XL)
      b(2) = ATAN2(A,D)
      b(3) = ATAN2(R,D+XL)
      b(4) = ATAN2(R,D)
      DO k = 1 , 9
         Q(k) = 0.0
         DO j = 1 , 3
            yl = b(j)
            yu = b(j+1)
            dl = (yu-yl)/100.
            DO m = 1 , 101
               xm = yl + dl*(m-1)
               IF ( j.EQ.2 ) THEN
                  ex = -Tau*XL/COS(xm)
               ELSEIF ( j.EQ.3 ) THEN
                  ex = Tau*(D*TAN(xm)-R)/SIN(xm)
               ELSE
                  ex = Tau*(A-(D+XL)*TAN(xm))/SIN(xm)
               ENDIF
               f(m) = SIN(xm)*(1-EXP(ex))*EXP(Thing/COS(xm))
               IF ( j.EQ.1 ) f(m) = f(m)*EXP(-Tau*(A/SIN(xm)-D/COS(xm)))
               IF ( k.EQ.1 ) THEN
               ELSEIF ( k.EQ.3 ) THEN
                  f(m) = f(m)*(1.5*COS(xm)**2-0.5)
               ELSEIF ( k.EQ.4 ) THEN
                  f(m) = f(m)*(2.5*COS(xm)**3-1.5*COS(xm))
               ELSEIF ( k.EQ.5 ) THEN
                  f(m) = f(m)*(4.375*COS(xm)**4-3.75*COS(xm)**2+.375)
               ELSEIF ( k.EQ.6 ) THEN
                  f(m) = f(m)*((63.*COS(xm)**5-70.*COS(xm)**3+15.)/8.)
               ELSEIF ( k.EQ.7 ) THEN
                  f(m) = f(m)
     &                   *((21.*COS(xm)**2*(11.*COS(xm)**4-15.*COS(xm)
     &                   **2+5.)-5.)/16.)
               ELSEIF ( k.EQ.8 ) THEN
                  f(m) = f(m)
     &                   *(429.*COS(xm)**7-693.*COS(xm)**5+315.*COS(xm)
     &                   **3-35.*COS(xm))/16.
               ELSEIF ( k.EQ.9 ) THEN
                  f(m) = f(m)
     &                   *(6435.*COS(xm)**8-12012.*COS(xm)**6+6930.*COS
     &                   (xm)**4-1260.*COS(xm)**2+35.)/128.
               ELSE
                  f(m) = f(m)*COS(xm)
               ENDIF
            ENDDO
            ev = 0.0
            od = 0.0
            DO m = 2 , 98 , 2
               ev = ev + f(m)
               od = od + f(m+1)
            ENDDO
            fint = dl/3.*(f(1)+4.*(ev+f(100))+2.*od+f(101))
            Q(k) = Q(k) + fint
         ENDDO
      ENDDO
      DO i = 1 , 8
         Q(i+1) = Q(i+1)/Q(1)
      ENDDO
      Q(1) = Q(1)/2.
      END
      COMPLEX*16 FUNCTION TCEXP(Z)
      IMPLICIT NONE
      REAL*8 a , b , c , d
      COMPLEX*16 Z
      a = DBLE(Z)
      b = IMAG(Z)
      a = EXP(a)
      c = a*COS(b)
      d = a*SIN(b)
      TCEXP = CMPLX(c,d)
      END
      FUNCTION TCABS(Z)
      IMPLICIT NONE
      REAL*8 a , b , TCABS
      COMPLEX*16 Z
      a = DBLE(Z)
      b = IMAG(Z)
      IF ( ABS(a).LT.1.E-16 ) a = 0.
      IF ( ABS(b).LT.1.E-16 ) b = 0.
      TCABS = SQRT(a*a+b*b)
      END
      FUNCTION TASIN(X)
      IMPLICIT NONE
      REAL*8 dol , TASIN , test , war , X
      test = ABS(X) - 1.
      IF ( ABS(test).LT.1.E-9 ) THEN
         TASIN = 1.570796327
         IF ( X.LT.0. ) TASIN = -1.570796327
         GOTO 99999
      ENDIF
      dol = SQRT(1.-X*X)
      war = X/dol
      TASIN = ATAN(war)
      RETURN
99999 END
      FUNCTION TACOS(X)
      IMPLICIT NONE
      REAL*8 TACOS , TASIN , X
      TACOS = 1.570796327 - TASIN(X)
      END
      SUBROUTINE OPENF
      IMPLICIT NONE
      INTEGER*4 i , IBPS , j , JZB , k
      CHARACTER name*60 , opt1*20 , opt2*20
      COMMON /SWITCH/ JZB , IBPS
 100  READ (JZB,*) i , j , k
      IF ( i.EQ.0 ) RETURN
      IF ( j.EQ.1 ) opt1 = 'OLD'
      IF ( j.EQ.2 ) opt1 = 'NEW'
      IF ( j.EQ.3 ) opt1 = 'UNKNOWN'
      IF ( k.EQ.1 ) opt2 = 'FORMATTED'
      IF ( k.EQ.2 ) opt2 = 'UNFORMATTED'
      READ (JZB,99001) name
99001 FORMAT (A)
      IF ( i.NE.25 .AND. i.NE.26 )
     &     OPEN (i,IOSTAT=k,FILE=name,STATUS=opt1,FORM=opt2)
      IF ( i.EQ.25 .OR. i.EQ.26 ) k = 0
c      IF (K.EQ.0) WRITE(6,1030) 'OPENED ',NAME
c 1030 FORMAT (1X,2A)
c      WRITE(6,1010) ' IO-num = ',I,OPT1,OPT2
c 1010 FORMAT (1X,A,I4,2(1x,A))
      IF ( k.EQ.0 ) GOTO 100
c      WRITE (6,1020) 'PROBLEMS OPENING ',NAME,K
c 1020 FORMAT(A,A,I6)
      END
      SUBROUTINE EFFIX(Ipd,En,Effi)
      IMPLICIT NONE
      REAL*8 ABC , AKAVKA , d , Effi , En , enek , enl , pw , s , t ,
     &     THICK , w , xx , yy
      INTEGER*4 i , Ipd , ixi , j , l , ll , n
      DIMENSION xx(51) , yy(51)
      COMMON /EFCAL / ABC(8,10) , AKAVKA(8,200) , THICK(200,7)
      Effi = 1.E-6
      En = En + 1.E-24
      enl = LOG(En)
      DO i = 1 , 10
         ll = 11 - i
         j = ll
         IF ( enl.GE.ABC(8,ll) ) GOTO 100
         j = -1
      ENDDO
 100  IF ( j.EQ.-1 ) Effi = 1.E-10
      IF ( j.EQ.-1 ) RETURN
      IF ( j.EQ.1 .OR. j.EQ.10 ) THEN
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THICK(Ipd,l)).GE.1.E-9 ) THEN
               t = EXP(ABC(l,j))
               d = THICK(Ipd,l)
               s = s + t*d
            ENDIF
         ENDDO
      ELSE
         IF ( j.EQ.9 ) THEN
            xx(1) = ABC(8,8)
            xx(2) = ABC(8,9)
            xx(3) = ABC(8,10)
         ELSE
            xx(1) = ABC(8,j)
            xx(2) = ABC(8,j+1)
            xx(3) = ABC(8,j+2)
         ENDIF
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THICK(Ipd,l)).GE.1.E-9 ) THEN
               IF ( j.EQ.9 ) THEN
                  yy(1) = ABC(l,8)
                  yy(2) = ABC(l,9)
                  yy(3) = ABC(l,10)
               ELSE
                  yy(1) = ABC(l,j)
                  yy(2) = ABC(l,j+1)
                  yy(3) = ABC(l,j+2)
               ENDIF
               CALL LAGRAN(xx,yy,3,0,enl,t,1,1)
               s = s + EXP(t)*THICK(Ipd,l)
            ENDIF
         ENDDO
      ENDIF
      Effi = EXP(-s)
c FITEFF or GREMLIN check
      IF ( AKAVKA(8,Ipd).LE.-999. ) THEN
C     LEUVEN CALIBRATION
         Effi = AKAVKA(1,Ipd)
         enek = 1000.*En
         w = LOG(enek)
         DO ixi = 1 , 6
            Effi = Effi + AKAVKA(ixi+1,Ipd)*w**ixi
         ENDDO
         Effi = EXP(Effi)
         GOTO 99999
      ELSEIF ( AKAVKA(5,Ipd).GT.0. .AND. AKAVKA(5,Ipd).LT.10. ) THEN
c FITEFF eff. calib. by P.Olbratowski use
c PJN@2000
         w = LOG(En/AKAVKA(5,Ipd))
         pw = AKAVKA(2,Ipd)*w
         IF ( En.LT.AKAVKA(5,Ipd) ) pw = pw + 
     &        w*w*(AKAVKA(3,Ipd)+w*AKAVKA(4,Ipd))
         Effi = Effi*EXP(pw)*AKAVKA(1,Ipd)
         RETURN
      ELSEIF ( AKAVKA(5,Ipd).LT.10. ) THEN
c GREMLIN
         w = LOG(20.*En)
         pw = AKAVKA(1,Ipd) + AKAVKA(2,Ipd)*w + AKAVKA(3,Ipd)
     &        *w*w + AKAVKA(4,Ipd)*w*w*w
         Effi = Effi*EXP(pw)
         IF ( ABS(AKAVKA(5,Ipd)).GE.1.E-9 ) THEN
            n = INT(AKAVKA(6,Ipd)+.1)
            pw = w**n
            w = AKAVKA(5,Ipd)/pw
            Effi = Effi*EXP(w)
         ENDIF
         IF ( ABS(AKAVKA(8,Ipd)).LT.1.E-9 ) RETURN
         w = (AKAVKA(7,Ipd)-1000.*En)/AKAVKA(8,Ipd)
         pw = EXP(w)
         IF ( ABS(pw-1.).LT.1.E-6 ) WRITE (22,99001)
99001    FORMAT (5x,'***** CRASH - EFFIX *****')
         Effi = Effi/(1.-pw)
         RETURN
      ENDIF
c     JAERI calibration - TC, Nov.2000
      w = LOG(En/.511)
      Effi = EXP(AKAVKA(1,Ipd)+AKAVKA(2,Ipd)
     &       *w-EXP(AKAVKA(3,Ipd)+AKAVKA(4,Ipd)*w))
      RETURN
99999 END
      SUBROUTINE ADHOC(Oph,Idr,Nfd,Ntap,Iyr)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AGELI , BRAT , CC , CORF , DELTA , DIPOL , 
     &       DIX , DMIX , DMIXE , DYEX , EAMX , EG , EN , ENDEC , ENZ , 
     &       EP , ODL , Q
      REAL*8 SPIN , TAU , TIMEL , TLBDG , UPL , VINF , wamx , wbra , 
     &       wdl , wlf , XA , XA1 , YEXP , YGN , YGP , YNRM , ZPOL
      INTEGER*4 IAMX , IAMY , iax , IBPS , IBRC , Idr , IDRN , iexp1 , 
     &          IFMO , ILE , ilft , IMIX , iosr , ipri , IPRM , ISO , 
     &          isrt1 , ITMA , ITS , iuf
      INTEGER*4 IVAR , IY , Iyr , IZ , IZ1 , jic , jicc , juf , JZB , 
     &          KSEQ , lb , li , licc , LIFCT , llia , LMAXE , lxt , 
     &          MAGEXC , MEM , MEMAX
      INTEGER*4 MEMX6 , n1 , n2 , NAMX , NANG , NBRA , ndas , NDL , 
     &          NDST , ndtp , NEXPT , Nfd , NICC , nistr , NLIFT , ns1 , 
     &          ns2 , ns3 , ns4 , Ntap
      INTEGER*4 nvare , NYLDE
      CHARACTER*4 Oph
      COMMON /SWITCH/ JZB , IBPS
      COMMON /CCCDS / NDST(50)
      COMMON /DIMX  / DIX(4) , ODL(200)
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /LIFE  / NLIFT
      COMMON /MIXD  / DMIXE(20,2) , DMIX(20) , IMIX(20) , NDL
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      COMMON /LIFE1 / LIFCT(50) , TIMEL(2,50)
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /YTEOR / YGN(500) , YGP(500) , IFMO
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /PRT   / IPRM(20)
      COMMON /TRB   / ITS
      iosr = 0
      READ (JZB,*) IFMO
      READ (JZB,*) NICC , nistr
      READ (JZB,*) (EG(jicc),jicc=1,NICC)
      Iyr = 1
      DO jic = 1 , nistr
         READ (JZB,*) isrt1
         IF ( isrt1.GT.6 ) isrt1 = isrt1 - 3
         READ (JZB,*) (CC(jicc,isrt1),jicc=1,NICC)
      ENDDO
      READ (JZB,*) (NANG(jicc),jicc=1,NEXPT)
      REWIND 9
      READ (9,*) Nfd
      DO jicc = 1 , Nfd
         READ (9,*) ODL(jicc)
         READ (9,*) ENZ(jicc)
         DO isrt1 = 1 , 8
            READ (9,*) (Q(licc,jicc,isrt1),licc=1,3)
         ENDDO
      ENDDO
      DO jic = 1 , NEXPT
         juf = NANG(jic)
         IF ( juf.LT.0 ) THEN
            juf = ABS(juf)
            DO jicc = 1 , juf
               AGELI(jic,jicc,1) = AGELI(jic-1,jicc,1)
               AGELI(jic,jicc,2) = AGELI(jic-1,jicc,2)
               ITMA(jic,jicc) = ITMA(jic-1,jicc)
            ENDDO
            IF ( Oph.NE.'GOSI' ) NANG(jic) = ABS(NANG(jic))
         ELSE
            READ (JZB,*) (ITMA(jic,jicc),jicc=1,juf)
            READ (JZB,*) (AGELI(jic,jicc,1),jicc=1,juf)
            READ (JZB,*) (AGELI(jic,jicc,2),jicc=1,juf)
         ENDIF
      ENDDO
      CALL SEQ(Idr)
      DO jic = 1 , NEXPT
         juf = NANG(jic)
         juf = ABS(juf)
         DO jicc = 1 , juf
            DO lxt = 1 , 2
               AGELI(jic,jicc,lxt) = AGELI(jic,jicc,lxt)*.0174532925
            ENDDO
         ENDDO
      ENDDO
      TAU(1) = 1.E+25
      READ (JZB,*) ns1 , ns2
      DO li = 1 , Idr
         IF ( KSEQ(li,3).EQ.ns1 .AND. KSEQ(li,4).EQ.ns2 ) GOTO 100
      ENDDO
 100  IDRN = li
      IF ( Oph.NE.'GOSI' ) RETURN
      DO li = 1 , NEXPT
         juf = NANG(li)
         IF ( juf.LT.0 ) THEN
            juf = ABS(juf)
            NANG(li) = juf
            NDST(li) = NDST(li-1)
            DO jicc = 1 , juf
               UPL(jicc,li) = UPL(jicc,li-1)
               YNRM(jicc,li) = YNRM(jicc,li-1)
            ENDDO
         ELSE
            READ (JZB,*) NDST(li)
            ndas = NDST(li)
            READ (JZB,*) (UPL(jicc,li),jicc=1,ndas)
            READ (JZB,*) (YNRM(jicc,li),jicc=1,ndas)
         ENDIF
      ENDDO
      READ (JZB,*) Ntap
      IF ( Ntap.NE.0 ) THEN
         ipri = IPRM(2)
         CALL READY(Idr,Ntap,ipri)
         ndtp = 0
         DO iexp1 = 1 , NEXPT
            juf = NDST(iexp1)
            DO iuf = 1 , juf
               ndtp = ndtp + NYLDE(iexp1,iuf)
            ENDDO
         ENDDO
         nvare = 0
         DO iexp1 = 1 , MEMAX
            IF ( IVAR(iexp1).EQ.1 .OR. IVAR(iexp1).EQ.2 )
     &           nvare = nvare + 1
         ENDDO
         WRITE (22,99001) ndtp , nvare
99001    FORMAT (1X//5X,1I4,1X,'EXPERIMENTAL YIELDS',10X,1I3,1X,
     &           'MATRIX ELEMENTS TO BE VARIED'///)
      ENDIF
      READ (JZB,*) NBRA , wbra
      IF ( ITS.EQ.2 ) THEN
         REWIND 18
         WRITE (18,*) MEMAX
      ENDIF
      IF ( NBRA.NE.0 ) THEN
         WRITE (22,99002)
99002    FORMAT (40X,'BRANCHING RATIOS',//5X,'NS1',5X,'NF1',5X,'NS2',5X,
     &           'NF2',5X,'RATIO(1:2)',9X,'ERROR')
         DO lb = 1 , NBRA
            READ (JZB,*) ns1 , ns2 , ns3 , ns4 , BRAT(lb,1) , BRAT(lb,2)
            BRAT(lb,2) = BRAT(lb,2)/(SQRT(wbra)+1.E-10)
            WRITE (22,99003) ns1 , ns2 , ns3 , ns4 , BRAT(lb,1) , 
     &                       BRAT(lb,2)
99003       FORMAT (5X,1I2,6X,1I2,6X,1I2,6X,1I2,5X,1F10.5,5X,1F10.5)
            DO li = 1 , Idr
               IF ( KSEQ(li,3).EQ.ns3 .AND. KSEQ(li,4).EQ.ns4 ) THEN
                  IBRC(2,lb) = li
               ELSEIF ( KSEQ(li,3).EQ.ns1 .AND. KSEQ(li,4).EQ.ns2 ) THEN
                  IBRC(1,lb) = li
               ENDIF
            ENDDO
            IF ( ITS.EQ.2 ) THEN
               n1 = IBRC(1,lb)
               n2 = IBRC(2,lb)
               WRITE (18,*) KSEQ(n1,1) , KSEQ(n2,1)
               WRITE (18,*) KSEQ(n1,1) , KSEQ(n2,2)
               WRITE (18,*) KSEQ(n1,1) , KSEQ(n1,2)
               WRITE (18,*) KSEQ(n2,1) , KSEQ(n1,2)
               WRITE (18,*) KSEQ(n2,1) , KSEQ(n2,2)
               IF ( KSEQ(n1,2).NE.0 .AND. KSEQ(n2,2).NE.0 ) WRITE (18,*)
     &              KSEQ(n1,2) , KSEQ(n2,2)
            ENDIF
         ENDDO
         WRITE (22,99004) wbra
99004    FORMAT (5X,'BRANCHING RATIOS ARE TAKEN WITH WEIGHT',2X,1E14.6)
      ENDIF
      READ (JZB,*) NLIFT , wlf
      IF ( NLIFT.NE.0 ) THEN
         WRITE (22,99005)
99005    FORMAT (1X///30X,'LIFETIMES(PSEC)'///5X,'LEVEL',9X,'LIFETIME',
     &           5X,'ERROR'/)
         DO ilft = 1 , NLIFT
            READ (JZB,*) LIFCT(ilft) , TIMEL(1,ilft) , TIMEL(2,ilft)
            TIMEL(2,ilft) = TIMEL(2,ilft)/(SQRT(wlf)+1.E-10)
            WRITE (22,99006) LIFCT(ilft) , TIMEL(1,ilft) , TIMEL(2,ilft)
99006       FORMAT (7X,1I2,6X,1F10.2,3X,1F10.2)
         ENDDO
         WRITE (22,99007) wlf
99007    FORMAT (1X/10X,'LIFETIMES ARE TAKEN WITH WEIGHT',2X,1E14.6)
      ENDIF
      READ (JZB,*) NDL , wdl
      IF ( NDL.NE.0 ) THEN
         WRITE (22,99008)
99008    FORMAT (1X//20X,'EXPERIMENTAL E2/M1 MIXING RATIOS'///10X,
     &           'TRANSITION',12X,'DELTA',10X,'ERROR'/)
         DO li = 1 , NDL
            READ (JZB,*) ns1 , ns2 , DMIXE(li,1) , DMIXE(li,2)
            DMIXE(li,2) = DMIXE(li,2)/(SQRT(wdl)+1.E-10)
            WRITE (22,99012) ns1 , ns2 , DMIXE(li,1) , DMIXE(li,2)
            DO lb = 1 , Idr
               IF ( KSEQ(lb,3).EQ.ns1 .AND. KSEQ(lb,4).EQ.ns2 ) THEN
                  IMIX(li) = lb
                  DMIX(li) = .8326*(EN(ns1)-EN(ns2))
                  IF ( ITS.EQ.2 ) WRITE (18,*) KSEQ(lb,1) , KSEQ(lb,2)
               ENDIF
            ENDDO
         ENDDO
         WRITE (22,99009) wdl
99009    FORMAT (/10X,'E2/M1 MIXING RATIOS ARE TAKEN WITH WEIGHT',2X,
     &           1E14.6)
      ENDIF
      IF ( ITS.EQ.2 ) WRITE (18,*) iosr , iosr
      READ (JZB,*) NAMX , wamx
      IF ( NAMX.EQ.0 ) RETURN
      WRITE (22,99010)
99010 FORMAT (1X//30X,'EXPERIMENTAL MATRIX ELEMENT(S)'///10X,
     &        'TRANSITION',10X,'MAT.EL.',10X,'ERROR'/)
      DO iax = 1 , NAMX
         READ (JZB,*) llia , ns1 , ns2 , EAMX(iax,1) , EAMX(iax,2)
         IAMY(iax,1) = ns1
         IAMY(iax,2) = ns2
         EAMX(iax,2) = EAMX(iax,2)/(SQRT(wamx)+1.E-10)
         WRITE (22,99012) ns1 , ns2 , EAMX(iax,1) , EAMX(iax,2)
         IAMX(iax) = MEM(ns1,ns2,llia)
      ENDDO
      WRITE (22,99011) wamx
99011 FORMAT (/10X,' MATRIX ELEMENT(S) ARE TAKEN WITH WEIGHT',2X,1E14.6)
99012 FORMAT (10X,1I2,'---',1I2,14X,1F9.4,8X,1F9.4)
      END
      REAL*8 FUNCTION ELMT(Xi1,Xi2,Lam,Nb1,Nb2,Xk1,Xk2,Xm1,Xm2,Xm3)
      IMPLICIT NONE
      REAL*8 addt , fac , fct , pha1 , pha2 , s1 , s2 , WTHREJ ,
     &     Xi1 , Xi2 , Xk1 , Xk2 , xlam , Xm1 , Xm2 , Xm3 , xn
      INTEGER*4 i1 , i2 , ipha , k1 , k2 , l , la , Lam , llam , n , 
     &        Nb1 ,Nb2
      la = Lam
      IF ( la.GT.6 ) la = la - 6
      xlam = FLOAT(la)
      i1 = INT(2.*Xi1)
      i2 = INT(2.*Xi2)
      llam = 2*la
      k1 = INT(2.*Xk1)
      k2 = INT(2.*Xk2)
      fac = SQRT(2.*Xi1+1.)*SQRT(2.*Xi2+1.)
C-----In-band matrix element
      IF ( Nb1.NE.Nb2 ) THEN
C-----Interband, K-allowed
C-----One K=0
         IF ( ABS(k1-k2).GE.llam ) THEN
C-----Forbidden and K1-K2=lambda, Mikhailov formula
            addt = 0.
            IF ( k1.EQ.1 ) addt = (-1.)**((i1+1)/2)*(i1+1)/2.*Xm3
            xn = ABS(Xk1-Xk2) - xlam
            n = INT(xn+.1)
            IF ( n.EQ.0 ) THEN
               fct = 1.
            ELSEIF ( n.EQ.1 ) THEN
               fct = SQRT((Xi1-Xk1)*(Xi1+Xk1+1.))
            ELSE
               s1 = Xi1 - Xk1
               s2 = Xi1 + Xk1 + 1.
               DO l = 1 , n
                  s1 = s1*(Xi1-Xk1-FLOAT(l))
                  s2 = s2*(Xi1+Xk2+1.+FLOAT(l))
               ENDDO
               fct = SQRT(s1*s2)
            ENDIF
            pha1 = (-1.)**INT((Xi1-xlam+Xk2)+.1)
            ELMT = fac*pha1*fct*WTHREJ(i1,llam,i2,k2-llam,llam,-k2)
     &             *(Xm1+Xm2*(Xi2*(Xi2+1.)-Xi1*(Xi1+1.))+addt)
         ELSEIF ( k1.NE.0 .AND. k2.NE.0 ) THEN
C-----Both K's non-zero
            pha1 = (-1.)**((i1-llam+k2)/2)
            pha2 = (-1.)**((i1+k1)/2)*pha1
            ELMT = fac*(pha1*WTHREJ(i1,llam,i2,k1,k2-k1,-k2)
     &             *Xm1+pha2*WTHREJ(i1,llam,i2,-k1,k1+k2,-k2)*Xm2)
            RETURN
         ELSE
            ipha = (i1-llam+k2)/2
            IF ( k2.EQ.0 ) ipha = ((i2-llam+k1)/2)
            pha1 = (-1.)**ipha
            ELMT = fac*pha1*WTHREJ(i1,llam,i2,0,k2,-k2)*Xm1
            IF ( k2.EQ.0 ) ELMT = fac*pha1*WTHREJ(i2,llam,i1,0,k1,-k1)
     &                            *Xm1
            IF ( k1.NE.0 .OR. k2.NE.0 ) ELMT = ELMT*SQRT(2.)
            RETURN
         ENDIF
C-----K=0
      ELSEIF ( k1.NE.0 ) THEN
C-----In band, K.ne.0
         pha1 = (-1.)**((i1-llam+k1)/2)
         pha2 = (-1.)**((k1+i1)/2+1)*pha1
         ELMT = fac*(pha1*WTHREJ(i1,llam,i2,k1,0,-k1)
     &          *Xm1+pha2*WTHREJ(i1,llam,i2,-k1,2*k1,-k1)*Xm2)
         RETURN
      ELSE
         ELMT = fac*WTHREJ(i1,llam,i2,0,0,0)*Xm1
         RETURN
      ENDIF
      END
      SUBROUTINE OPENF1
      IMPLICIT NONE
      INTEGER*4 i , j , k
      CHARACTER name*60 , opt1*20 , opt2*20
 100  READ * , i , j , k
      IF ( i.EQ.0 ) RETURN
      IF ( j.EQ.1 ) opt1 = 'OLD'
      IF ( j.EQ.2 ) opt1 = 'NEW'
      IF ( j.EQ.3 ) opt1 = 'UNKNOWN'
      IF ( k.EQ.1 ) opt2 = 'FORMATTED'
      IF ( k.EQ.2 ) opt2 = 'UNFORMATTED'
      READ 99001 , name
99001 FORMAT (A)
      OPEN (i,IOSTAT=k,FILE=name,STATUS=opt1,FORM=opt2)
c      IF (K.EQ.0) WRITE(6,1030) 'OPENED ',NAME
c 1030 FORMAT (1X,2A)
c      WRITE(6,1010) ' IO-num = ',I,OPT1,OPT2
c 1010 FORMAT (1X,A,I4,2(1x,A))
      IF ( k.EQ.0 ) GOTO 100
c      WRITE (6,1020) 'PROBLEMS OPENING ',NAME,K
c 1020 FORMAT(A,A,I6)
      END
C-----------------------------------------------------------------------
