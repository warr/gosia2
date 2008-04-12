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
      LOGICAL ERR
      COMPLEX ARM,EXPO
      character*4 oph,op1,opcja,op2
      character*1 prp
      DIMENSION IHLM(32),ESP(20),DEDX(20),BTEN(1200)
     *,FIEX1(11,20,2),TITLE(20),PFI(101),ZMIR(6,2,50)
     *,IECD(50),WPI(11,2),TAU1(10),ENG(10),TAU2(10,7),XL1(7),
     *QUI(8,10),CF(8,2),IVARH(500),liscl(200),dsxm(100,20,20),levl(50)
     *,xlevb(50,2),bm(8,20,20,3),mlt(500),ivari(500),jpin(50)
      common/clust/iclust(50,200),lastcl(50,20),sumcl(20,500)
     *,irawex(50)
      common/cccds/ndst(50)
	common/switch/jzb,ibps
      common/inhi/inhb
      COMMON/IDENT/BEQ
      common/efcal/abc(8,10),akavka(8,200),thick(200,7)
      COMMON/PCOM/PSI(500)
      COMMON/TCM/TETACM(50),TREP(50),DSIGS(50)
      COMMON/BREC/BETAR(50)
      COMMON/ADBXI/EXPO(500)
      COMMON/DIMX/DIX(4),ODL(200)
      COMMON/TRA/DELTA(500,3),ENDEC(500),ITMA(50,200),ENZ(200)
      COMMON/CINIT/cnor(32,75),innr,mcfix,cnor1(32,75),cnor2(32,75)
      COMMON/XRA/SE
      COMMON/HHH/HLM(500)
      COMMON/VAC/VACDP(3,75),QCEN,DQ,XNOR,AKS(6,75),IBYP
      COMMON/LIFE/NLIFT
      COMMON/MIXD/DMIXE(20,2),DMIX(20),IMIX(20),NDL
      COMMON/ME2D/NAMX,IAMX(100),IAMY(100,2),EAMX(100,2)
      COMMON/LIFE1/LIFCT(50),TIMEL(2,50)
      COMMON/DFTB/DEVD(500),DEVU(500)
      COMMON/ERRAN/KFERR
	common/resc/elm25(500),elm26(500)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/SECK/ISKIN(50)
      COMMON/VLIN/XV(51),YV(51),ZV(20),DSG(20),DSE(20),DS
      COMMON/DUMM/GRAD(500),HLMLM(500),ELMH(500)
      COMMON/BRNCH/BRAT(50,2),IBRC(2,50),NBRA
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      COMMON/YTEOR/YGN(500),YGP(500),IFMO
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/MAP/PARX(50,12,5),PARXM(50,4,10,6),XIR(6,50)
      COMMON/CCC/eg(50),cc(50,5),NANG(200),Q(3,200,8),NICC,
     *AGELI(50,200,2)
      COMMON/GGG/G(7)
	common/tutaj/part(32,50,2),partl(32,50,2)
      COMMON/AZ/ARM(600,7)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/CXI/XI(500)
      COMMON/ALLC/LOCQ(8,7)
      COMMON /CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON /COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/MINNI/IMIN,LNORM(50)
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON /CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/PRT/IPRM(20)
      COMMON /CCOUP/ZETA(50000),LZETA(8)
      COMMON /CB/B(20)
      COMMON /CLM/LMAX
      COMMON/CLCOM0/IFAC(75)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/CLCOM9/ERR
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/CEXC9/INTERV(50)
      COMMON/CAUX0/NCM,EMMA(75)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      COMMON/WARN/SGW,SUBCH1,SUBCH2,IWF
      COMMON/THTAR/ITTE(50)
      COMMON/FIT/LOCKF,NLOCK,IFBFL,LOCKS,DLOCK
      COMMON/APRX/LERF,IDIVE(50,2)
      COMMON/SKP/JSKIP(50)
      COMMON/TRB/ITS
      COMMON/SEL/KVAR(500)
      COMMON/ERCAL/JENTR,ICS
      COMMON/LOGY/LNY,INTR,IPS1
      COMMON/FAKUL/IP(26),IPI(26),KF(101,26),PILOG(26)
      DATA(ENG(K),K=1,10)/.05,.06,.08,.1,.15,.2,.3,.5,1.,1.5/
      DATA(TAU1(K),K=1,10)/17.656,10.726,5.076,2.931,1.3065,
     *.8828,.5959,.4357,.3041,.2472/
      DATA(TAU2(K,1),K=1,10)/.9883,.7473,.5442,.4592,.3718,
     *.3302,.2814,.2278,.1657,.1350/
      DATA(TAU2(K,2),K=1,10)/1.014,.7443,.5195,.4261,.3362,
     *.2967,.2518,.2038,.1479,.1204/
      DATA(TAU2(K,3),K=1,10)/15.167,9.405,4.652,2.889,1.525,
     *1.135,.8643,.6592,.4703,.3830/
      DATA(TAU2(K,4),K=1,10)/23.184,14.182,6.777,4.059,1.970,
     *1.384,.9936,.7473,.5274,.4297/
      DATA(TAU2(K,5),K=1,10)/84.351,51.445,23.822,13.070,
     *4.774,2.605,1.339,.7925,.5005,.4032/
      DATA(TAU2(K,6),K=1,10)/93.364,58.559,125.96,70.713,
     *25.302,12.541,5.193,2.215,1.077,.8176/
      DATA(TAU2(K,7),K=1,10)/89.809,56.338,27.009,62.966,
     *22.933,11.334,4.540,1.813,.8020,.5900/
	read*,ibps
	read 997,op1,op2
	call openf1
	ibps=ibps-1
	jzb=25
	nawr=0
	if(ibps.eq.1)jzb=26
      DO 8367 I=1,75
      DO 8367 J=1,32
	cnor1(j,i)=1.
	cnor2(j,i)=1.
 8367 CNOR(J,I)=1.
      mcfix=1
	mres1=0
	mrepf=0
	mret=-1
	chir=0.
	mawr=0
 1408 continue
      do 8812 i=1,50
 8812 lifct(i)=0
      chp=chir 
	if(mret.eq.1.and.jzb.eq.26)chir=0.
      rewind 25
	rewind 26
	read(jzb,*)ibps
	ibps=ibps-1
      jzb=25
      if(ibps.eq.1)jzb=26
      IBYP=0
      IP(1)=2
      IP(2)=3
      IP(3)=5
      IP(4)=7
      IP(5)=11
      IP(6)=13
      IP(7)=17
      IP(8)=19
      IP(9)=23
      IP(10)=29
      IP(11)=31
      IP(12)=37
      IP(13)=41
      IP(14)=43
      IP(15)=47
      IP(16)=53
      IP(17)=59
      IP(18)=61
      IP(19)=67
      IP(20)=71
      IP(21)=73
      IP(22)=79
      IP(23)=83
      IP(24)=89
      IP(25)=97
      IP(26)=101
      inhb=0
      BEQ=-983872.
      ipinf=0
      IYR=0
      PI=3.141592654
      INNR=0
      ITNO=0
      CHISQ=0.
      CHILO=0.
      IWF=1
      IFM=0
      IPS1=11
      IFWD=-1
      INTR=0
      LNY=0
      JENTR=0
      LP0=50000
      ICS=0
      LP1=50
      LP2=500
      LP3=75
      LP4=1500
      LP6=32
      LP7=LP0-4900
      LP8=LP3*28+1
      LP9=LP0-LP3*28
      LP10=600
      LP11=LP8-1
      LP12=365
      LP13=LP9+1
      LP14=4900
      DO 1676 I=1,LP1
      jpin(i)=0
 1676 IECD(I)=0
      TXX=0.
      SGW=3.
      SUBCH1=0.
      SUBCH2=0.
      ITS=0
      IOSR=0
      LOCKS=0
      DLOCK=1.1
      KERF=0
      IFBFL=0
      NLOCK=0
      LOCKF=0
      DO 3709 I=1,LP4
      DO 3709 J=1,LP6
 3709 CORF(I,J)=1.
      DO 3509 I=1,20
      IPRM(I)=1
      DO 3509 J=1,5
 3509 CC(I,J)=0.
      iprm(2)=0
      IPRM(4)=0
      IPRM(5)=1
      IPRM(6)=11111
      IPRM(7)=0
      IPRM(16)=0
      IPRM(17)=0
      IPRM(18)=0
      IPRM(19)=0
      IPRM(20)=0
	iprm(11)=0
	iprm(12)=0
	iprm(14)=0
      DO 3007 I=1,LP1
      DO 3008 J=1,5
      IF(J.EQ.5)GO TO 3031
      DO 3009 K=1,10
      DO 3009 KUKU=1,6
 3009 PARXM(I,J,K,KUKU)=0.
 3031 CONTINUE
      DO 3010 K=1,12
 3010 PARX(I,K,J)=0.
 3008 CONTINUE
 3007 CONTINUE
      DO 3011 K=1,LP1
      IDIVE(K,1)=1
      IDIVE(K,2)=1
      DO 3011 IUY=1,6
      XIR(IUY,K)=0.
 3011 CONTINUE
      IOBL=0
      LFAGG=0
      IZCAP=12800
      KFERR=0
      NDIM=LP3
      ISO=1
      B(1)=1.
      DO 80 I=2,20
 80   B(I)=B(I-1)*(I-1)
      LMAXE=0
      CALL FAKP
      CALL FHIP
      NCM=2
      DO 1282 IJX=1,LP1
 1282 INTERV(IJX)=1
      LA=0
      IPO3=1
      INDX=0
      ACCUR=.00001
      ICG=1
      IENT=1
      JPHD=1
      DIPOL=0.005
      MAGEXC=0
      LAMMAX=0
      DO 160 LAM=1,8
      DO 150 LEXP=1,LP3
 150  LDNUM(LAM,LEXP)=0
      MULTI(LAM)=0
 160  LAMDA(LAM)=0
      DO 200 J=1,LP2
      EXPO(J)=(1.,0.)
      KVAR(J)=1
 200  ELM(J)=0.
      DO 201 J=1,LP1
      JSKIP(J)=1
 201   ISKIN(J)=0
      DO 3810 J=1,LP3
 3810  ISEX(J)=1111
      ISEX(1)=0
      ACCA=.00001
       OPH='    '
      NMEMX=LP2+9
      IEXP=1
      IMIN=0
      I122=0
      DO 1741 J=1,LP2
      DO 1741 K=1,2
      DO 1741 L=1,7
 1741 QAPR(J,K,L)=0.
      ERR=.FALSE.
      INTEND=0
 994  FORMAT(1H1/1X,125(1H*)/1X,125(1H*)/1X,50(1H*),25X,50(1H*)/
     *1X,50(1H*),10X,6HGOSIA2,10X,50(1H*)/1X,50(1H*),25X,50(1H*)/
     *1X,125(1H*)/1X,125(1H*)////)
 2713 FORMAT(1X/20X,43HROCHESTER COULOMB EXCITATION DATA ANALYSIS
     *,37HCODE BY T.CZOSNYKA,D.CLINE AND C.Y.WU /50X,
     *27HLATEST REVISION- JUNE  2006//////)
 1    READ(jzb,997)OP1,OP2
      IF(OP1.NE.'OP, ')GO TO 7893
      IF(OP2.EQ.'GOSI')OPH=OP2
      IF(OP2.EQ.'GOSI')OPCJA=OP2
	if(OP2.eq.'FILE')call openf
      IF(OP2.EQ.'FILE')GO TO 1
      IF(JPHD.EQ.1)WRITE(22,994)
      IF(JPHD.EQ.1)WRITE(22,2713)
      JPHD=0
      IF(OP2.EQ.'GDET')GO TO 5200
      IF(OP2.EQ.'RAND')GO TO 331
      IF(OP2.EQ.'TROU')GO TO 9998
      IF(OP2.EQ.'REST')GO TO 330
      IF(OP2.EQ.'RE,A')GO TO 318
      IF(OP2.EQ.'RE,F')GO TO 318
      IF(OP2.EQ.'ERRO')GO TO 2
      IF(OP2.EQ.'RE,C')GO TO 319
      IF(OP2.EQ.'TITL')GOTO 899
      IF(OP2.EQ.'GOSI')GO TO 900
      IF(OP2.EQ.'COUL')GO TO 900
      IF(OP2.EQ.'EXIT')GO TO 430
      IF(OP2.EQ.'MINI')GO TO 673
      if(OP2.EQ.'THEO')go to 2913
      IF(OP2.EQ.'YIEL')GO TO 950
      IF(OP2.EQ.'INTG')GO TO 471
      IF(OP2.EQ.'CORR')GO TO 310
      IF(OP2.EQ.'POIN')GO TO 805
      IF(OP2.EQ.'MAP ')IOBL=1
      IF(OP2.EQ.'STAR')GO TO 805
      IF(OP2.EQ.'SIXJ')GO TO 1410
      IF(OP2.EQ.'RAW ')GO TO 30
      IF(OP2.EQ.'MAP ')GO TO 805
 7893 WRITE(22,7894)OP1,OP2
 7894 FORMAT(5X,19HUNRECOGNIZED OPTION,1X,1A3,1A4)
      GO TO 9999
 900  READ(jzb,996)OP1
      IF(OP1.EQ.'    ')GO TO 1
      IF(OP1.EQ.'LEVE')GO TO 800
      IF(OP1.EQ.'ME  ')GO TO 801
      IF(OP1.EQ.'CONT')GO TO 803
      IF(OP1.EQ.'EXPT')GO TO 804
      WRITE(22,7895)OP1
 7895 FORMAT(5X,22HUNRECOGNIZED SUBOPTION,1X,1A4)
      GO TO 9999
 2913 rewind(12)
      ibaf=1
      do 2915 jb=1,lp1
      do 2915 lb=1,2
      xlevb(jb,lb)=0
 2915 continue     
      read(jzb,*)nbands
      if(nbands.le.0)ibaf=0
      nbands=iabs(nbands) 
      do 2920 nl=1,8
      do 2920 jb=1,nbands
      do 2920 jl=1,nbands
      do 2920 kl=1,3
 2920 bm(nl,jb,jl,kl)=0.     
      do 2914 jb=1,nbands
      read(jzb,*)bk,ilevls
      read(jzb,*)(levl(ib),ib=1,ilevls)
      do 2916 kb=1,ilevls
      inva=levl(kb)
      xlevb(inva,2)=bk
 2916 xlevb(inva,1)=float(jb)               
 2914 continue
      do 2955 nl=1,8
      read(jzb,*)nnl
 2949 continue
      if(nnl.le.0)go to 2917 
      read(jzb,*)jb1,jb2
      if(jb1.eq.0)go to 2955 
      read(jzb,*)(bm(nnl,jb1,jb2,j),j=1,3)
      do 2981 j=1,3
 2981 bm(nnl,jb2,jb1,j)=bm(nnl,jb1,jb2,j)
      go to 2949
 2955 continue
 2917 continue
      do 2918 kb=1,memax
      if(ibaf.eq.0)go to 2918
      ind1=lead(1,kb)
      ind2=lead(2,kb)
      xi1=spin(ind1)
      xi2=spin(ind2)
      lamd=mlt(kb)
      nb1=ifix(xlevb(ind1,1)+.1)            
      nb2=ifix(xlevb(ind2,1)+.1)
      xk1=xlevb(ind1,2)
      xk2=xlevb(ind2,2)
      xm1=bm(lamd,nb1,nb2,1)
      xm2=bm(lamd,nb1,nb2,2)
      xm3=bm(lamd,nb1,nb2,3)
      elm(kb)=elmt(xi1,xi2,lamd,nb1,nb2,xk1,xk2,xm1,xm2,xm3)
      if(abs(elm(kb)).lt.1e-6)elm(kb)=1.e-6
	irix=12
	if(ibps.eq.1)irix=32
      write(irix,*)elm(kb)
 2918 continue
      go to 1     
 30   rewind 8
      do 31 l=1,8
      read(8,*)(abc(l,j),j=1,10)
      do 32 j=1,10
  32  abc(l,j)=alog(abc(l,j))
  31  continue
      do 33 l=1,nfd
  33  read(8,*)(thick(l,j),j=1,7)
      do 34 l=1,lp1
      do 38 j=1,200
  38  iclust(l,j)=0
      do 39 j=1,20
  39  lastcl(l,j)=0
  34  irawex(l)=0
      do 35 l=1,lp1
      read(jzb,*)mexl
      if(mexl.eq.0)go to 1
      irawex(mexl)=1
      n=nang(mexl)
      do 36 j=1,n
      jj=itma(mexl,j)
      read(jzb,*)(akavka(k,jj),k=1,8)
  36  continue
      read(jzb,*)kclust
      if(kclust.eq.0)go to 35
      do 37 j=1,kclust
      read(jzb,*)numcl
      read(jzb,*)(liscl(k),k=1,numcl)
      lastcl(l,j)=liscl(numcl)
      do 41 k=1,numcl
      kk=liscl(k)
      iclust(l,kk)=j
   41 continue
   37 continue
   35 continue
      go to 1
 5200 CONTINUE
      NL=7
      read(jzb,*)NFDD
      nfd=iabs(nfdd)
      if(nfdd.gt.0)go to 5208
      rewind 8
      do 5209 i=1,nl
 5209 write(8,*)(tau2(l,i),l=1,10)
      write(8,*)(eng(l),l=1,10)
 5208 REWIND 9
      WRITE(9,*)NFD
      DO 5201 I=1,NFD
      read(jzb,*)(DIX(K),K=1,4)
      read(jzb,*)(XL1(K),K=1,NL)
      IF(DIX(1).le.0.)DIX(1)=.01
      WRITE(9,*)DIX(4)
      if(nfdd.gt.0)go to 5207
      write(8,*)(xl1(k),k=1,nl)
 5207 continue
      IND=1
      IF(XL1(5).gt.0.)IND=3
      IF(XL1(6).gt.0.)IND=4
      IF(XL1(7).gt.0.)IND=5
      WRITE(9,*)ENG(IND)
      CALL QFIT(QUI,TAU1,TAU2,ENG,XL1,CF,NL,IND)
      WRITE(22,5206)I
      DO 5202 K=1,8
      WRITE(22,5203)K,CF(K,1),CF(K,2)
      WRITE(9,*)CF(K,1),CF(K,2),QUI(K,IND)
      DO 5204 L=1,10
      ARG=(ENG(L)-ENG(IND))**2
      QC=(QUI(K,IND)*CF(K,2)+CF(K,1)*ARG)/(CF(K,2)+ARG)
 5204 WRITE(22,5205)ENG(L),QC,QUI(K,L),100.*(QC-QUI(K,L))/QUI(K,L)
 5202 CONTINUE
 5201 CONTINUE
 5206 FORMAT(10X,8HDETECTOR,1X,1I2)
 5203 FORMAT(1X,//5X,2HK=,1I1,2X,3HC1=,1E14.6,2X,3HC2=,1E14.6/
     *5X,11HENERGY(MEV),5X,9HFITTED QK,5X,7HCALC.QK,5X,
     *8HPC.DIFF./)
 5205 FORMAT(8X,1F4.2,6X,1F9.4,5X,1F9.4,3X,1E10.2)
      GO TO 1
 1410 CONTINUE
      DO 1411 K=1,2
      L=4*K
      DO 1412 J=1,80
      IXJ=J-1
      DO 1413 MS=1,5
      MEND=2*(MS-3)+IXJ
      WRITE(14,*)WSIXJ(L,4,4,IXJ,MEND,IXJ-4),
     *WSIXJ(L,4,4,IXJ,MEND,IXJ-2),
     *WSIXJ(L,4,4,IXJ,MEND,IXJ),WSIXJ(L,4,4,IXJ,MEND,IXJ+2),
     *WSIXJ(L,4,4,IXJ,MEND,IXJ+4)
 1413 CONTINUE
 1412 CONTINUE
 1411 CONTINUE
      GO TO 9999
 2    read(jzb,*)IDF,MS,MEND,IREP,IFC,REMAX
      mcfix=0
	do 7321 mmmm=1,32
	do 7321 kkkk=1,50
 7321 read(13,*)cnor1(mmmm,kkkk)	
      REM=ALOG(REMAX)
      LOCKS=0
      LOCKF=0
      JENTR=1
      SH=1.
      IFBP=0
      INPO=1
      INKO=1
      IF(IOSR.EQ.0.OR.IDF.EQ.0)GO TO 3212
      INN=0
      IJ=MULTI(1)
      IF(IJ.EQ.0)GO TO 3213
      DO 3214 IJ=1,NMAX
      LXD=LDNUM(1,IJ)
      IF(LXD.EQ.0)GO TO 3214
      DO 3215 IJK=1,LXD
 3215 INN=INN+1
 3214 CONTINUE
      INPO=INN+1
 3213 CONTINUE
      DO 3216 IJ=1,NMAX
      LXD=LDNUM(2,IJ)
      IF(LXD.EQ.0)GO TO 3216
      DO 3218 IJK=1,LXD
 3218 INN=INN+1
 3216 CONTINUE
      INKO=INN
      IF(IREP.EQ.2)GO TO 3212
      WRITE(33,*)NMAX,MEMAX,INPO,INKO
      DO 3273 INN=1,NMAX
 3273 WRITE(33,*)INN,SPIN(INN),EN(INN)
      DO 3219 INN=1,MEMAX
 3219 WRITE(33,*)INN,LEAD(1,INN),LEAD(2,INN)
      DO 3220 INN=1,MEMAX
 3220 WRITE(33,*)INN,ELM(INN)
 3212 CONTINUE
      IF(IREP.NE.0)GO TO 4222
      DO 91 KH1=1,MEMAX
      DEVD(KH1)=ELML(KH1)-ELM(KH1)
 91   DEVU(KH1)=ELMU(KH1)-ELM(KH1)
      GO TO 4223
 4222 REWIND 15
      READ(15,*)(DEVD(KH1),DEVU(KH1),KH1=1,MEMAX)
 4223 CONTINUE
      IF(IMIN.EQ.0)CALL CMLAB(0,DSIG,TTTTT)
      IF(ERR)GO TO 9999
      IF(IMIN.EQ.0)GO TO 610
  7   CONTINUE
      IF(ICS.EQ.1)GO TO 7227
      CALL FTBM(0,CHISS,IDR,0,CHILO,BTEN)
      REWIND 11
      DO 3712 KH1=1,LP4
 3712 write(11)(CORF(KH1,KH2),KH2=1,LP6)
      GO TO 3917
 7227 REWIND 11
      DO 3711 KH1=1,LP4
 3711 read(11)(CORF(KH1,KH2),KH2=1,LP6)
 3917 CALL FTBM(3,CHISS,IDR,1,CHILO,BTEN)
      CHIS0=CHISS
      WRITE(22,1912)CHIS0
      inhb=1
 1912 FORMAT(1X///10X,20H***** CENTRAL CHISQ=,1E12.4,1X,5H*****//)
      CHISL=CHISS
      DO 8 KH=1,MEMAX
 8    HLM(KH)=ELM(KH)
      IF(IDF.EQ.1)GO TO 4000
      NAXFL=0
      IF(MS.EQ.0)MEND=MEMAX
      IF(MS.EQ.0)MS=1
      DO 4002 KH=MS,MEND
      DO 4271 IJ=1,2
       PV=(ELMU(KH)-ELML(KH))/100.
       IF(IJ.EQ.1.AND.(ELM(KH)-ELML(KH)).LT.PV)GO TO 4271
       IF(IJ.EQ.2.AND.(ELMU(KH)-ELM(KH)).LT.PV)GO TO 4271
      DO 4003 KH1=1,MEMAX
 4003 SA(KH1)=0.
      IF(IVAR(KH).EQ.0)GO TO 4002
      SA(KH)=1.*(-1)**IJ
      KH1=KH
      CALL KONTUR(IDR,CHIS0,CHISL,IFBP,-1,KH1,SH,BTEN,REM)
      ELM(KH)=HLM(KH)
 4271 CONTINUE
      REWIND 15
      WRITE(15,*)(DEVD(IJ),DEVU(IJ),IJ=1,MEMAX)
 4002 CONTINUE
 3235 IF(IFBP.NE.1)GO TO 65
      REWIND 17
      DO 1645 LKJ=1,MEMAX
      READ(17,*)ELM(LKJ)
 1645 CONTINUE
      WRITE(22,6828)
 6828 FORMAT(1X///20X,33H*** BEST POINT FOUND (TAPE17) ***///)
      CALL PRELM(3)
  65  CONTINUE
      IF(NAXFL.EQ.0)WRITE(22,6662)
      IF(NAXFL.NE.0)WRITE(22,6673)
 6673 FORMAT(1X///44X,7HOVERALL)
 6662 FORMAT(1X///43X,8HDIAGONAL)
      WRITE(22,66)
 66   FORMAT(40X,16HESTIMATED ERRORS//5X,
     *5HINDEX,5X,2HNI,5X,2HNF,5X,13HME AND ERRORS//)
      DO 67 KH1=1,MEMAX
      IF(IVAR(KH1).EQ.0.OR.IVAR(KH1).GT.999)GO TO 67
      WRITE(22,68)KH1,LEAD(1,KH1),LEAD(2,KH1),HLM(KH1),DEVD(KH1),
     *DEVU(KH1),DEVD(KH1)*100./ABS(HLM(KH1)),
     *DEVU(KH1)*100./ABS(HLM(KH1))
 68   FORMAT(6X,1I3,6X,1I2,5X,1I2,5X,1F9.5,2X,1H(,1F9.5,2H ,,1F9.5,
     *1H),6H......,1F7.1,2H ,,1F7.1,1X,2HPC)
 67   CONTINUE
      IF(NAXFL.NE.0)WRITE(22,6673)
      IF(NAXFL.EQ.0)WRITE(22,6662)
      WRITE(22,6692)
 6692 FORMAT(40X,16HESTIMATED ERRORS,//5X,
     *5HINDEX,5X,2HNI,5X,2HNF,5X,29HB(E,ML)(OR QUADRUPOLE MOMENT)
     *,11H AND ERRORS//)
      DO 6695 KH2=1,MEMAX
      IF(IVAR(KH2).EQ.0.OR.IVAR(KH2).GT.999) GO TO 6695
      ISPA=LEAD(2,KH2)
      IF(LEAD(1,KH2).NE.LEAD(2,KH2)) GO TO 6697
      ISPB=IFIX(SPIN(ISPA))*2
      QFAC=3.170662*WTHREJ(ISPB,4,ISPB,-ISPB,0,ISPB)
      WRITE(22,6694) KH2,LEAD(2,KH2),LEAD(1,KH2),HLM(KH2)*QFAC,
     *DEVD(KH2)*QFAC,DEVU(KH2)*QFAC
      GO TO 6695
 6697 SBE=2.*SPIN(ISPA)+1.
      BE2=HLM(KH2)*HLM(KH2)/SBE
      BE2A=HLM(KH2)+DEVD(KH2)
      BE2B=HLM(KH2)+DEVU(KH2)
      BE2C=BE2B
      IF(ABS(BE2A).GT.ABS(BE2B)) BE2B=BE2A
      IF(abs(be2a-be2c).lt.1.e-6) BE2A=BE2C
      IF(BE2A/HLM(KH2).LE.0..OR.BE2B/HLM(KH2).LE.0.) BE2A=0.
      BE2A=BE2A**2/SBE
      BE2B=BE2B**2/SBE
      WRITE(22,6694) KH2,LEAD(2,KH2),LEAD(1,KH2),BE2,BE2A-BE2,BE2B-BE2
 6694 FORMAT(6X,1I3,6X,1I2,5X,1I2,5X,1F10.5,2X,1H(,1F10.5,2H ,,1F10.5,
     *1H))
 6695 CONTINUE
      GO TO 9999
 4000 CONTINUE
      IFBFL=1
      IF(IREP.NE.2)GO TO 4057
      IF(IOSR.EQ.0)GO TO 4057
      REWIND 33
      READ(33,*)LL,MM,KK,INN
      DO 4058 INN=1,LL
 4058 READ(33,*)MM,YYY,ZZ
      DO 4059 INN=1,MEMAX
 4059 READ(33,*)MM,LL,KK
      DO 4060 INN=1,MEMAX
 4060 READ(33,*)MM,YYY
 4063 READ(33,*)MM,LL
      IF(MM.EQ.0)GO TO 4061
      READ(33,*)KK,LL,YYY
      READ(33,*)(SA(MM),MM=1,MEMAX)
      GO TO 4063
 4061 BACKSPACE 33
 4057 CONTINUE
      IREA=0
      IF(MS.LT.0)IREA=1
      IF(MS.EQ.0)MEND=MEMAX
      IF(MS.EQ.0)MS=1
 4127 NAXFL=1
      IF(IREA.EQ.1)read(jzb,*)MS,MEND
      IF(MS.EQ.0)GO TO 4128
      DO 4510 KH=MS,MEND
      IF(IFC.EQ.1)GO TO 4552
      REWIND 18
      DO 4821 KH1=1,KH
 4821 READ(18,*)(KVAR(JYI),JYI=1,MEMAX)
      DO 4621 KH1=1,MEMAX
      IVRH=IVAR(KH1)
      IF(KVAR(KH1).EQ.0)IVAR(KH1)=0
 4621 KVAR(KH1)=IVRH
 4552 CONTINUE
      DO 4553 IJ=1,2
      SH=DEVU(KH)
      IF(IJ.EQ.1)SH=DEVD(KH)
      IF(abs(SH).lt.1.e-6)SH=(-1)**IJ*ABS(HLM(KH))/10.
      ELM(KH)=HLM(KH)+1.5*SH
      MM=0
      DO 4129 KH1=1,MEMAX
      if(ifc.eq.1)kvar(kh1)=ivar(kh1)
 4129 MM=MM+IVAR(KH1)
      IF(MM.EQ.0)WRITE(22,4131)KH
      IF(MM.EQ.0)GO TO 3232
 4131 FORMAT(10X,3HME=,1I3,5X,
     *23HNO FREE MATRIX ELEMENTS)
      KFERR=1
      IF(IOSR.EQ.1)WRITE(33,*)KH,KH
      IF(IOSR.EQ.1)WRITE(33,*)KH,IJ,ELM(KH)
      LOCKS=1
      DLOCK=.05
      CALL MINI(CHISS,-1.,2,.0001,1000,IDR,100000.,0,IOSR,KH,BTEN)
      DO 4511 KH1=1,MEMAX
 4511 SA(KH1)=(ELM(KH1)-HLM(KH1))/ABS(SH)
      CALL KONTUR(IDR,CHIS0,CHISL,IFBP,INPO,KH,SH,BTEN,REM)
 3232 CONTINUE
      DO 4512 KH1=1,MEMAX
      if(ifc.eq.1)ivar(kh1)=kvar(kh1)
 4512 ELM(KH1)=HLM(KH1)
 4553 CONTINUE
      IF(IFC.EQ.1)GO TO 4555
      DO 4554 KH1=1,MEMAX
 4554 IVAR(KH1)=KVAR(KH1)
 4555 CONTINUE
      REWIND 15
      WRITE(15,*)(DEVD(KH1),DEVU(KH1),KH1=1,MEMAX)
 4510 CONTINUE
      IF(IREA.EQ.1)GO TO 4127
 4128 CONTINUE
      IF(IOSR.EQ.0)GO TO 3234
      IM=0
      WRITE(33,*)IM,IM
 3234 GO TO 3235
 318  JFRE=0
      IRFIX=0
      IF(OP2.EQ.'RE,F')IRFIX=1
      GO TO 321
 319  JFRE=1
      IRFIX=0
 321  DO 322 JRLS=1,MEMAX
      IF(IVAR(JRLS).EQ.0.AND.IRFIX.EQ.1)GO TO 322
      IF(IVAR(JRLS).GT.999)GO TO 323
 324  IVAR(JRLS)=2
      ELML(JRLS)=-ABS(ELML(JRLS))
      ELMU(JRLS)=ABS(ELMU(JRLS))
      IF(JRLS.GT.MEMX6)IVAR(JRLS)=1
      GO TO 322
 323  IF(JFRE.EQ.1)GO TO 322
      GO TO 324
 322  CONTINUE
      DO 1895 JRLS=1,MEMAX
 1895 IVARH(JRLS)=IVAR(JRLS)
      GO TO 1
 331  read(jzb,*)SE
      CALL MIXUP
      WRITE(22,317)
 317  FORMAT(1X///5X,29HMATRIX ELEMENTS RANDOMIZED...///)
      CALL PRELM(2)
      GO TO 1
 330  irix=12
      if(ibps.eq.1)irix=32
      REWIND irix
      memax1=memax+1
      DO 5296 LKJ=1,memax
 5296 READ(irix,*)ELM(LKJ)
      DO 5298 LKJ=1,memax1
      read(jzb,*)LKJ1,XLK
      IF(LKJ1.EQ.0)GO TO 5299
      if(mres1.eq.0)ELM(LKJ1)=XLK
 5298 CONTINUE
 5299 CONTINUE
      do 8331 lkj=1,memax
	if(mres1.eq.0.and.jzb.eq.25)elm25(lkj)=elm(lkj)
	if(mres1.eq.0.and.jzb.eq.26)elm26(lkj)=elm(lkj)
 8331 continue	
      if(jzb.eq.26)mres1=1
      WRITE(22,5523)
 5523  FORMAT(1X///5X,5H*****,2X,
     *35HRESTART-MATRIX ELEMENTS OVERWRITTEN,2X,5H*****///)
       do 7933 kk=1,memax
       la=mlt(kk)
       IF(IVARI(KK).LT.10000)GO TO 7933
      KK1=IVARI(KK)/10000
      KK2=IVARI(KK)-10000*KK1
       LA1=LA
       IF(KK2.LT.100)GO TO 8428
       LA1=KK2/100
       KK2=KK2-100*LA1
 8428  continue
      INX1=MEM(KK1,KK2,LA1)
C      ELML(KK)=ELML(INX1)*ELM(KK)/ELM(INX1)
C      ELMU(KK)=ELMU(INX1)*ELM(KK)/ELM(INX1)
      SA(KK)=ELM(KK)/ELM(INX1)
      IVAR(KK)=1000+INX1
      IF(ELMU(KK).GT.ELML(KK))GO TO 7933
      ELMI=ELMU(KK)
      ELMU(KK)=ELML(KK)
      ELML(KK)=ELMI
 7933  CONTINUE
      CALL PRELM(4)
      GO TO 1
 471  CONTINUE
      REWIND 14
      LFAGG=1
      IF(SPIN(1).lt..25)ISO=0
      DO 472 LX=1,NEXPT
      lpin=1
      if(ipinf.eq.0)go to 1472
      if(jpin(lx).eq.0)go to 1472
      lpin=jpin(lx)
 1472 continue
      IEXP=LX
      TTH=TLBDG(LX)
      ENH=EP(LX)
      do 1473 mpin=1,lpin
      IF(IECD(LX).EQ.1) GO TO 1677
      read(jzb,*)NE,NTT,EMN,EMX,TMN,TMX
      GO TO 1678
 1677 read(jzb,*)NE,NTT,EMN,EMX,WTH,WPH,WTHH
      MFLA=1
      CALL COORD(WTH,WPH,WTHH,NTT,0,PFI,WPI,TTH,LX,TMN,TMX)
      GO TO 1679
 1678 CONTINUE
       MFLA=0
       IF(NTT.LT.0)MFLA=1
 1679 NTT=IABS(NTT)
      JAN=NANG(LX)
      jan1=ndst(lx)
      if(irawex(lx).eq.0)jan1=jan
      IF(IECD(LX).EQ.1) GO TO 1690
      WRITE(14,*)NE,NTT,EMN,EMX,TMN,TMX,JAN1,TMX,TMX,TMX
      GO TO 1691
 1690 WRITE(14,*)NE,NTT,EMN,EMX,TMN,TMX,JAN1,WTH,WPH,WTHH
 1691 CONTINUE
      read(jzb,*)(XV(I),I=1,NE)
      IF(IECD(LX).EQ.1) GO TO 1692
      read(jzb,*)(YV(I),I=1,NTT)
 1692 CONTINUE
      IF(TTH.LT.0.)ELMH(2*LX-1)=YV(1)
      IF(TTH.LT.0.)ELMH(2*LX)=YV(NTT)
      DO 474 KLOOP=1,NE
      ENB=XV(KLOOP)
      EP(LX)=ENB
      DO 474 KTT=1,NTT
      TTA=SIGN(YV(KTT),TTH)
      IF(IAXS(LX).EQ.0)GO TO 552
      IF(IECD(LX).EQ.1) GO TO 552
      IF(KLOOP.NE.1) GO TO 552
      read(jzb,*)NFI
      read(jzb,*)(FIEX1(KTT,JFI,1),FIEX1(KTT,JFI,2),JFI=1,NFI)
      IF(TTH.GE.0.) GO TO 552
      DO 555 JFI=1,NFI
      FIEX1(KTT,JFI,1)=FIEX1(KTT,JFI,1)+180.
 555  FIEX1(KTT,JFI,2)=FIEX1(KTT,JFI,2)+180.
 552  TLBDG(LX)=TTA
      IF(KLOOP.NE.1) GO TO 1680
      IF(IECD(LX).EQ.0) GO TO 1680
      NFI=1
      FIEX1(KTT,1,1)=WPI(KTT,1)
      FIEX1(KTT,1,2)=WPI(KTT,2)
 1680 CALL CMLAB(LX,DSIG,TETRC)
      IF(ERR)GO TO 9999
      TTING=TLBDG(LX)
      IF(ERR)GO TO 999
      CALL LOAD(LX,1,1,0.,JJ)
      CALL ALLOC(ACCUR)
      CALL SNAKE(LX,ZPOL)
      CALL SETIN
      DO 475 J=1,LMAX
      POLM=FLOAT(J-1)-SPIN(1)
      CALL LOAD(LX,2,1,POLM,JJ)
      CALL STING(JJ)
      CALL PATH(JJ)
      CALL INTG(IEXP)
      CALL TENB(J,BTEN,LMAX)
 475  CONTINUE
      CALL TENS(BTEN)
      CALL DECAY(CCD,0,CCC)
      do 7257 j=1,lp2
      do 7257 ijan=1,20
 7257 sumcl(ijan,j)=0.
      ija0=0
      DO 497 IJAN=1,JAN
      IF(IAXS(LX).EQ.0)NFI=1
      DO 477 JYI=1,IDR
 477  GRAD(JYI)=0.
      TODFI=0.
      DO 476 JFI=1,NFI
      FI0=FIEX1(KTT,JFI,1)/57.2957795
      FI1=FIEX1(KTT,JFI,2)/57.2957795
      GTH=AGELI(IEXP,IJAN,1)
      fm=(fi0+fi1)/2.
      FIGL=AGELI(IEXP,IJAN,2)
      CALL ANGULA(YGN,IDR,1,FI0,FI1,TETRC,GTH,FIGL,IJAN,OP2)
      IF(IFMO.EQ.0)GO TO 7232
      ID=ITMA(IEXP,IJAN)
      D=ODL(ID)
      RX=D*SIN(GTH)*COS(FIGL-FM)-.25*SIN(TETRC)*COS(FM)
      RY=D*SIN(GTH)*SIN(FIGL-FM)-.25*SIN(TETRC)*SIN(FM)
      RZ=D*COS(GTH)-.25*COS(TETRC)
      RL=SQRT(RX*RX+RY*RY+RZ*RZ)
      SF=D*D/RL/RL
      THC=TACOS(RZ/RL)
      FIC=ATAN2(RY,RX)
      CALL ANGULA(YGP,IDR,1,FI0,FI1,TETRC,THC,FIC,IJAN,OP2)
      DO 7233 IXL=1,IDR
      IXM=KSEQ(IXL,3)
      TFAC=TAU(IXM)
 7233 YGN(IXL)=YGN(IXL)+.01199182*TFAC*BETAR(IEXP)*
     *(SF*YGP(IXL)-YGN(IXL))
 7232 CONTINUE
      if(irawex(lx).eq.0)go to 7252
      ipd=itma(lx,ijan)
      do 7253 jyi=1,idr
      ni=kseq(jyi,3)
      nf=kseq(jyi,4)
      decen=en(ni)-en(nf)
      cocos=sin(tetrc)*sin(gth)*cos(fm-figl)+cos(tetrc)*
     *cos(gth)
      decen=decen*(1.+betar(lx)*cocos)
      call effix(ipd,decen,effi)
      ygn(jyi)=ygn(jyi)*effi
 7253 continue
      inclus=iclust(lx,ijan)
      if(inclus.eq.0)go to 7252
      do 7254 jyi=1,idr
 7254 sumcl(inclus,jyi)=sumcl(inclus,jyi)+ygn(jyi)
      if(ijan.ne.lastcl(lx,inclus))go to 497
      do 7255 jyi=1,idr
 7255 ygn(jyi)=sumcl(inclus,jyi)
 7252 continue
      if(jfi.eq.1)ija0=ija0+1
      DO 498 JYI=1,IDR
 498  GRAD(JYI)=GRAD(JYI)+YGN(JYI)
       TODFI=TODFI+ABS(FI1-FI0)
 476  CONTINUE
      IF(IAXS(LX).EQ.0) TODFI=6.283185
      AX=1.
       IF(MFLA.EQ.1)AX=1./TODFI
      DSX=DSIG
      IF(MFLA.NE.1) DSX=DSIG*TODFI
      dsxm(mpin,kloop,ktt)=dsx
	write(17,*)lx,mpin,kloop,ktt,dsx
      WRITE(14,*)LX,ENB,TTING,IJA0,DSX,(GRAD(JYI)*DSIG*AX,JYI=1
     *,IDR)
      IF(IPRM(11).NE.1)GO TO 497
      WRITE(22,478)LX,IJA0,ENB,TTA
      IF(TTA.LT.0.)WRITE(22,387)TTING
 387  FORMAT(5X,35HRESPECTIVE TARGET SCATTERING ANGLE=,1F7.3
     *,1X,3HDEG/)
      DO 479 JYI=1,IDR
      NI=KSEQ(JYI,3)
      NF=KSEQ(JYI,4)
 479  WRITE(22,480)NI,NF,SPIN(NI),SPIN(NF),GRAD(JYI)*DSIG*AX
     *,GRAD(JYI)/GRAD(IDRN)
 478  FORMAT(1X//50X,17HCALCULATED YIELDS//
     *5X,11HEXPERIMENT ,1I2,2X,9HDETECTOR ,1I2/
     *5X,7HENERGY ,1F10.3,1X,3HMEV,2X,6HTHETA ,1F7.3,1X,3HDEG
     *//5X,2HNI,5X,2HNF,5X,2HII,5X,2HIF,5X,5HYIELD,
     *5X,16HNORMALIZED YIELD/)
4478  format(1x//50x,17hCALCULATED YIELDS//
     *5x,11hEXPERIMENT ,1i2,2x,9hDETECTOR ,1i2/
     *5x,7hENERGY ,1f10.3,1x,3hMEV,2x,6hTHETA ,1f7.3,1x,3hDEG
     *//5x,2hNI,5x,2hNF,5x,2hII,5x,2hIF,5x,6hE(MeV),
     *5x,'EFFICIENCY'/)
 480  FORMAT(5X,1I2,5X,1I2,3X,1F4.1,3X,1F4.1,3X,1E11.5,3X
     *,1E11.5)
 497  CONTINUE
 474  CONTINUE
 1473 continue
      EP(LX)=ENH
      TLBDG(LX)=TTH
 472  CONTINUE
      REWIND 14
      REWIND 15
      ISKE=0
      DO 4328 NA=1,LP6
 4328 ILE(NA)=1
      ilx=0
      DO 530 LX=1,NEXPT
	rewind 17
	do 1728 ijaja=1,300000
	read(17,*,END=1729)jjlx,jmpin,jkloo,jktt,dsx
	if(jjlx.eq.lx)dsxm(jmpin,jkloo,jktt)=dsx
 1728 continue
 1729 continue	
      NA=NANG(LX)
      IF(LX.EQ.1)GO TO 4329
      DO 4330 NA1=1,LP6
 4330 ILE(NA1)=ILE(NA1)+NYLDE(LX-1,NA1)
 4329 read(jzb,*)NPTX
      IF(NPTX.EQ.0)GO TO 531
      read(jzb,*)(ESP(I),I=1,NPTX)
      read(jzb,*)(DEDX(I),I=1,NPTX)
      NPT=NPTX
 531  read(jzb,*)NPCE,NPCT
       MFLA=0
      IF(NPCT.LT.0)MFLA=1
      IF(IECD(LX).EQ.1) MFLA=1
      NPCT=IABS(NPCT)
      NPCE=NPCE+MOD(NPCE,2)
      NPCT=NPCT+MOD(NPCT,2)
      mpin=1
      if(ipinf.eq.0)go to 1530
      if(jpin(lx).eq.0)go to 1530
      mpin=jpin(lx)
 1530 continue
      dst=0.
      do 1531 lpin=1,mpin
      ilx=ilx+1
      IF(iLX.EQ.1)GO TO 532
      CALL TAPMA(LX,ISKE,ISKO,ISKF,NFLR,IDR,0,NFT,ENB)
 532  READ(14,*)NE,NTT,EMN,EMX,TMN,TMX,JAN,WTH,WPH,WTHH
      IOCC=(NE+NTT)*IDR
      IF(IOCC.GT.IZCAP)GO TO 598
      HEN=(EMX-EMN)/NPCE
      NPCE1=NPCE+1
      HET=(TMX-TMN)/NPCT
      NPCT1=NPCT+1
      IF(IECD(LX).EQ.1) CALL COORD(WTH,WPH,WTHH,NPCT1,1,PFI,WPI,
     *TLBDG(LX),LX,TMN,TMX)
      IF(IECD(LX).EQ.1) GO TO 1681
       IF(MFLA.EQ.1)read(jzb,*)(PFI(J),J=1,NPCT1)
 1681 HET=HET/57.2957795
      DO 533 J=1,NPCE1
      XX=(J-1)*HEN+EMN
      CALL LAGRAN(ESP,DEDX,NPT,1,XX,YY,3,1)
 533  HLMLM(J)=1./YY
      naa=ndst(lx)
      if(irawex(lx).eq.0)naa=nang(lx)
      ISKF=NAA-1
      DO 534 JA=1,NAA
      ICLL=3
      DO 535 JE=1,NE
      LU=ILE(JA)
      ISKO=(JE-1)*NAA*NTT+JA-1
      CALL TAPMA(LX,ISKE,ISKO,ISKF,NTT,IDR,1,NFT,ENB)
      IF(NFT.EQ.1)GO TO 999
      DO 536 JD=1,IDR
      DO 537 JTP=1,NTT
      IF(JD.EQ.1.AND.JA.EQ.1)DSG(JTP)=dsxm(lpin,je,jtp)
      JYV=(JTP-1)*IDR+JD
 537  YV(JTP)=ZETA(JYV)
      DO 5000 JT=1,NPCT1
      XX=(JT-1)*HET+TMN/57.2957795
      CALL LAGRAN(XV,YV,NTT,JT,XX,YY,2,ICLL)
      CALL LAGRAN(XV,DSG,NTT,JT,XX,ZZ,2,ICLL)
      IF(MFLA.EQ.1)YY=YY*PFI(JT)/57.2957795
      IF(YY.le.0.)YY=1.E-15
      IF(MFLA.EQ.1)ZZ=ZZ*PFI(JT)/57.2957795
      XI(JT)=YY*SIN(XX)
      IF(JD.EQ.1.AND.JA.EQ.1)HLM(JT)=ZZ*SIN(XX)
 5000 CONTINUE
      ICLL=4
      LOCAT=NTT*IDR+(JE-1)*IDR+JD
      ZETA(LOCAT)=SIMIN(NPCT1,HET,XI)
      IF(JD.EQ.1.AND.JA.EQ.1)DSE(JE)=SIMIN(NPCT1,HET,HLM)
      ZV(JE)=ENB
 536  CONTINUE
 535  CONTINUE
      ICLL=3
      DO 539 JD=1,IDR
      DO 540 JTP=1,NE
      JYV=(JTP-1)*IDR+JD+NTT*IDR
 540  YV(JTP)=ZETA(JYV)
      DO 541 JT=1,NPCE1
      XX=(JT-1)*HEN+EMN
      CALL LAGRAN(ZV,YV,NE,JT,XX,YY,2,ICLL)
      IF(JD.EQ.1.AND.JA.EQ.1)CALL LAGRAN(ZV,DSE,NE,JT,XX,ZZ,2,ICLL)
      IF(JD.EQ.1.AND.JA.EQ.1)HLM(JT)=ZZ*HLMLM(JT)
 541  XI(JT)=YY*HLMLM(JT)
      ICLL=4
      IF(JD.EQ.1.AND.JA.EQ.1)DS=SIMIN(NPCE1,HEN,HLM)
 539  GRAD(JD)=SIMIN(NPCE1,HEN,XI)
      if(ja.eq.1)dst=dst+ds
      IF(JA.EQ.1)WRITE(22,5001)DS,LX
      WRITE(22,483)LX,JA,EMN,EMX,TMN,TMX
 5001 FORMAT(1X/////5X,36HINTEGRATED RUTHERFORD CROSS SECTION=,1E9.4,
     *2X,8HFOR EXP.,1I2///)
      DO 7218 JD=1,IDR
 7218 WRITE(15,*)GRAD(JD)
      DO 542 JD=1,IDR
      NI=KSEQ(JD,3)
      NF=KSEQ(JD,4)
 542  WRITE(22,480)NI,NF,SPIN(NI),SPIN(NF),GRAD(JD),
     *GRAD(JD)/GRAD(IDRN)
  534 CONTINUE
      IF(IECD(LX).NE.1) GO TO 1683
      if(jpin(lx).ne.0)go to 1683
      CALL COORD(WTH,WPH,WTHH,1,2,PFI,WPI,TLBDG(LX),LX,TXX,TXX)
      WRITE(22,1682) FIEX(LX,1)*57.2957795,FIEX(LX,2)*57.2957795,LX
      IF(TLBDG(LX).GE.0) GO TO 1683
      FIEX(LX,1)=FIEX(LX,1)+3.14159265
      FIEX(LX,2)=FIEX(LX,2)+3.14159265
 1682 FORMAT(//5X,38HWARNING: THE PHI ANGLE WAS REPLACED BY,1X,
     *F8.3,1X,2HTO,F8.3,3X,14HFOR EXPERIMENT,2X,I3)
 1683 ISKE=ISKE+NE*NTT*NAA
 1531 continue 
      if(mpin.gt.1)write(22,1711)dst,lx
 1711 format(1x//2x,'Total integrated Rutherford cross section=',1e8.3,
     *' for exp. ',1i2/)	
 530  CONTINUE
      if(ipinf.eq.0)go to 1
	ngpr=0
      do 1249 lx=1,nexpt
	nged=ndst(lx)
	if(irawex(lx).eq.0)nged=nang(lx)
	if(lx.eq.1)go to 1277
      ngpr=ngpr+idr*jpin(lx-1)*ndst(lx-1)
 1277 lpin=jpin(lx)
      if(lpin.eq.0)lpin=1
      do 1251 jgd=1,nged
      do 1250 jd=1,idr
 1250 grad(jd)=0.
      do 1257 mpin=1,lpin
	rewind 15
      ndum=ngpr+(jgd-1)*idr+(mpin-1)*jgd*idr
      if(ndum.eq.0)go to 1260
      do 1261 jd=1,ndum
 1261 read(15,*)xx
 1260 continue
      do 1259 jd=1,idr
      read(15,*)xx
 1259 grad(jd)=grad(jd)+xx
 1257 continue
      write(17,*)(grad(jd),jd=1,idr)
 1251 continue
 1249 continue
      rewind 15
      rewind 17
      do 1252 lx=1,nexpt
      nged=ndst(lx)
      if(irawex(lx).eq.0)nged=nang(lx)
      do 1252 ija0=1,nged
      read(17,*)(grad(jdy),jdy=1,idr)
      do 1253 jd=1,idr
 1253 write(15,*)grad(jd)
 1252 continue
 483   FORMAT(1X,//50X,17HINTEGRATED YIELDS//
     *5X,11HEXPERIMENT ,1I2,2X,9HDETECTOR ,1I2/
     *5X,13HENERGY RANGE ,1F8.3,3H---,1F8.3,1X,3HMEV,
     *3X,23HSCATTERING ANGLE RANGE ,1F7.3,3H---,1F7.3,1X,
     *3HDEG//5X,2HNI,5X,2HNF,5X,2HII,5X,2HIF,5X,
     *5HYIELD,5X,16HNORMALIZED YIELD/)
      GO TO 1
 673  read(jzb,*)IMODE,NPTL,CHIOK,CONU,XTEST,LOCKF,NLOCK,
     *IFBFL,LOCKS,DLOCK
      OP2=OPCJA
      IMIN=IMIN+1
      IF(IMIN.EQ.1)GO TO 805
      GO TO 674
 800  NMAX=0
      IF(IABS(IPRM(1)).EQ.1)WRITE(22,428)
 428  FORMAT(1X/40X,6HLEVELS,//5X,5HINDEX,5X,6HPARITY
     *,9X,4HSPIN,11X,11HENERGY(MEV))
      NDIMA=NDIM+1
      DO 810 K=1,NDIMA
      read(jzb,*)IPO1,IPO2,PO2,PO1
      IF(IPO1.EQ.0)GO TO 900
      IF(IPO1.EQ.1.AND.abs(PO2).lt.1.e-6)ISO=0
      NMAX=NMAX+1
      SPIN(IPO1)=PO2
      IF(K.EQ.1)IPH=IPO2
      IPRC=IPO2-IPH
      IF(IPRC.NE.0)IPRC=1
      IFAC(IPO1)=(-1)**(IPRC-IFIX(PO2-SPIN(1)))
      EN(IPO1)=PO1
      PRP='+'
      IF(IPO2.EQ.-1)PRP='-'
      IF(IABS(IPRM(1)).EQ.1)WRITE(22,995)IPO1,PRP,SPIN(IPO1),EN(IPO1)
 810  CONTINUE
 995  FORMAT(6X,1I2,11X,1A1,10X,1F4.1,8X,1F10.4)
      GO TO 900
 801  DO 820 K=1,NMEMX
      IF(OP2.EQ.'GOSI')GO TO 608
      IOPRI=1
      read(jzb,*)IPO1,IPO2,PO1
      GO TO 607
 608  read(jzb,*)IPO1,IPO2,PO1,BL,BU
      IOPRI=2
      ICG=2
 607  IF(IPO1.EQ.0)GO TO 234
      IF(IPO2.EQ.0)GO TO 821
      MULTI(LA)=MULTI(LA)+1
      INDX=INDX+1
      IF(IPO1.GT.IABS(IPO2))GO TO 822
      IF(IPO1.EQ.IPO3)GO TO 824
      IF(IPO1.LT.IPO3)GO TO 825
      IPO3=IPO1
 824  ELM(INDX)=PO1
      mlt(indx)=la
      LEAD(1,INDX)=IPO1
      LEAD(2,INDX)=IABS(IPO2)
       LDNUM(LA,IPO1)=LDNUM(LA,IPO1)+1
      IF(OP2.NE.'GOSI')GO TO 820
      IF(IPO2.LT.0)GO TO 334
      ELMU(INDX)=BU
      ELML(INDX)=BL
      IF(abs(bl-bu).lt.1.e-6)GO TO 333
      IVAR(INDX)=2
      IF(LA.GT.4)IVAR(INDX)=1
      GO TO 335
 333  IVAR(INDX)=0
      GO TO 335
 334  IVAR(INDX)=10000*IFIX(BL)+IFIX(BU)
 335  ISIP=ISEX(IPO1)+1
      ISEX(IABS(IPO2))=MIN0(ISIP,ISEX(IABS(IPO2)))
      GO TO 820
 821  IF(IPO1.LE.LA)GO TO 823
      LAMMAX=LAMMAX+1
      LAMDA(LAMMAX)=IPO1
      IPO3=0
      IF(INDX.EQ.0)GO TO 235
 234  DO 233 KK=1,INDX
       IF(abs(ELM(KK)).le.1.e-6)ELM(KK)=1.E-6
      IF(IVAR(KK).LT.10000)GO TO 233
      KK1=IVAR(KK)/10000
      KK2=IVAR(KK)-10000*KK1
       LA1=LA
       IF(KK2.LT.100)GO TO 7428
       LA1=KK2/100
       KK2=KK2-100*LA1
 7428 INX1=MEM(KK1,KK2,LA1)
      ELML(KK)=ELML(INX1)*ELM(KK)/ELM(INX1)
      ELMU(KK)=ELMU(INX1)*ELM(KK)/ELM(INX1)
      SA(KK)=ELM(KK)/ELM(INX1)
      ivari(kk)=ivar(kk)
      IVAR(KK)=1000+INX1
      IF(ELMU(KK).GT.ELML(KK))GO TO 233
      ELMI=ELMU(KK)
      ELMU(KK)=ELML(KK)
      ELML(KK)=ELMI
 233  CONTINUE
       IF(IPO1.EQ.0)GO TO 870
 235  LA=IPO1
      IF(LA.GT.LMAXE.AND.LA.LE.6)LMAXE=LA
 820  CONTINUE
 870  MEMAX=INDX
      if(jzb.eq.25)m25=memax
	if(jzb.eq.26)m26=memax
      IF(LA.GT.6)MAGEXC=1
      MEMX4=MULTI(1)+MULTI(2)+MULTI(3)+MULTI(4)
      MEMX6=MEMX4+MULTI(5)+MULTI(6)
      IF(IABS(IPRM(1)).EQ.1)CALL PRELM(IOPRI)
      DO 1213 KH=1,NMAX
      IF(ISEX(KH).EQ.1111)ISEX(KH)=1
 1213 CONTINUE
      DO 2123 KH=1,MEMAX
 2123 IVARH(KH)=IVAR(KH)
      GO TO 900
 803  READ(jzb,907)OP1,FIPO1
      IPO1=IFIX(FIPO1)
      IF(OP1.EQ.'ACP,')ACCA=10.**(-FIPO1)
      IF(OP1.EQ.'SEL,')ITS=2
      IF(OP1.EQ.'SMR,')IOSR=1
      IF(OP1.EQ.'FMI,')IFM=1
      IF(OP1.EQ.'TEN,')ITNO=1
      IF(OP1.EQ.'NCM,')NCM=IPO1
      IF(OP1.EQ.'WRN,')SGW=FIPO1
      IF(OP1.EQ.'INT,')GO TO 1637
      IF(OP1.EQ.'VAC,')GO TO 1327
      IF(OP1.EQ.'DIP,')DIPOL=0.001*FIPO1
      IF(OP1.EQ.'ACC,')ACCUR=10.**(-FIPO1)
      IF(OP1.EQ.'PRT,')GO TO 440
      IF(OP1.EQ.'FIX,')GO TO 1785
      IF(OP1.EQ.'SKP,')GO TO 1329
      IF(OP1.EQ.'CRF,')ICS=1
      IF(OP1.EQ.'LCK,')GO TO 1602
      IF(OP1.EQ.'INR,')INNR=1
      IF(OP1.EQ.'CRD,') GO TO 1673
      IF(OP1.EQ.'CCF,')IPS1=IPO1
      if(OP1.eq.'PIN,')ipine=IPO1
      if(OP1.eq.'PIN,')ipinf=1 
      if(OP1.eq.'PIN,')go to 2347
      IF(OP1.EQ.'END,')GO TO 900
      GO TO 803
 2347 do 2348 ipp=1,ipine
      read(jzb,*)ig1,ig2
 2348 jpin(ig1)=ig2
      go to 803
 1602 read(jzb,*)LCK1,LCK2
      IF(LCK1.EQ.0)GO TO 803
      DO 1603 JJX=LCK1,LCK2
      IVARH(JJX)=0
 1603 IVAR(JJX)=0
      GO TO 1602
 1673 DO 1674 JJX=1,IPO1
      read(jzb,*)IPO2
 1674 IECD(IPO2)=1
      GO TO 803
 1637 DO 1638 JJX=1,IPO1
      read(jzb,*)IPO2,IJX
 1638 INTERV(IPO2)=IJX
      GO TO 803
 1329 DO 1330 JJX=1,IPO1
      read(jzb,*)IJX
 1330 JSKIP(IJX)=0
      GO TO 803
 1327 CONTINUE
      DO 1328 JJX=1,7
      read(jzb,*)IJX,VAL
      IF(IJX.EQ.0)GO TO 803
      G(IJX)=VAL
 1328 CONTINUE
 1785 read(jzb,*)NALLOW
      DO 1786 JJX=1,NALLOW
      read(jzb,*)IJK
      IVAR(IJK)=-IVAR(IJK)
 1786 CONTINUE
      DO 1787 JJX=1,MEMAX
      IF(IVAR(JJX).LT.0)GO TO 1787
      IF(IVAR(JJX).GT.999)GO TO 1787
      IVAR(JJX)=0
 1787 CONTINUE
      DO 1888 JJX=1,MEMAX
      IF(IVAR(JJX).LT.0)IVAR(JJX)=-IVAR(JJX)
      IVARH(JJX)=IVAR(JJX)
 1888 CONTINUE
      GO TO 803
 440  CONTINUE
      DO 444 JJX=1,20
      read(jzb,*)INM1,INM2
      IF(INM1.EQ.0)GO TO 803
 444  IPRM(INM1)=INM2
      GO TO 803
 804  read(jzb,*)NEXPT,IZ,XA
      G(1)=3.
      G(2)=.02
      G(3)=.0345
      G(4)=3.5
      G(5)=FLOAT(IZ)/XA
      G(6)=6.E-06
      G(7)=.6
      DO 840 K=1,NEXPT
      read(jzb,*)IZ1(K),XA1(K),EP(K),TLBDG(K),EMMA(K),MAGA(K),IAXS(K)
     *,FI0,FI1,ISKIN(K),LNORM(K)
      ITTE(K)=0
      IF(XA1(K).LT.0.)ITTE(K)=1
      XA1(K)=ABS(XA1(K))
      FIEX(K,1)=FI0/57.2957795
      FIEX(K,2)=FI1/57.2957795
      IF(TLBDG(K).GE.0.) GO TO 840
      FIEX(K,1)=FIEX(K,1)+3.14159265
      FIEX(K,2)=FIEX(K,2)+3.14159265
 840  CONTINUE
      GO TO 900
 310  CALL READY(IDR,NTAP,0)
      REWIND 3
      REWIND 15
      REWIND 4
 805  CALL CMLAB(0,DSIG,TTTTT)
      IF(ERR)GO TO 9999
      IF(OP2.EQ.'POIN')read(jzb,*)IFWD,SLIM
      IENT=1
      ICG=1
      IF(SPIN(1).lt.1.e-6)ISO=0
      IF(IOBL.GE.1)GO TO 610
      IF(OP2.EQ.'GOSI')GO TO 610
      IAPX=0
       DO 4739 II=1,LP6
 4739 ILE(II)=1
      NCH=0
      DO 901 JEXP=1,NEXPT
      IEXP=JEXP
      TTTTT=TREP(IEXP)
      DSIG=DSIGS(IEXP)
      IF(OP2.EQ.'STAR')GO TO 4740
      JMM=IEXP
       IF(IEXP.EQ.1)GO TO 4740
      DO 4741 LLI=1,LP6
 4741 ILE(LLI)=ILE(LLI)+NYLDE(IEXP-1,LLI)
 4740 CONTINUE
      FI0=FIEX(IEXP,1)
      FI1=FIEX(IEXP,2)
      CALL LOAD(IEXP,1,ICG,0.,JJ)
      CALL ALLOC(ACCUR)
      CALL SNAKE(IEXP,ZPOL)
      CALL SETIN
      DO 799 J=1,LMAX
      POLM=FLOAT(J-1)-SPIN(1)
      CALL LOAD(IEXP,2,ICG,POLM,JJ)
      CALL STING(JJ)
      CALL PATH(JJ)
      CALL INTG(IEXP)
      CALL TENB(J,BTEN,LMAX)
      PR=0.
      IF(OP2.EQ.'STAR'.OR.IPRM(19).EQ.1)WRITE(22,300)(FLOAT(J)-1.
     *-SPIN(1)),IEXP
      DO 788 K=1,ISMAX
      PR=PR+REAL(ARM(K,5))**2+AIMAG(ARM(K,5))**2
 788  IF(OP2.EQ.'STAR'.OR.IPRM(19).EQ.1)WRITE(22,301)
     *IFIX(CAT(K,1)),CAT(K,2),CAT(K,3),REAL(ARM(K,5)),AIMAG(ARM(K,5))
      IF(OP2.EQ.'STAR' .OR.IPRM(19).EQ.1)WRITE(22,302)PR
 799  CONTINUE
      CALL TENS(BTEN)
      IF(ITNO.EQ.0)GO TO 882
      DO 883 K=2,NMAX
      WRITE(17,*)K
      DO 884 KK=1,4
      IN1=(K-1)*28+1+(KK-1)*7
      IN2=IN1+2*KK-2
      WRITE(17,*)(ZETA(KKK),KKK=IN1,IN2)
  884 CONTINUE
  883 CONTINUE
  882 CONTINUE
      SUMMM=0.
      DO 6329 JGL=2,NMAX
      LOCT=(JGL-1)*28+1
 6329 SUMMM=SUMMM+ZETA(LOCT)
      POP1=1.-SUMMM
      JGL=1
      IF(OP2.EQ.'STAR'.OR.IPRM(19).EQ.1)WRITE(22,303)JGL,POP1
      DO 305 JGL=2,NMAX
      LOCT=(JGL-1)*28+1
 305  IF(OP2.EQ.'STAR'.OR.IPRM(19).EQ.1)WRITE(22,303)JGL,ZETA(
     *LOCT)
 303  FORMAT(2X,5HLEVEL,1X,1I2,10X,10HPOPULATION,1X
     *,1E14.6)
      IF(OP2.EQ.'STAR')GO TO 901
      CALL DECAY(CCD,0,CCC)
      NOGELI=NANG(IEXP)
      jgl1=0
      do 7258 js=1,lp2
      do 7258 jgl=1,20
 7258 sumcl(jgl,js)=0.
      do 962 jgl=1,nogeli
      if(irawex(iexp).eq.0)go to 12
      if(op2.eq.'POIN'.and.iprm(20).eq.1)write(23,4478)
     *iexp,jgl,ep(iexp),tlbdg(iexp)
 12   GTH=AGELI(IEXP,JGL,1)
      FIGL=AGELI(IEXP,JGL,2)
      fm=(fi0+fi1)/2.
      CALL ANGULA(YGN,IDR,1,FI0,FI1,TTTTT,GTH,FIGL,JGL,OP2)
      IF(IFMO.EQ.0)GO TO 8232
      ID=ITMA(IEXP,JGL)
      D=ODL(ID)
      RX=D*SIN(GTH)*COS(FIGL-FM)-.25*SIN(TTTTT)*COS(FM)
      RY=D*SIN(GTH)*SIN(FIGL-FM)-.25*SIN(TTTTT)*SIN(FM)
      RZ=D*COS(GTH)-.25*COS(TTTTT)
      RL=SQRT(RX*RX+RY*RY+RZ*RZ)
      THC=TACOS(RZ/RL)
      SF=D*D/RL/RL
      FIC=ATAN2(RY,RX)
      CALL ANGULA(YGP,IDR,1,FI0,FI1,TTTTT,THC,FIC,JGL,OP2)
      DO 8233 IXL=1,IDR
      IXM=KSEQ(IXL,3)
      TFAC=TAU(IXM)
 8233 YGN(IXL)=YGN(IXL)+.01199182*TFAC*BETAR(IEXP)*
     *(SF*YGP(IXL)-YGN(IXL))
 8232 CONTINUE
      if(irawex(iexp).eq.0)go to 7262
      ipd=itma(iexp,jgl)
      do 7263 jyi=1,idr
      ni=kseq(jyi,3)
      nf=kseq(jyi,4)
      decen=en(ni)-en(nf)
      cocos=sin(ttttt)*sin(gth)*cos(fm-figl)+cos(ttttt)*
     *cos(gth)
      decen=decen*(1.+betar(iexp)*cocos)
      call effix(ipd,decen,effi)
      if(op2.eq.'POIN'.and.iprm(20).eq.1)write(23,480)
     *ni,nf,spin(ni),spin(nf),decen,effi
      ygn(jyi)=ygn(jyi)*effi
 7263 continue
      inclus=iclust(iexp,jgl)
      if(inclus.eq.0)go to 7262
      do 7264 jyi=1,idr
 7264 sumcl(inclus,jyi)=sumcl(inclus,jyi)+ygn(jyi)
      if(jgl.ne.lastcl(iexp,inclus))go to 962
      do 7265 jyi=1,idr
 7265 ygn(jyi)=sumcl(inclus,jyi)
 7262 continue
      jgl1=jgl1+1
      lu=ile(jgl1)
      if(op2.eq.'POIN'.or.iprm(11).eq.1)
     *write(22,478)
     *iexp,jgl1,ep(iexp),tlbdg(iexp)
      JMM=0
      DO 961 JYI=1,IDR
      NI=KSEQ(JYI,3)
      NF=KSEQ(JYI,4)
      IF(OP2.EQ.'POIN'.OR.IPRM(11).EQ.1)
     *WRITE(22,480)
     *NI,NF,SPIN(NI),SPIN(NF),YGN(JYI),YGN(JYI)/YGN(IDRN)
      IF(IFWD.NE.1)GO TO 3211
      IF((YGN(JYI)/YGN(IDRN)).LT.SLIM)GO TO 3211
      IF(JGL1.EQ.1)SH1=YGN(IDRN)
      JMM=JMM+1
      CORF(JMM,1)=FLOAT(NI)
      CORF(JMM,2)=FLOAT(NF)
      CORF(JMM,3)=YGN(JYI)
      IF(YGN(JYI).GE.YGN(IDRN))CORF(JMM,4)=CORF(JMM,3)/20.
      IF(YGN(JYI).LT.YGN(IDRN))CORF(JMM,4)=CORF(JMM,3)*
     *(.05+.2*(1.-YGN(JYI)/YGN(IDRN)))
 3211 CONTINUE
      IF(OP2.NE.'CORR')GO TO 961
      READ(15,*)YYDD
      NCH=NCH+1
      JJJJ=IY(LU,JGL1)/1000
      JYI1=IY(LU,JGL1)-JJJJ*1000
      IF(IY(LU,JGL1).NE.JYI.AND.JJJJ.NE.JYI.AND.JYI1.NE.JYI)GO TO 961
      IF(IY(LU,JGL1).LT.1000)GO TO 3743
      JYI2=JYI1-JJJJ
      IF(JYI2.LE.0)GO TO 961
      DO 3744 IHUJ=1,JYI2
 3744 READ(15,*)YYD1
      YYDD=YYDD+YYD1
      YGN(JYI)=YGN(JYI)+YGN(JYI1)
      REWIND 15
      DO 3745 IHUJ=1,NCH
 3745 READ(15,*)YYD1
 3743 IF(IEXP.EQ.1.AND.LU.EQ.NYLDE(1,1).AND.JGL1.EQ.1)
     *CNST=YYDD/YGN(JYI)
      CORF(LU,JGL1)=YEXP(JGL1,LU)
      YEXP(JGL1,LU)=YEXP(JGL1,LU)/YYDD*YGN(JYI)
      DYEX(JGL1,LU)=DYEX(JGL1,LU)/YYDD*YGN(JYI)
      LU=LU+1
 961  CONTINUE
      IF(IFWD.NE.1)GO TO 962
      XW=1.
      WRITE(4,*)IEXP,jgl1,IABS(IZ1(IEXP)),ABS(XA1(IEXP)),
     *ABS(EP(IEXP)),JMM,XW
      DO 3209 JYI=1,JMM
 3209 WRITE(4,*)IFIX(CORF(JYI,1)),IFIX(CORF(JYI,2)),CORF(JYI,3)
     *,CORF(JYI,4)
 962  CONTINUE
      IF(OP2.NE.'CORR')GO TO 901
      jgl1=0
      DO 4093 JGL=1,NOGELI
      if(irawex(jexp).eq.0)go to 4095
      inclus=iclust(jexp,jgl)
      if(inclus.eq.0)go to 4095
      if(jgl.ne.lastcl(jexp,inclus))go to 4093
 4095 continue
      jgl1=jgl1+1
      READ(3,*)NE,NA,ZP,AP,XEP,NVAL,WAGA
      WRITE(4,*)NE,NA,ZP,AP,EP(IEXP),NVAL,WAGA
      WRITE(22,9990) IEXP,JGL1
 9990 FORMAT(///10X,10HEXPERIMENT,1X,I2,8X,8HDETECTOR,1X,I2,
     1       //9X,2HNI,5X,2HNF,5X,4HYEXP,8X,4HYCOR,8X,5HCOR.F/)
       ILE1=ILE(JGL1)
      DO 313 ITP=1,NVAL
      READ(3,*)NS1,NS2,FIEX1(1,1,1),FIEX1(1,1,2)
      LTRN=IY(ILE1+ITP-1,JGL1)
      IF(LTRN.LT.1000)GO TO 3747
      LTRN1=LTRN/1000
      NS1=KSEQ(LTRN1,3)*100
      NS2=KSEQ(LTRN1,4)*100
      LTRN2=LTRN-LTRN1*1000
      NS1=NS1+KSEQ(LTRN2,3)
      NS2=NS2+KSEQ(LTRN2,4)
      GO TO 312
 3747 NS1=KSEQ(LTRN,3)
      NS2=KSEQ(LTRN,4)
 312  YCORR=YEXP(JGL1,ILE1+ITP-1)
      WRITE(4,*)NS1,NS2,YCORR,DYEX(JGL1,ILE1+ITP-1)
 313  WRITE(22,9991) NS1,NS2,CORF(ILE1+ITP-1,JGL1),YCORR,
     1              YCORR/CORF(ILE1+ITP-1,JGL1)
 9991 FORMAT(5X,I4,5X,I4,3X,E8.3,4X,E8.3,4X,E8.3)
 4093 CONTINUE
 901   CONTINUE
      IF(OP2.EQ.'STAR')OPH=OP2
      IF(OP2.EQ.'STAR')GO TO 1
      IF(OP2.NE.'CORR')GO TO 1
      NTAP=4
      CALL READY(IDR,NTAP,IPRI)
      REWIND NTAP
      GO TO 1
 610  CONTINUE
      IF(IOBL.GE.1)GO TO 2001
      REWIND 7
	rewind 27
	mdupa=7
	if(ibps.eq.1)mdupa=27
      DO 2929 IUY=1,6
      READ(mdupa,*)(XIR(IUY,JJ),JJ=1,NEXPT)
2929  READ(mdupa,*)(ZMIR(IUY,1,JJ),ZMIR(IUY,2,JJ),JJ=1,NEXPT)
      DO 2002 JJ=1,NEXPT
      DO 2003 JK=1,4
      DO 2003 KUKU=1,6
 2003 READ(mdupa,*)(PARXM(JJ,JK,JL,KUKU),JL=1,10)
      DO 2017 JK=1,12
 2017 READ(mdupa,*)(PARX(JJ,JK,JL),JL=1,5)
 2002 CONTINUE
      DO 6342 JGS=1,MEMAX
      DO 6342 JGR=1,7
 6342 QAPR(JGS,1,JGR)=0.
      GO TO 2000
 2001 IENT=1
      ICG=2
      NMAXH=NMAX
      LMAX1=LMAX
      SH1=SPIN(1)
      SH2=SPIN(2)
      IH1=IFAC(1)
      IH2=IFAC(2)
      MAGH=MAGEXC
      LMAXH=LMAXE
      ISOH=ISO
      ISO=0
      EH1=ELM(1)
      LH1=LEAD(1,1)
      LH2=LEAD(2,1)
      LAMH=LAMMAX
      MEMH=MEMAX
      DO 699 KH=1,8
      IHLM(KH)=MULTI(KH)
      IHLM(KH+24)=LDNUM(KH,2)
      IHLM(KH+8)=LAMDA(KH)
 699  IHLM(KH+16)=LDNUM(KH,1)
      DO 611 JEXP=1,NEXPT
      IEXP=JEXP
      INTVH=INTERV(IEXP)
      DO 6241 JGS=1,MEMAX
      DO 6241 JGR=1,7
 6241 QAPR(JGS,1,JGR)=0.
      DO 2930 IUY=1,6
2930  XIR(IUY,IEXP)=0.
      EMHL1=EMMA(IEXP)
      EMMA(IEXP)=FLOAT(MAGA(IEXP))
      JDE=2
      IF(MAGA(IEXP).EQ.0)JDE=1
      DO 2931 IUY=1,6
      ZMIR(IUY,1,IEXP)=0.
 2931 ZMIR(IUY,2,IEXP)=0.
      CALL LOAD(IEXP,1,2,0.,JJ)
      DO 631 JGS=1,LMAX
      POLM=FLOAT(JGS-1)-SPIN(1)
      CALL LOAD(IEXP,3,2,POLM,JJ)
      CALL PATH(JJ)
      CALL LOAD(IEXP,2,2,POLM,JJ)
      ICTL=1
      DO 612 KK=1,6
      LL=IHLM(KK)
      IF(LL.EQ.0)GO TO 612
      LFINI=LL+ICTL-1
      ICT=ICTL
      DO 613 LLL=ICT,LFINI
      ICTL=ICTL+1
      IF(JGS.NE.1)GO TO 632
      XIR(KK,IEXP)=AMAX1(XIR(KK,IEXP),ABS(XI(LLL)))
 632  R1=ABS(QAPR(LLL,1,1))
      R2=ABS(QAPR(LLL,1,4))
      R3=ABS(QAPR(LLL,1,7))
      RM=AMAX1(R1,R2,R3)
      BMX=AMAX1(ABS(ELMU(LLL)),ABS(ELML(LLL)))
      ZMIR(KK,2,IEXP)=AMAX1(ZMIR(KK,2,IEXP),RM*BMX/ABS(ELM(LLL)),RM)
      R1=ABS(QAPR(LLL,1,2))
      R2=ABS(QAPR(LLL,1,3))
      R3=ABS(QAPR(LLL,1,5))
      R4=ABS(QAPR(LLL,1,6))
      RM=AMAX1(R1,R2,R3,R4)
      ZMIR(KK,1,IEXP)=AMAX1(ZMIR(KK,1,IEXP),RM*BMX/ABS(ELM(LLL)),RM)
 613  CONTINUE
      IF(ZMIR(KK,1,IEXP).LT..5)ZMIR(KK,1,IEXP)=.5
      IF(ZMIR(KK,2,IEXP).LT..5)ZMIR(KK,2,IEXP)=.5
 612  CONTINUE
 631  CONTINUE
      DO 633 KK=1,6
      XIR(KK,IEXP)=XIR(KK,IEXP)*1.01
      DO 698 KH=1,8
      MULTI(KH)=0
      LAMDA(KH)=0
      LDNUM(KH,2)=0
 698  LDNUM(KH,1)=0
      NMAX=2
      ELM(1)=1.
      LEAD(1,1)=1
      LEAD(2,1)=2
      SPIN(1)=0.
      IFAC(1)=1
      LAMMAX=1
      MEMAX=1
      MAGEXC=0
      KKK=0
      ICG=1
      IF(IHLM(KK).EQ.0)GO TO 633
      MULTI(KK)=1
      LAMDA(1)=KK
      SPIN(2)=FLOAT(KK)
      IFAC(2)=1
      LDNUM(KK,1)=1
      ICG=1
      CALL LOAD(IEXP,1,ICG,0.,JJ)
      CALL LOAD(IEXP,2,ICG,0.,JJ)
      CALL PATH(1)
      SZ1=AMIN1(ZMIR(KK,1,IEXP),10.)
      SZ2=ZMIR(KK,2,IEXP)/50.
      ACOF=2.4009604E-3/ZMIR(KK,2,IEXP)
      BCOF=8.163265E-4
      DO 614 JD=1,JDE
      NKSI=5
      IF(JD.EQ.2)NKSI=10
      IF(MAGA(IEXP).EQ.0)NKSI=10
      DO 625 JK=1,3
 625  ZETA(JK)=0.
      NZ=50
      IF(JD.EQ.1.AND.MAGA(IEXP).NE.0)NZ=1
      DO 606 JK=1,NKSI
      XI(1)=XIR(KK,IEXP)*(JK-1)/(NKSI-1)
      IF(JK.EQ.1)XI(1)=.02
      S11=0.
      S21=0.
      S12=0.
      S22=0.
      PH1=0.
      PH2=0.
      DO 615 JZ=1,NZ
      ZETA(JD)=SZ2*JZ
      IF(JD.EQ.1.AND.MAGA(IEXP).NE.0)ZETA(JD)=SZ1
      IF(ZETA(JD).LT..1)INTERV(IEXP)=1000
      IF(ZETA(JD).GE..1)INTERV(IEXP)=INTVH
      CALL ALLOC(ACCUR)
      CALL SNAKE(IEXP,ZPOL)
      CALL SETIN
      CALL STING(1)
      IF(KK.LE.2)GO TO 1082
      ARM(1,5)=(.9999999,0.)
      ARM(2,5)=(1.2E-6,0.)
      ARM(1,6)=(.9999998,0.)
      ARM(2,6)=(.9E-6,0.)
      DO 1083 KH=1,4
      ARM(1,KH)=(-1.E-6,0.)
      ARM(2,KH)=(1.E-6,0.)
 1083 CONTINUE
 1082 CALL INTG(IEXP)
      JP=2
      IF(MAGA(IEXP).NE.0.AND.JD.EQ.2)JP=3
      P=REAL(ARM(1,5))
      R=AIMAG(ARM(1,5))
      QR=REAL(ARM(JP,5))
      S=AIMAG(ARM(JP,5))
      TEST=P*P+R*R+QR*QR+S*S
      P=P/SQRT(TEST)
      S=ABS(R/S)
      IF(JK.NE.1)GO TO 620
      IF(MAGA(IEXP).EQ.0)GO TO 2418
      IF(JD.NE.2.AND.MAGA(IEXP).NE.0)GO TO 620
 2418 Q1=0.
      GO TO 621
 620  Q1=ARCTG(S,PH1,PI)
      PH1=Q1
 621  IF(JK.NE.1)GO TO 622
      IF(JD.NE.1.OR.MAGA(IEXP).EQ.0)GO TO 622
      Q2=0.
      GO TO 623
 622  Q2=ARCCOS(P,PH2,PI)
      PH2=Q2
 623  CONTINUE
      Q1=Q1/ZETA(JD)/2.
      Q2=Q2/ZETA(JD)
      IF(JD.EQ.1.AND.MAGA(IEXP).NE.0)Q2=-Q2
      IF(JD.EQ.1.AND.MAGA(IEXP).NE.0)GO TO 615
      S11=S11+Q1
      S12=S12+Q1*JZ
      S21=S21+Q2
      S22=S22+JZ*Q2
 615  CONTINUE
      IF(JD.EQ.1.AND.MAGA(IEXP).NE.0)GO TO 616
      PARXM(IEXP,1,JK,KK)=ACOF*(2.*S12-51.*S11)
      PARXM(IEXP,2,JK,KK)=BCOF*(101.*S11-3.*S12)
      PARXM(IEXP,3,JK,KK)=ACOF*(2.*S22-51.*S21)
      PARXM(IEXP,4,JK,KK)=BCOF*(101.*S21-3.*S22)
      GO TO 606
 616  PARX(IEXP,2*KK-1,JK)=Q1
      PARX(IEXP,2*KK,JK)=Q2
 606  CONTINUE
 614  CONTINUE
 633  CONTINUE
      EMMA(IEXP)=EMHL1
      NMAX=NMAXH
      SPIN(1)=SH1
      SPIN(2)=SH2
      IFAC(1)=IH1
      IFAC(2)=IH2
      MAGEXC=MAGH
      ISO=ISOH
      ELM(1)=EH1
      LEAD(1,1)=LH1
      LEAD(2,1)=LH2
      LAMMAX=LAMH
      MEMAX=MEMH
      DO 696 KH=1,8
      LDNUM(KH,2)=IHLM(KH+24)
      MULTI(KH)=IHLM(KH)
      LAMDA(KH)=IHLM(KH+8)
 696  LDNUM(KH,1)=IHLM(KH+16)
      INTERV(IEXP)=INTVH
 611  CONTINUE
      REWIND 7
	rewind 27
	jkmx=7
	if(ibps.eq.1)jkmx=27
      DO 2932 IUY=1,6
      WRITE(jkmx,*)(XIR(IUY,JJ),JJ=1,NEXPT)
 2932 WRITE(jkmx,*)(ZMIR(IUY,1,JJ),ZMIR(IUY,2,JJ),JJ=1,NEXPT)
      DO 4005 JJ=1,NEXPT
      DO 4006 JK=1,4
      DO 4006 KUKU=1,6
 4006 WRITE(jkmx,*)(PARXM(JJ,JK,JL,KUKU),JL=1,10)
      DO 4008 JK=1,12
 4008 WRITE(jkmx,*)(PARX(JJ,JK,JL),JL=1,5)
 4005 CONTINUE
      jkmy=12
	if(ibps.eq.1)jkmy=32
      do 8112 kuku=1,memax
 8112 write(jkmy,*)elm(kuku)	
      DO 1739 JJ=1,2
      DO 1739 JJ1=1,LP1
 1739 IDIVE(JJ1,JJ)=1
 2000 IF(IPRM(12).EQ.0.and.op2.ne.'MAP ')GO TO 19
      IPRM(12)=0
      DO 21 JEX=1,NEXPT
      DO 22 LEX=1,6
      IF(MULTI(LEX).EQ.0)GO TO 22
      WRITE(22,8300)JEX,XIR(LEX,JEX)
 8300 FORMAT(1X//30X,10HEXPERIMENT,1X,1I2,10X,7HMAX.XI=,1F6.4)
      WRITE(22,8301)LEX,ZMIR(LEX,2,JEX)
 8301 FORMAT(1X/30X,1HE,1I1,8X,4HMI=0,5X,9HMAX.ZETA=,1F6.3//)
      WRITE(22,8302)
 8302 FORMAT(5X,2HXI,13X,2HQ1,22X,2HQ2///13X,5HSLOPE,2X,
     *9HINTERCEPT,7X,5HSLOPE,5X,9HINTERCEPT//)
      DO 8401 KEX=1,10
      XXI=XIR(LEX,JEX)*(KEX-1)/9.
 8401 WRITE(22,8303)XXI,(PARXM(JEX,ILX,KEX,LEX),ILX=1,4)
 8303 FORMAT(2X,1F6.4,3X,1E8.2,2X,1E8.2,6X,1E8.2,2X,1E8.2)
      IF(MAGA(JEX).EQ.0)GO TO 22
      WRITE(22,8405)LEX,ZMIR(LEX,1,JEX)
 8405 FORMAT(1X//30X,1HE,1I1,8X,7HMI=+/-1,5X,9HMAX.ZETA=,1F6.3//)
      WRITE(22,8302)
      DO 8406 KEX=1,5
      XXI=XIR(LEX,JEX)*(KEX-1)/4.
      U=0.
 8406 WRITE(22,8303)XXI,U,PARX(JEX,2*LEX-1,KEX),U,PARX(JEX,
     *2*LEX,KEX)
 22   CONTINUE
 21   CONTINUE
 19   CONTINUE
       IF(OP2.NE.'GOSI'.AND.OP2.NE.'ERRO')GO TO 1
       IF(OP2.EQ.'ERRO')GO TO 7
 674  CONTINUE
      mret=iabs(mret)
      DO 1487 KH1=1,MEMAX
 1487 HLM(KH1)=ELM(KH1)
      LFAGG=0
      DO 2127 KH1=1,MEMAX
 2127 IVAR(KH1)=IVARH(KH1)
      if(nawr.eq.1)go to 3074
      mcfix=1
      CALL MINI(CHISQ,1.e+38,NPTL,CONU,IMODE,IDR,XTEST,0,0,0,BTEN)
	mawr=mawr+1
	if(jzb.eq.25)ccch1=chisq
	if(jzb.eq.26)ccch2=chisq
	chir=chir+chisq
	if(jzb.eq.26)cht=abs(chir-chp)
	if(jzb.eq.26.and.cht.lt.0.1)mret=0
	if(jzb.eq.26.and.mawr.ge.20)mret=0
	if(jzb.eq.26)nawr=1
	do 3073 kh1=1,32
	do 3073 kh2=1,lp3
      if(jzb.eq.26)cnor2(kh1,kh2)=cnor(kh1,kh2)
      if(jzb.eq.25)cnor1(kh1,kh2)=cnor(kh1,kh2)	
 3073 continue
      if(jzb.ne.26)go to 3076
       do 3075 kh1=1,32
	do 3075 kh2=1,75
c	cnor1(kh1,kh2)=cnor1(kh1,kh2)*cnor2(kh1,kh2)*2/(cnor1(kh1,kh2)
c     *+cnor2(kh1,kh2))
c      cnor1(kh1,kh2)=(cnor1(kh1,kh2)*ccch1+cnor2(kh1,kh2)*ccch2)/
c     *(ccch1+ccch2)	
      cnor1(kh1,kh2)=(cnor1(kh1,kh2)+cnor2(kh1,kh2))/2.
 3075 continue
 3076 continue
      if(ibps.eq.0.and.jzb.eq.25)jzb=26
      if(ibps.eq.0)go to 1408
	jzb=25
	go to 1408
 3074	mcfix=0
      CALL MINI(CHISQ,CHIOK,NPTL,CONU,IMODE,IDR,XTEST,0,0,0,BTEN)
	if(ibps.eq.0)jzb=25
	if(ibps.eq.0)mrepf=1
	if(ibps.eq.0)go to 3911
	ibps=ibps+1
	if(ibps.eq.1)jzb=26
	if(ibps.eq.1)go to 1408
 3911 IF(IPS1.EQ.0)GO TO 9999
      IMIN=IMIN+1
      DO 1336 IVA=1,LP1
 1336 JSKIP(IVA)=1
      irix=12
	if(ibps.eq.2)irix=32
      REWIND irix
      DO 5789 LKJ=1,MEMAX
 5789 WRITE(irix,*)ELM(LKJ)
      IF(IFM.EQ.1)CALL PRELM(3)
      IF(IFM.EQ.1)GO TO 9999
	go to 430
  950 CONTINUE
      CALL ADHOC(OPH,IDR,NFD,NTAP,IYR)
      GO TO 1
 822  WRITE(22,522)
      GO TO 999
 899   READ(jzb,998)(TITLE(K),K=1,20)
      WRITE(22,429)(TITLE(K),K=1,20)
 429  FORMAT(10X,20A4/10X,100(1H-))
      GO TO 1
 522  FORMAT(5X,48HERROR-M.E. DOES NOT BELONG TO THE UPPER TRIANGLE)
 823  WRITE(22,523)
      GO TO 999
 523  FORMAT(5X,39HERROR-WRONG SEQUENCE OF MULTIPOLARITIES)
 825  WRITE(22,525)
      GO TO 999
 598  WRITE(22,599)
 599  FORMAT(1X///10X,44HERROR-INSUFFICIENT SPACE FOR E-THETA INTEGR
     *,5HATION)
      GO TO 999
 525  FORMAT(5X,38HERROR-REPEATED APPEARANCE OF THE STATE)
 998  FORMAT(20A4)
 997   FORMAT(1A3,1A4)
 907  FORMAT(1A4,1F7.1)
 996   FORMAT(1A4)
 430  CONTINUE
      if(mret.eq.1.and.jzb.eq.26)nawr=0
      if(mret.eq.1)go to 3912
      IF(IPRM(18).NE.0)CALL PTICC(IDR)
      IF(OPH.NE.'GOSI')GO TO 999
      IF(LFAGG.EQ.1)GO TO 999
      IF(IMIN.EQ.0)GO TO 45
      IF(IPRM(4).EQ.-1)GO TO 452
      GO TO 455
 452  CONTINUE
      IPRM(4)=111111
 455  ISKOK=IPRM(7)+IPRM(8)+IPRM(13)+IPRM(14)
      IF(ISKOK.EQ.0.AND.IPRM(4).EQ.111111)GO TO 45
      IF(ISKOK.EQ.0)GO TO 456
      IF(IPRM(7).EQ.1)IPRM(7)=-1
      IF(IPRM(8).EQ.1)IPRM(8)=-1
      IF(IPRM(3).EQ.1.AND.NBRA.NE.0)IPRM(3)=-1
      IF(IPRM(13).EQ.1)IPRM(13)=-1
      IF(IPRM(14).EQ.1)IPRM(14)=-1
 456  CONTINUE
      CALL MINI(CHISQ,CHIOK,+1,CONU,2000,IDR,XTEST,2,0,0,BTEN)
 45   CONTINUE
      CALL MIXR(IVA,1,CHISQ,CHILO)
      IF(IPRM(15).EQ.0.OR.KFERR.EQ.1.OR.IYR.EQ.0)GO TO 468
      WRITE(22,465)
      DO 466 IVA=2,NMAX
      DO 442 IVA1=1,10
      IF(LIFCT(IVA1).EQ.IVA)GO TO 2782
 442  CONTINUE
      GO TO 443
 2782 WRITE(22,2781)IVA,TAU(IVA),TIMEL(1,IVA1),TIMEL(2,IVA1)
      GO TO 4666
 443  WRITE(22,469)IVA,TAU(IVA)
 4666 continue
      if(iva.ne.nmax)go to 466
      if(namx.lt.1)go to 1479
      write(22,1469)
 1469 format(5x,//,'CALCULATED AND EXPERIMENTAL MATRIX ELEMENTS',//)
      write(22,1470)
 1470 format(5x,'NI ','NF ',' EXP. ME   ', 'CURRENT ME','   SIGMA')
      do 1471 kq=1,namx
      ni=iamy(kq,1)
      nf=iamy(kq,2)
      ind=iamx(kq)
      ess=elm(ind)
      esd=eamx(kq,1)
      dsd=eamx(kq,2)
 1471 write(22,1477)ni,nf,esd,ess,(ess-esd)/dsd
 1477 format(5x,1i2,1x,1i2,1x,1f9.4,1x,1f9.4,1x,1f9.4)
 1479 continue
 466  CONTINUE
 468  CONTINUE
      IF(IMIN.NE.0)CALL PRELM(3)
 999  CONTINUE
      IF(ITS.EQ.0)GO TO 9999
      IVA=0
      WRITE(18,*)IVA,IVA,IVA,CHISQ
      IF(ITS.EQ.2)GO TO 9999
      WRITE(15,*)IVA,CHISQ,CHISQ,CHISQ,CHISQ
      CALL KLOPOT(KMAT,RLR)
      GO TO 9999
 9998 ITS=1
      read(jzb,*)KMAT,RLR
      GO TO 1
 9999 WRITE(22,3728)
 3912 continue
      if(mrepf.ne.1)go to 1911
	mrepf=2
	jzb=26
	go to 1408
 3728 FORMAT(15X,'********* END OF EXECUTION **********')
 465  FORMAT(1X//20X,20HCALCULATED LIFETIMES//5X,5HLEVEL
     *,5X,14HLIFETIME(PSEC),5X,3HEXP,8X,5HERROR/)
 2781 FORMAT(7X,1I2,7X,1E10.4,5X,1E10.4,4X,1E10.4)
 469  FORMAT(7X,1I2,7X,1E10.4)
 300  FORMAT(1X//40X,21HEXCITATION AMPLITUDES//10X,2HM=
     *,1F5.1,5X,10HEXPERIMENT,1X,1I2//5X,5HLEVEL,2X,4HSPIN
     *,2X,1HM,5X,14HREAL AMPLITUDE,2X,19HIMAGINARY AMPLITUDE
     *//)
 301  FORMAT(7X,1I2,3X,1F4.1,2X,1F5.1,2X,1E14.6,2X,1E14.6)
 302  FORMAT(1X/5X,21HSUM OF PROBABILITIES=,1E14.6)
 1911 if(mret.eq.1)jzb=25
      if(mret.eq.1)go to 1408
	do 7621 mmmm=1,32
	do 7621 kkkk=1,50
 7621 write(13,*)cnor1(mmmm,kkkk) 	
      STOP
      END
      FUNCTION ARCCOS(A,F,PI)
      Q=TACOS(A)
      QA=Q
       QAP=Q
      IF(Q.GT.F)GO TO 1
      DO 2 J=1,20
      AN=2*J*PI
      DO 2 K=1,2
      QAP=QA
      IE=(-1)**K
      QA=AN+IE*Q
      IF(QA.GT.F)GO TO 1
 2    CONTINUE
 1    ARCCOS=QA
      IF((QA-F).GT.PI/2.)ARCCOS=QAP
      RETURN
      END
      FUNCTION ARCTG(A,F,PI)
      Q=ATAN(A)
      QA=Q
      QAP=Q
      IF(Q.GT.F)GO TO 1
      DO 2 J=1,40
      AN=J*PI
      DO 2 K=1,2
      QAP=QA
      IE=(-1)**K
      QA=AN+IE*Q
      IF(QA.GT.F)GO TO 1
 2    CONTINUE
 1    ARCTG=QA
      IF((QA-F).GT.PI/4.)ARCTG=QAP
      RETURN
      END
      SUBROUTINE LOAD(IEXP,IENT,ICG,POLM,JOJ)
      LOGICAL ERR
      COMMON /CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON /COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/PSPIN/ISHA
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON /PCOM/PSI(500)
      COMMON /CCOUP/ZETA(50000),LZETA(8)
      COMMON /CLM/LMAX
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/CLCOM9/ERR
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/CXI/XI(500)
      COMMON/CAUX0/NCM,EMMA(75)
      COMMON/CAUX4/DISTA
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      COMMON/PTH/IPATH(75),MAGA(75)
      DIMENSION ETAN(75),CPSI(8)
      LMAX=IFIX(SPIN(1)+1.1)
      IF(IENT.NE.1)GO TO 271
      ISHA=0
      ISPI=IFIX(SPIN(1)+.51)
      ISPO=IFIX(SPIN(1)+.49)
      IF(ISPI.NE.ISPO)ISHA=1
      Z1=FLOAT(IABS(IZ1(IEXP)))
      Z2=FLOAT(IZ)
      A1=XA1(IEXP)
      A2=XA
      ZPOL=DIPOL*EP(IEXP)*A2/(Z2*Z2*(1.+A1/A2))
      IF(IZ1(IEXP).LT.0)ZPOL=DIPOL*EP(IEXP)
     **A1/(Z1*Z1*(1.+A2/A1))
      IF(IZ1(IEXP).GT.0)GO TO 11
      AH=A1
      A1=A2
      A2=AH
 11   ETA=Z1*Z2*SQRT(A1/EP(IEXP))/6.349770
      DO 110 M=1,NMAX
      DEP=(1.0+A1/A2)*EN(M)
      ZET=DEP/EP(IEXP)
      SZET=SQRT(1.0-ZET)
 110  ETAN(M)=ETA/SZET
      DO 150 N=1,MEMAX
      I1=LEAD(1,N)
      I2=LEAD(2,N)
 150  XI(N)=ETAN(I1)-ETAN(I2)
      AAZZ=1./(1.+A1/A2)/Z1/Z2
      CPSI(1)=5.169286*AAZZ
      IF(LMAXE.EQ.1)GO TO 200
      AAZ2=AAZZ*AAZZ
      CPSI(2)=14.359366*AAZ2
      IF(LMAXE.EQ.2)GO TO 200
      AAZ3=AAZZ*AAZ2
      CPSI(3)=56.982577 *AAZ3
      IF(LMAXE.EQ.3)GO TO 200
      AAZZ=AAZ2*AAZ2
      CPSI(4)=263.812653*AAZZ
      IF(LMAXE.EQ.4)GO TO 200
      AAZ2=AAZ3*AAZ2
      CPSI(5)=1332.409500*AAZ2
      IF(LMAXE.EQ.5)GO TO 200
      AAZZ=AAZ3*AAZ3
      CPSI(6)=7117.691577*AAZZ
 200   IF(MAGEXC.EQ.0)GO TO 210
      AAZZ=VINF(IEXP)/95.0981942
      CPSI(7)=AAZZ*CPSI(1)
      IF(LAMMAX.EQ.8)GO TO 210
      CPSI(8)=AAZZ*CPSI(2)
 210   ZSQA=Z1*SQRT(A1)
      I3=1
      PPP=1.+A1/A2
      DO 270 I1=1,LAMMAX
      LAM=LAMDA(I1)
      LAM1=LAM
      IF(LAM.GT.6)LAM1=LAM-6
      DO 280 N2=1,NMAX
      NN=LDNUM(LAM,N2)
      IF(NN.EQ.0)GO TO 280
      N3=LEAD(1,I3)
      PP1=EP(IEXP)-PPP*EN(N3)
      DO 260 M1=1,NN
      M2=LEAD(2,I3)
      I2=I3
      I3=I3+1
      PP2=EP(IEXP)-PPP*EN(M2)
 260  PSI(I2)=CPSI(LAM)*ZSQA*(PP1*PP2)**((2.*FLOAT(LAM1)-1.)/4.)
 280  CONTINUE
 270  CONTINUE
      IF(IENT.EQ.1)RETURN
 271  CONTINUE
      DO 803 N=1,NMAX
      NSTART(N)=0
  803 NSTOP(N)=0
      IS=1
      NSTART(1)=1
      DO 380 N=1,NMAX
      WRT=POLM-EMMA(IEXP)
      WRTM=POLM+EMMA(IEXP)
      IF(ICG.EQ.2)WRT=POLM-FLOAT(MAGA(IEXP))
      IF(ICG.EQ.2)WRTM=POLM+FLOAT(MAGA(IEXP))
      IF(WRTM.LT.-SPIN(N))GO TO 371
      IF(ABS(WRT).GT.SPIN(N))WRT=-SPIN(N)
      IF(WRTM.GT.SPIN(N))WRTM=SPIN(N)
      MSTOP=IFIX(WRTM-WRT+1.01)
      DO 370 I=1,MSTOP
      CAT(IS,1)=N
      CAT(IS,2)=SPIN(N)
      CAT(IS,3)=WRT+FLOAT(I-1)
      IF(N.EQ.1.AND.abs(CAT(IS,3)-POLM).lt.1.e-6)JOJ=IS
      IS=IS+1
 370  CONTINUE
      GO TO 372
  371 CONTINUE
      NSTART(N)=0
  372 CONTINUE
      NSTART(N+1)=IS
      NSTOP(N)=IS-1
 380  CONTINUE
      ISMAX=IS -1
      IF(ISMAX.LE.LP10)GOTO 400
      WRITE(22,916) LP10
      GOTO 850
 400  CONTINUE
      IF(IENT.EQ.3)RETURN
      NZ=0
      DO 2 JJ=1,7
      DO 2 JJJ=1,MEMAX
      QAPR(JJJ,1,JJ)=0.
  2   QAPR(JJJ,2,JJ)=0.
      DO 717 I=1,8
 717  LZETA(I)=0
      DO 700 I1=1,LAMMAX
      LAM=LAMDA(I1)
      IF(ICG.EQ.2.AND.LAM.GT.6)GO TO 700
      LA=LAM
      IF(LAM.GT.6)LAM=LAM-6
      RLAM=FLOAT(LAM)
      SSQRT=SQRT(2.*RLAM+1.)
      LZETA(LA)=NZ
      IR=0
 440  IR=IR+1
      IF(IR.GT.ISMAX)GO TO 700
      N=CAT(IR,1)
      IF(ICG.EQ.1)GO TO 1
      IF(MAGA(IEXP).EQ.0.AND.IR.NE.IPATH(N))GO TO 440
      IF(IABS(IR-IPATH(N)).GT.1)GO TO 440
 1    LD=LDNUM(LA,N)
      IF(LD.EQ.0)GO TO 441
      CALL LSLOOP(IR,N,NZ,LD,LAM,LA,SSQRT,ICG,IEXP)
      GO TO 440
 441  IR=IR+NSTOP(N)-NSTART(N)
      GO TO 440
 700   CONTINUE
      IF ( NZ.GT.LP7 ) GOTO 540
      RETURN
 540  WRITE(22,936) LP7
 850  ERR=.TRUE.
      RETURN
 916  FORMAT(27H ERROR-ISMAX EXCEEDS MAGMAX,5X,9H(MAGMAX =,I4,1H))
 936  FORMAT(1x,'ERROR - NUMBER OF ELEMENTS IN ZETA ARRAY EXCEEDS'
     *,5HZEMAX,5X,8H(ZEMAX =,I6,1H))
      END
      SUBROUTINE LSLOOP(IR,N,NZ,LD,LAM,LA,SSQRT,ICG,IEXP)
      COMMON /CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON /COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON /PCOM/PSI(500)
      COMMON /CCOUP/ZETA(50000),LZETA(8)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/CLCOM0/IFAC(75)
      LAM2=2*LAM
      INR=CAT(IR,2)*2.
      RMIR=CAT(IR,3)
      JRMIR=2.*RMIR
      DO 300 I2=1,LD
      M=LEADF(N,I2,LA)
      INDX=MEM(N,M,LA)
      IAPR(INDX,1)=N
      IAPR(INDX,2)=M
      ISMIN=0
      INS=SPIN(M)*2.
      IS1=NSTART(M)
       IF(IS1.EQ.0)GO TO 300
      ISPLUS=IFIX(RMIR-CAT(IS1,3))-LAM
      IF(ISPLUS.GE.0)GOTO 100
      ISMIN=ISPLUS
      ISPLUS=0
 100  IS2=IS1+ISPLUS-1
      MRANGE=2*LAM+1+ISMIN
      IF(IS2+MRANGE.GT.NSTOP(M)) MRANGE=NSTOP(M)-IS2
      IF(MRANGE.LE.0)GOTO 300
      DO 200 I3=1,MRANGE
      IS=IS2+I3
      RMIS=CAT(IS,3)
      IF(ISO.EQ.0.AND.RMIS.GT..1.AND.RMIR.GT..1)GO TO 200
      JG1=-RMIS*2.
      JG2=(RMIS-RMIR)*2.
      IF(ICG.EQ.2.AND.IABS(JG2).GT.2*MAGA(IEXP))GO TO 200
      IF(LA.GT.6.AND.JG2.EQ.0)GO TO 200
      NZ=NZ+1
      IF(NZ.GT.LP7) GOTO 200
      IIEX=(INS+JG1)/2
      PHZ=(-1.0)**IIEX
      ZETA(NZ)=PHZ*PSI(INDX)*SSQRT*WTHREJ(INS,LAM2,INR,JG1,JG2,JRMIR)
      IF(ICG.EQ.1)GO TO 200
      MT=CAT(IS,1)
      CALL CODE7(IR,IS,N,MT,INQA,INDX)
      IF(abs(ELM(INDX)).lt.1.e-6)ELM(INDX)=1.E-6
      IF(INQA.EQ.-1)GO TO 200
      QAPR(INDX,1,INQA)=ZETA(NZ)*ELM(INDX)
      IF(ISO.EQ.0.AND.INQA.EQ.1)QAPR(INDX,1,7)=QAPR(INDX,1,1)*IFAC(M)
 200  CONTINUE
 300  CONTINUE
      RETURN
      END
      FUNCTION LEADF(N1,N2,N3)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      LSUM=0
      N3M=N3-1
      IF(N3M.EQ.0)GO TO 5
      DO 1 K=1,N3M
 1    LSUM=LSUM+MULTI(K)
 5    N1M=N1-1
      IF(N1M.EQ.0)GO TO 6
      DO 2 K=1,N1M
 2    LSUM=LSUM+LDNUM(N3,K)
 6    N1M=LSUM+N2
      LEADF=LEAD(2,N1M)
      RETURN
      END
      FUNCTION MEM(N1,N2,N3)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      MSUM=0
      IF(N3.EQ.1)GO TO 5
      N3M=N3-1
      DO 1 K=1,N3M
 1    MSUM=MSUM+MULTI(K)
 5    N1M=N1-1
      IF(N1M.EQ.0)GO TO 6
      DO 2 K=1,N1M
 2    MSUM=MSUM+LDNUM(N3,K)
 6    N1M=MSUM+1
      N3M=N1M+LDNUM(N3,N1)
      DO 3 K=N1M,N3M
      MSUM=MSUM+1
 3    IF(LEAD(2,K).EQ.N2)GO TO 7
 7    MEM=MSUM
      RETURN
      END
      SUBROUTINE CMLAB(II,DSIG,TETRN)
      LOGICAL ERR
      COMMON/CLCOM9/ERR
      COMMON/SECK/ISKIN(50)
      COMMON/PRT/IPRM(20)
      COMMON/TCM/TETACM(50),TREP(50),DSIGS(50)
      COMMON/BREC/BETAR(50)
      COMMON/CAUX0/NCM,EMMA(75)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      LEXP0=1
      LEXP1=NEXPT
      IF(II.NE.0)LEXP0=II
      IF(II.NE.0)LEXP1=II
      DO 1 LEXP=LEXP0,LEXP1
      IFLAA=0
      IF(TLBDG(LEXP).LT.0)IFLAA=1
      IF(IPRM(1).NE.1)GO TO 1234
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,906)LEXP
 1234 TLBDG(LEXP)=ABS(TLBDG(LEXP))
      A1=XA1(LEXP)
      IF(IZ1(LEXP).LT.0)A1=XA
      A2=XA
      IF(IZ1(LEXP).LT.0)A2=XA1(LEXP)
      Z1=FLOAT(IABS(IZ1(LEXP)))
      Z2=FLOAT(IZ)
      IF(IPRM(1).NE.1)GO TO 2345
      IF(IZ1(LEXP).LT.0.AND.(II.EQ.0.AND.IPRM(10).EQ.1))WRITE(22,907)
     *IZ,XA,IABS(IZ1(LEXP)),XA1(LEXP)
      IF(IZ1(LEXP).GT.0.AND.(II.EQ.0.AND.IPRM(10).EQ.1))WRITE(22,908)
     *IZ,XA,IZ1(LEXP),XA1(LEXP)
 2345 DISTS=1.44*(A1+A2)*Z1*Z2/((A1**.33333+A2**.33333)*1.25+5.)/A2
      DISTA=0.0719949*(1.0+A1/A2)*Z1*Z2/EP(LEXP)
      D2A=20.0*DISTA
      VINF(LEXP)=0.0463365*SQRT(EP(LEXP)/A1)
      IF(IPRM(1).NE.1)GO TO 1111
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,909)EP(LEXP),VINF(LEXP)
      IF(EP(LEXP).GT.DISTS.AND.(II.EQ.0.AND.IPRM(10).EQ.1))
     *WRITE(22,905)(EP(LEXP)/DISTS-1.)*100.
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,910)D2A
 1111 TLBRAD=TLBDG(LEXP)/57.2957795
      ARED=1.0+A1/A2
      EMAX=EP(LEXP)/ARED
      DO 200 N=1,NMAX
 200  IF (EN(N).GT.EMAX) GOTO 100
      GO TO 150
 100  CONTINUE
      WRITE(22,900)EMAX,LEXP
      GOTO 850
 150  EPMIN=EP(LEXP)- EN(NCM)*ARED
      TAUP=SQRT(EP(LEXP)/EPMIN)
      TAU=TAUP*A1/A2
      IF (TAU.LE.1.0) GOTO 120
      TMXDG=TASIN(1.0/TAU)*57.2957795
      IF (TMXDG.GE.TLBDG(LEXP)) GOTO 120
      WRITE(22,902)TMXDG,LEXP
      GOTO 850
 120  TCMRAD=TLBRAD+ TASIN(TAU*SIN(TLBRAD))
      TCMDG=TCMRAD*57.2957795
      IF (TAU.LE.1.0) GOTO 140
      IF(IPRM(1).NE.1)GO TO 4444
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,904) TCMDG,LEXP
 4444 CONTINUE
      IF(ISKIN(LEXP).EQ.1)GO TO 140
      TCMDG=180.+2.*TLBDG(LEXP)-TCMDG
      TCMRAD=TCMDG/57.2957795
 140  CONTINUE
      EPS(LEXP)=1./SIN(TCMRAD/2.)
       TETACM(LEXP)=TCMRAD
      IF(IPRM(1).NE.1)GO TO 5555
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,911)TCMDG,EPS(LEXP)
 5555 CONTINUE
      IF(IZ1(LEXP).GT.0)BETAR(LEXP)=A1*A2/(A1+A2)
     ***2*(1.+TAUP*TAUP-2.*TAUP*COS(TCMRAD))*EPMIN
      IF(IZ1(LEXP).LT.0)BETAR(LEXP)=(A2/(A1+A2))
     ***2*(1.+TAU*TAU+2.*TAU*COS(TCMRAD))*EPMIN
      IF(IPRM(1).NE.1)GO TO 6666
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,920)BETAR(LEXP)
 6666 BETAR(LEXP)=.0463365*SQRT(BETAR(LEXP)/XA)
      IF(IPRM(1).NE.1)GO TO 7777
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,921)BETAR(LEXP)
      IF(II.EQ.0.AND.IPRM(10).EQ.1)WRITE(22,931)
     > EP(LEXP)/(DISTS*.5*(1.+EPS(LEXP)))
 7777 CONTINUE
      IF(IFLAA.EQ.1)GO TO 505
      IF(ABS(TCMDG-180.).LT.1.E-5)GO TO 500
      R3=SIN(TLBRAD)/SIN(TCMRAD)
      R3=R3*R3*ABS(COS(TCMRAD-TLBRAD))
      R3=1./R3
      GO TO 505
 500  R3=(1.-TAU)**2
 505  CONTINUE
      ZCMDG=180.-TCMDG
      ZCMRAD=ZCMDG/57.2957795
      ZLBRAD=ATAN(SIN(ZCMRAD)/(COS(ZCMRAD)+TAUP))
      IF(IFLAA.EQ.0)GO TO 515
      IF(ABS(TCMDG-180.).LT.1.E-5)GO TO 511
      R3=SIN(ZLBRAD)/SIN(ZCMRAD)
      R3=R3*R3
      R3=R3*ABS(COS(ZCMRAD-ZLBRAD))
      R3=1./R3
      TLBDG(LEXP)=ZLBRAD*57.2955795
      GO TO 515
 511  R3=(1.+TAUP)**2
      TLBDG(LEXP)=0.
 515  DSIG=250.*R3*SQRT(EP(LEXP)/(EP(LEXP)-ARED*EN(NCM)))*
     *DISTA*DISTA*(EPS(LEXP))**4
      EROOT(LEXP)=SQRT(EPS(LEXP)*EPS(LEXP)-1.)
      DSIGS(LEXP)=DSIG
      TETRN=ZLBRAD
      IF(IZ1(LEXP).LT.0.)TETRN=TLBRAD
  1   TREP(LEXP)=TETRN
      IPRM(10)=0
      RETURN
 850  ERR=.TRUE.
      RETURN
 900  FORMAT(1X,36HERROR- MAXIMUM EXCITATION ENERGY IS ,F8.4,4H MEV,
     *16H FOR EXPERIMENT ,1I2)
 902  FORMAT(1X,35HERROR- MAXIMUM SCATTERING ANGLE IS ,F7.2,8H DEGREES,
     *16H FOR EXPERIMENT ,1I2)
 904  FORMAT(5X,38HSECOND POSSIBLE CM SCATTERING ANGLE IS ,F7.2,
     *24H DEGREES FOR EXPERIMENT ,1I2)
 905  FORMAT(5X,6H***** ,
     *20HBE CAREFUL-ACCORDING,
     *29H TO D.CLINE BOMBARDING ENERGY,1X,1F6.2,1X,2HPC,1X,
     *39H TOO HIGH FOR HEAD-ON COLLISIONS! *****)
 906  FORMAT(1X,///10X,13H** EXPERIMENT,1X,1I2,1X,2H**//)
 907  FORMAT(5X,25HPROJECTILE EXCITATION OF(,1I3,1H,,1F7.3,
     *5H) ON(,1I3,1H,,1F7.3,1H))
 908  FORMAT(5X,21HTARGET EXCITATION OF(,1I3,1H,,1F7.3,
     *5H) BY(,1I3,1H,,1F7.3,1H))
 909  FORMAT(5X,6HENERGY,1X,1F10.3,1X,3HMEV,5X,4HBETA,1X,
     *1E14.6)
 910  FORMAT(5X,51HDISTANCE OF CLOSEST APPROACH FOR HEAD-ON COLLISIONS
     *,1X,1F10.4,1X,2HFM)
 931  FORMAT(5X,18HBOMBARDING ENERGY=,1F10.3,1X,
     *39HOF SAFE BOMBARDING ENERGY AT THIS ANGLE)
 911  FORMAT(5X,19HCM SCATTERING ANGLE,1X,1F10.3
     *,1X,3HDEG,5X,7HEPSILON,1X,1F10.4)
 920  FORMAT(5X,18HRECOIL ENERGY(MEV),2X,1F10.4)
 921   FORMAT(5X,11HRECOIL BETA,2X,1E14.6)
      END
       SUBROUTINE QE(C,D,B2,C2,D2,B4,B6,D3,B8,C4,D4,B10,D5,
     *B12,D6,LMDA,POL,CQ)
      DIMENSION CQ(7)
      GO TO(110,120,130,140,150,160)LMDA
 110    CQ(1)=0.5*C/B2
      CQ(2)=-0.35355339*D/B2
      RETURN
 120    CQ(1)=0.75*(2.0*C2-D2)/B4*POL
      CQ(2)=-1.83711730*C*D/B4*POL
      CQ(3)=-0.91855865*D2/B4*POL
      RETURN
 130    CQ(1)=1.875*C*(2.0*C2-3.0*D2)/B6
      CQ(2)=-1.62379763*(4.0*C2-D2)*D/B6
      CQ(3)=-5.13489890*C*D2/B6
      CQ(4)=2.09631373*D3/B6
      RETURN
 140    CQ(1)=1.09375000*(8.0*C4-24.0*C2*D2+3.0*D4)/B8
            CQ(2)=-4.89139867*C*(4.0*C2-3.0*D2)*D/B8
      CQ(3)=-3.45874113*(6.0*C2-D2)*D2/B8
      CQ(4)=12.9414244 *C*D3/B8
      CQ(5)=4.57548440*D4/B8
      RETURN
 150    CQ(1)=1.230468*C*(-14.*C2*(9.*D2+B2)+30.*B4)/B10
      CQ(2)=-1.347911*D*(35.*C2*(-3.*D2+B2)+5.*B4)/B10
      CQ(3)=-35.662372*D2*C*(-3.*D2+2.*B2)/B10
      CQ(4)=7.279552*D3*(9.*C2-B2)/B10
      CQ(5)=30.884521*D4*C/B10
      CQ(6)=-9.766543*D5/B10
      RETURN
 160    CQ(1)=2.707031*(21.*C2*(-C2*(11.*D2+4.*B2)+5.*B4)-
     1 5.*B6)/B12
      CQ(2)=-17.543567*D*C*(3.*C2*(-11.*D2+B2)+5.*B4)/
     1 B12
      CQ(3)=-13.869408*D2*(3.*C2*(-11.*D2+5.*B2)+B4)/B12
      CQ(4)=-27.738815*D3*C*(-11.*D2+8.*B2)/B12
      CQ(5)=15.193177*D4*(11.*C2-B2)/B12
      CQ(6)=-71.262308*D5*C/B12
      CQ(7)=-20.571656*D6/B12
      RETURN
      END
      SUBROUTINE QM(C,D,B2,B4,ERT,LMDA,CQ)
      DIMENSION CQ(7)
      IF(LMDA.EQ.8)GO TO 220
      CQ(1)=-.3535533905*ERT/B2
      RETURN
 220  CQ(1)=-.9185586536*C*ERT/B4
      CQ(2)=-CQ(1)*D/C
      RETURN
      END
      SUBROUTINE SNAKE(NEXP,ZPOL)
      DIMENSION LLOC(8),CQ(7),IRL(8)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/ALLC/LOCQ(8,7)
      COMMON/HIPER/SH(365),CH(365)
      ICNT=0
  1    ICNT=ICNT+1
      CALL QRANGE(ICNT,NLM,LLOC,IBM,ICM,IDM,IRL)
      IF(NLM.EQ.0)RETURN
      CHI=CH(ICNT)
      SHI=SH(ICNT)
      B2=EPS(NEXP)*CHI+1.
      POL=1.-ZPOL/B2
      B2=B2*B2
      IF(IBM.EQ.2)GO TO 3
      B4=B2*B2
      IF(IBM.EQ.4)GO TO 3
      B6=B4*B2
      IF(IBM.EQ.6)GO TO 3
      B8=B4*B4
      IF(IBM.EQ.8)GO TO 3
      B10=B6*B4
      IF(IBM.EQ.10)GO TO 3
      B12=B6*B6
 3    IF(ICM.EQ.0)GO TO 4
      C=CHI+EPS(NEXP)
      IF(ICM.EQ.1)GO TO 4
      C2=C*C
      IF(ICM.EQ.2)GO TO 4
      C4=C2*C2
      IF(ICM.EQ.4)GO TO 4
      C6=C2*C4
 4    IF(IDM.EQ.0)GO TO 5
      D=EROOT(NEXP)*SHI
      IF(IDM.EQ.1)GO TO 5
      D2=D*D
      IF(IDM.EQ.2)GO TO 5
      D3=D*D2
      IF(IDM.EQ.3)GO TO 5
      D4=D2*D2
      IF(IDM.EQ.4)GO TO 5
      D5=D3*D2
      IF(IDM.EQ.5)GO TO 5
      D6=D3*D3
 5    DO 7 J=1,NLM
      LMDA=LLOC(J)
      IF(LMDA.GT.6)GO TO 11
      CALL QE(C,D,B2,C2,D2,B4,B6,D3,B8,C4,D4,B10,D5,B12,D6,LMDA,POL,CQ)
      MIMX=LMDA+1
      DO 8 K=1,MIMX
      NIND=LOCQ(LMDA,K)+ICNT
 8    ZETA(NIND+LP7)=CQ(K)
      GO TO 7
 11   LMD=LMDA
      LMDA=LMDA-6
      ERT=EROOT(NEXP)
      CALL QM(C,D,B2,B4,ERT,LMDA,CQ)
      MIMX=LMDA
      DO 9 K=1,MIMX
      NIND=LOCQ(LMD,K)+ICNT
 9    ZETA(NIND+LP7)=CQ(K)
 7    CONTINUE
      GO TO 1
      END
      SUBROUTINE FHIP
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/HIPER/SH(365),CH(365)
      W=-.03
      DO 1 J=1,LP12
      W=W+.03
      EX=EXP(W)
      ER=1./EX
      SH(J)=(EX-ER)/2.
 1    CH(J)=(EX+ER)/2.
      RETURN
      END
      SUBROUTINE ALLOC(ACCUR)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/ALLC/LOCQ(8,7)
      COMMON/RNG/IRA(8),MAXLA
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      CALL RANGEL(ACCUR)
      LOAD=0
      IFLAG=0
      DO 10 J=1,8
      DO 10 K=1,7
 10   LOCQ(J,K)=0
      DO 1 K=1,6
      K1=K+1
      DO 1 J=1,K1
      LOCQ(K,J)=LOAD
 1    LOAD=LOAD+IRA(K)
      DO 2 K=7,8
      K1=K-6
      DO 2 J=1,K1
      LOCQ(K,J)=LOAD
 2    LOAD=LOAD+IRA(K)
      IF(LOAD.LE.LP14)RETURN
      WRITE(22,12)
 12   FORMAT(5X,35HNO SPACE FOR Q FUNCTIONS TABULATION//
     *5X,38HSORRY,JOB WILL BE BRUTALLY TERMINATED!)
      V=-1.
      U=ALOG10(V)
      u=sin(u)
      RETURN
      END
      SUBROUTINE RANGEL(ACC1)
      COMMON/A50/ACC50
      COMMON/RNG/IRA(8),MAXLA
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      ACL=-ALOG(ACC1)
      ACC50=ACC1/50.
      DO 110 I=1,8
      IF(MULTI(I).NE.0)GO TO 100
      IRA(I)=0
      GO TO 110
 100  GO TO(1,2,3,4,5,6,2,3)I
 1    W=ACL-.693
      GO TO 120
 2    W=ACL/2.+.203
      GO TO 120
 3    W=ACL/3.+.536
      GO TO 120
 4    W=ACL/4.+.716
      GO TO 120
 5    W=ACL/5.+.829
      GO TO 120
 6    W=ACL/6.+.962
 120  W=W/.03
      IRA(I)=IFIX(W+1.5)
 110  CONTINUE
      IF(IRA(7).NE.0)IRA(7)=IRA(7)+1
      IF(IRA(8).NE.0)IRA(8)=IRA(8)+1
      RETURN
      END
      SUBROUTINE QRANGE(ICNT,NLM,LLOC,IBM,ICM,IDM,IRL)
      DIMENSION LLOC(8),IRL(8)
      COMMON/RNG/IRA(8),MAXLA
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      IF(ICNT.EQ.1)GO TO 1
      IF(IRL(NLM).GE.ICNT)RETURN
      LD=LLOC(NLM)
      LLOC(NLM)=0
      NLM=NLM-1
      IF(NLM.EQ.0)RETURN
      IF(LD.GT.6)RETURN
      GO TO 10
 1    NLM=0
      DO 2 L=1,8
      LLOC(L)=0
 2    IRL(L)=0
      DO 3 K=1,6
      KE=7-K
      KM=13-K
      IF(KM.GT.8)GO TO 4
      IF(MULTI(KM).EQ.0)GO TO 4
      NLM=NLM+1
      LLOC(NLM)=KM
      IRL(NLM)=IRA(KM)
 4    IF(MULTI(KE).EQ.0)GO TO 3
      NLM=NLM+1
      LLOC(NLM)=KE
      IRL(NLM)=IRA(KE)
 3    CONTINUE
      NLEND=IFIX((FLOAT(NLM)+1.1)/2.)
      DO 5 K=1,NLEND
      KE=NLM-K+1
      LS=LLOC(KE)
      IS=IRL(KE)
      LLOC(KE)=LLOC(K)
      IRL(KE)=IRL(K)
      LLOC(K)=LS
 5    IRL(K)=IS
      L=0
      DO 6 K=1,6
      IF(MULTI(K).EQ.0)GO TO 6
      L=K
 6    CONTINUE
      ICM=MIN0(4,L)
      IBM=2*L
      IDM=L
      L=0
      DO 7 K=7,8
      KE=K-6
      IF(MULTI(K).EQ.0)GO TO 7
      L=KE
 7    CONTINUE
      IBM=MAX0(IBM,2*L)
      IDM=MAX0(IDM,L)
      IF(ICM.EQ.1.AND.L.GT.1)ICM=2
      MAXLA=LLOC(1)
      RETURN
 10   L=LLOC(NLM)
      IF(L.GT.6)L=L-6
      ICM=MIN0(2,L)
      IBM=2*L
      IDM=L
      RETURN
      END
      SUBROUTINE AMPDER(I57)
      COMPLEX ARM,EXPO
      COMMON /COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/AZ/ARM(600,7)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON /CAUX/NPT,NDIV,KDIV,LAMR(8),ISG,D2W,NSW,ISG1
      COMMON/PINT/ISSTAR(76),ISSTO(75),MSTORE(2,75)
      COMMON/ADBXI/EXPO(500)
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      DO 100 K=1,ISMAX
      ARM(K,6)=(0.,0.)
 100  ARM(K,4)=(0.,0.)
      ISG1=ISG
      IF(NPT.EQ.1)ISG1=IABS(ISG1)
      RSG=FLOAT(ISG)
      DO 1 I1=1,LAMMAX
      LAM=LAMDA(I1)
      LAX=LAM
      NZ=LZETA(LAM)
      IF(LAMR(LAM).EQ.0)GO TO 1
      IFLG=1
      NHOLD=1
      GO TO 400
 401  NHOLD=NHOLD+1
      IF(NSTART(NHOLD).EQ.0)GO TO 401
 400  CALL NEWLV(NHOLD,LD,LAM)
      IF(LD.EQ.0)GO TO 401
      IR=NSTART(NHOLD)-1
 200  IR=IR+1
      IF(IR.GT.ISMAX)GO TO 1
      N=CAT(IR,1)
      IF(N.NE.NHOLD)GO TO 203
 204  CALL LAISUM(IR,N,RSG,LAX,LD,NZ,I57)
      GO TO 200
 201  IR=IR+NSTOP(N)-NSTART(N)+1
      N=N+1
      IF(N.GT.NMAX)GO TO 1
      GO TO 202
 203  DO 300 MM=1,LD
      M=MSTORE(1,MM)
      IF(M.EQ.NHOLD)GO TO 300
      INDX=MSTORE(2,MM)
      IBG=ISSTAR(MM)
      IEND=ISSTO(MM)
      DO 301 IS2=IBG,IEND
      ARM(IS2,4)=ARM(IS2,4)+ARM(IS2,6)*ELM(INDX)/EXPO(INDX)
 301   ARM(IS2,6)=(0.,0.)
 300  CONTINUE
 202  CALL NEWLV(N,LD,LAM)
      IF(LD.EQ.0)GO TO 201
      NHOLD=N
      GO TO 204
 1    CONTINUE
      RETURN
      END
      SUBROUTINE LAISUM(IR,N,RSG,LAM,LD,NZ,I57)
      COMPLEX ARM,FAZA,PAMP,EXPO,PAMP1
      COMMON/PSPIN/ISHA
      COMMON/AZ/ARM(600,7)
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON /CAUX/NPT,NDIV,KDIV,LAMR(8),ISG,D2W,NSW,ISG1
      COMMON/PINT/ISSTAR(76),ISSTO(75),MSTORE(2,75)
      COMMON/ADBXI/EXPO(500)
      COMMON /CCOUP/ZETA(50000),LZETA(8)
       COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/ALLC/LOCQ(8,7)
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      RMIR=CAT(IR,3)
      III=0
      IF(LAM.GT.6)III=1
      LA=LAM
      IF(LAM.GT.6)LAM=LAM-6
      DO 560 I2=1,LD
      PAMP=(0.,0.)
      M=MSTORE(1,I2)
      INDX=MSTORE(2,I2)
      ISMIN=0
      IS1=NSTART(M)
      IF(IS1.EQ.0)GO TO 560
      ISPLUS=IFIX(RMIR-CAT(IS1,3))-LAM
      IF(ISPLUS.GE.0)GOTO 150
      ISMIN=ISPLUS
      ISPLUS=0
 150  IS2=IS1+ISPLUS-1
      MRANGE=2*LAM+1+ISMIN
      IF(IS2+MRANGE.GT.NSTOP(M)) MRANGE=NSTOP(M)-IS2
      IF(MRANGE.LE.0)GOTO 560
      DO 540 I3=1,MRANGE
      IS=IS2+I3
      RMIS=CAT(IS,3)
      IF(ISO.EQ.0.AND.RMIR.GT..1.AND.RMIS.GT..1)GO TO 540
      RMU=RMIS-RMIR
      MUA=ABS(RMU)+1.1
      IF(LA.GT.6.AND.MUA.EQ.1)GO TO 540
      INDQ=LOCQ(LAM,MUA)+NPT
      NZ=NZ+1
      Z=ZETA(NZ)
      Q=ZETA(INDQ+LP7)
      IF(NDIV.NE.0)Q=ZETA(INDQ+LP7)+FLOAT(KDIV)*(ZETA(INDQ+LP7+ISG1)-
     *ZETA(INDQ+LP7))/FLOAT(NDIV)
      PAMP1=FAZA(LA,MUA,RMU,RSG)*Q*Z
      IF(ISO.EQ.0.AND.RMIR.GT..1)GO TO 545
      PAMP=PAMP1*ARM(IS,I57)+PAMP
      IF(ISO.EQ.0.AND.RMIS.GT..1)GO TO 540
 545   IF(N.EQ.M)GO TO 540
      IRS=(-1)**(IFIX(RMIR+RMIS)-ISHA+III)
      ARM(IS,6)=ARM(IS,6)+IRS*PAMP1*ARM(IR,I57)
      ISSTAR(I2)=MIN0(IS,ISSTAR(I2))
      ISSTO(I2)=MAX0(IS,ISSTO(I2))
 540  CONTINUE
      IF(N.EQ.M)GO TO 570
      ARM(IR,4)=ARM(IR,4)+PAMP*ELM(INDX)*EXPO(INDX)
      GO TO 560
 570  ARM(IR,4)=ARM(IR,4)+PAMP*ELM(INDX)
 560  CONTINUE
      LAM=LA
      RETURN
      END
       COMPLEX FUNCTION EXPON(INX,NPT,ISG,ISG1,NDIV,KDIV)
       COMPLEX EXPO1,CI,EXPOX,TCEXP
       COMMON/ADX/ADB(365)
       COMMON/CXI/XI(500)
       DATA CI/(0.,1.)/
      EXPOX=TCEXP(CI*XI(INX)*ADB(NPT)*ISG)
      EXPON=EXPOX
      IF(NDIV.EQ.0)GO TO 1
      EXPO1=TCEXP(CI*XI(INX)*ADB(NPT+ISG1)*ISG)
      EXPON=EXPOX+FLOAT(KDIV)*(EXPO1-EXPOX)/FLOAT(NDIV)
  1   CONTINUE
      RETURN
      END
      COMPLEX FUNCTION FAZA(LA,MI,RMU,RSG)
      COMPLEX CI
      DATA CI/(0.,1.)/
      IF(LA.GT.6)GO TO 100
      IEVEN=(-1)**MI
      IF(IEVEN)2,2,3
 2    FAZA=-CI
      RETURN
  3   FAZA=CMPLX(RSG,0.)
      RETURN
 100  FAZA=-CI
      IF(RMU.LT.0.)FAZA=-FAZA
      IF(LA.EQ.7)RETURN
      IF(MI.EQ.2)RETURN
      FAZA=CMPLX(RSG,0.)
      IF(RMU.LT.0.)FAZA=-FAZA
      RETURN
      END
      SUBROUTINE SETIN
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/HIPER/SH(365),CH(365)
      COMMON/ADX/ADB(365)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      DO 10 K=1,LP12
 10   ADB(K)=EPS(IEXP)*SH(K)+.03*(K-1)
      RETURN
      END
      SUBROUTINE STING(IRLD)
      COMPLEX ARM,EXPO
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/AZ/ARM(600,7)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/ADBXI/EXPO(500)
      COMMON/FLA/IFLG
      COMMON/PINT/ISSTAR(76),ISSTO(75),MSTORE(2,75)
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/CAUX/NPT,NDIV,KDIV,LAMR(8),ISG,D2W,NSW,ISG1
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/RNG/IRA(8),MAXLA
      MAXH=MAXLA
 20   ISG=-1
      N=1
      RSG=-1.
      IFLG=1
      W0=IRA(MAXLA)*.03+.03
      DO 1 J=1,ISMAX
      DO 1 JJ=1,6
   1  ARM(J,JJ)=(0.,0.)
      ARM(IRLD,5)=(1.,0.)
      DO 2 J=1,8
 2    LAMR(J)=0
      LAMR(MAXLA)=1
      NPT=IRA(MAXLA)+1
      IF(MAXLA.EQ.7.AND.IRA(2).NE.0)GO TO 111
      GO TO 112
 111  LAMR(2)=1
      NPT=NPT-1
      W0=W0-.03
 112  NDIV=0
      KDIV=0
      DO 5 J=1,4
      NPT=NPT-1
      DO 7 J1=1,LAMMAX
      LAM=LAMDA(J1)
      IF(LAMR(LAM).EQ.0)GO TO 7
      CALL NEWLV(N,LD,LAM)
      IF(LD.EQ.0)GO TO 77
      NZ=LZETA(LAM)
      LD=LDNUM(LAM,1)
      I57=5
      CALL LAISUM(IRLD,N,RSG,LAM,LD,NZ,I57)
      DO 300 MM=1,LD
      INDX=MSTORE(2,MM)
      IBG=ISSTAR(MM)
      IEND=ISSTO(MM)
      DO 301 IS2=IBG,IEND
      ARM(IS2,4)=ARM(IS2,4)+ARM(IS2,6)*ELM(INDX)/EXPO(INDX)
 301  ARM(IS2,6)=(0.,0.)
 300  CONTINUE
      GO TO 7
 77   IF(J1.NE.MAXLA)GO TO 7
      IRA(MAXLA)=-IRA(MAXLA)
      DO 88 JJ=1,LAMMAX
      LAM=LAMDA(JJ)
      IF(IRA(LAM).GT.0)GO TO 99
 88   CONTINUE
 99   MAXLA=LAMDA(JJ)
      GO TO 20
 7    CONTINUE
      IF(J.EQ.4)GO TO 4
      DO 5 I=1,ISMAX
      ARM(I,J)=ARM(I,4)
 5    ARM(I,4)=(0.,0.)
 4    CALL LAIAMP(IRLD,W0)
      MAXLA=MAXH
      DO 101 JJ=1,8
 101  IRA(JJ)=IABS(IRA(JJ))
      RETURN
      END
      SUBROUTINE LAIAMP(IR,W0)
      COMPLEX ARM,STAMP,DIS,UHUJ
      COMMON /CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/AZ/ARM(600,7)
      COMMON /CAUX/NPT,NDIV,KDIV,LAMR(8),ISG,D2W,NSW,ISG1
      COMMON /CCOUP/ZETA(50000),LZETA(8)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/CXI/XI(500)
       PPP=0.
       EPSI=EPS(IEXP)
      ERRT=EROOT(IEXP)
      RMIR=CAT(IR,3)
      DO 580 I1=1,LAMMAX
      LAM=LAMDA(I1)
      NZ=LZETA(LAM)
      IF(LAMR(LAM).EQ.0)GO TO 580
      LA=LAM
      IF(LAM.GT.6) LAM=LAM-6
      LD=LDNUM(LA,1)
      IF(LD.EQ.0)GOTO 580
      DO 560 I2=1,LD
      M=LEADF(1,I2,LA)
      INDX=MEM(1,M,LA)
      XIV=XI(INDX)
      ISMIN=0
      IS1=NSTART(M)
      IF(NSTART(M).EQ.0)GO TO 560
      ISPLUS=IFIX(RMIR-CAT(IS1,3))-LAM
      IF(ISPLUS.GE.0)GOTO 150
      ISMIN=ISPLUS
      ISPLUS=0
 150  IS2=IS1+ISPLUS-1
      MRANGE=2*LAM+1+ISMIN
      IF(IS2+MRANGE.GT.NSTOP(M)) MRANGE=NSTOP(M)-IS2
      IF(MRANGE.LE.0)GOTO 560
      DO 540 I3=1,MRANGE
      IS=IS2+I3
      NZ=NZ+1
      Z=ZETA(NZ)
      RMIS=CAT(IS,3)
      RMU=RMIS-RMIR
      MUA=ABS(RMU)+1.1
      IF(LAM.GT.6.AND.MUA.EQ.1)GO TO 540
      CALL FAZA1(LA,MUA,RMIR,RMIS,DIS,RMU)
      PM=ELM(INDX)*Z
      UHUJ=STAMP(EPSI,ERRT,XIV,.03,W0,LAM,MUA)
      ARM(IS,5)=DIS*PM*UHUJ
      PPP=PPP+TCABS(ARM(IS,5))*TCABS(ARM(IS,5))
 540  CONTINUE
 560  CONTINUE
 580  CONTINUE
      ARM(IR,5)=CMPLX(SQRT(1.-PPP),0.)
      RETURN
      END
       SUBROUTINE FAZA1(LA,MI,RMIR,RMIS,DIS,RMU)
       COMPLEX DIS,CI
       DATA CI/(0.,1.)/
      IRS=(-1)**IFIX(RMIR+RMIS)
      IF(LA.EQ.7)GO TO 4
      IEVEN=(-1)**MI
      IF(IEVEN)2,2,3
  2    DIS=-CI*IRS
      RETURN
 3    DIS=CMPLX(-FLOAT(IRS),0.)
      RETURN
 4    DIS=-CI*IRS
      IF(RMU.LT.0.)DIS=-DIS
      RETURN
      END
      SUBROUTINE TRINT(ARG,SI,CI)
      A=ARG*ARG
      IF(ARG.LT.1.)GO TO 1
      S=SIN(ARG)
      C=COS(ARG)
      F=POL4(1.,38.027246,265.187033,335.67732,38.102495,A)
      F=F/POL4(1.,40.021433,322.624911,570.23628,157.105423,A)/ARG
      G=POL4(1.,42.242855,302.757865,352.018498,21.821899,A)
      G=G/POL4(1.,48.196927,482.485984,1114.978885,449.690326,A)/A
      SI=F*C+G*S
      CI=G*C-F*S
      RETURN
 1    SI=POL4(0.,2.83446712E-5,-1.66666667E-3,.055555555,-1.,A)
      SI=SI*ARG
      CI=POL4(-3.100198413E-6,2.314814815E-4,-.0104166667,.25,0.,A)
      CI=CI-ALOG(ARG)
      RETURN
      END
      FUNCTION POL4(C0,C1,C2,C3,C4,A)
      pol4=1.
      if(abs(a).gt.1.e+9)return
      POL4=C4+A*(C3+A*(C2+A*(C1+A*C0)))
      RETURN
      END
      COMPLEX FUNCTION STAMP(EPSI,ERRT,XIV,DW,W0,LMDA,MUA)
      MI=MUA-1
      AXI=ABS(XIV)
      LA=LMDA
      IF(LMDA.EQ.7)LA=3
      IF(AXI.LT.1.E-5)GO TO 100
      EX=EXP(W0)/2.
      B=AXI*(EPSI*EX+W0)
      CALL TRINT(B,SIB,CIB)
      SB=SIN(B)/B
      CB=COS(B)/B
      BIS=SB+CIB
      BIC=CB-SIB
      BIS2=-SB/B
      BIC2=-CB/B
      DWI=-3.*DW
      EXA=EXP(DWI)
      A=AXI*(EPSI*EX*EXA+W0+DWI)
      SA=SIN(A)/A
      CA=COS(A)/A
      CALL TRINT(A,SIA,CIA)
      CIS=SA+CIA-BIS
      CIC=CA-SIA-BIC
      IF(LA.EQ.1)GO TO 125
      DWI=(BIC2-CIS+CA/A)/2.
      EXA=(BIS2+CIC+SA/A)/2.
      STAMP=CMPLX(DWI,EXA)
      GO TO 101
 125  STAMP=CMPLX(CIC,CIS)
 101  GO TO 200
 100  A=-2.*W0
      IF(LA.EQ.3)A=-W0
      EXA=EXP(A)
      DWI=3*DW
      CIC=EXA*(EXP(DWI)-1.)
      STAMP=CMPLX(CIC,0.)
      GO TO 300
 200  GO TO(51,52,53)LA
 51   IF(MI.EQ.0)FCT=.5*AXI/EPSI
      IF(MI.EQ.1)FCT=.3535533907*ERRT*AXI/EPSI
      GO TO 400
 52   IF(MI.EQ.0)FCT=.75*(3.-EPSI*EPSI)*AXI*AXI/EPSI/EPSI
      IF(MI.EQ.1)FCT=1.837117307*ERRT*AXI*AXI/EPSI/EPSI
      IF(MI.EQ.2)FCT=-.9185586535*ERRT*ERRT*AXI*AXI/EPSI/EPSI
      GO TO 400
 53   FCT=-.3535533905*ERRT*AXI*AXI
      GO TO 400
 300  GO TO(61,62,63)LA
 61   IF(MI.EQ.0)FCT=1./EPSI/EPSI
      IF(MI.EQ.1)FCT=1.414213562*ERRT/EPSI/EPSI
      GO TO 400
 62   IF(MI.EQ.0)FCT=3.*(3.-EPSI*EPSI)/EPSI/EPSI/EPSI/EPSI
      IF(MI.EQ.1)FCT=1.837117307*ERRT/EPSI/EPSI/EPSI/EPSI
      IF(MI.EQ.2)FCT=-3.674234613*ERRT*ERRT/EPSI/EPSI/EPSI/EPSI
      GO TO 400
 63   FCT=-1.414213562*ERRT/EPSI/EPSI
 400  STAMP=STAMP*FCT
      STAMP=CONJG(STAMP)
      RETURN
      END
      SUBROUTINE RESET(ISO)
      COMPLEX ARM
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/AZ/ARM(600,7)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/COEX2/NMAX,NDIM,NMAX1
      IF(ISO.EQ.0)GO TO 1
      DO 2 J=1,ISMAX
      ARM(J,1)=ARM(J,2)
      ARM(J,2)=ARM(J,3)
 2    ARM(J,3)=ARM(J,4)
      RETURN
 1    DO 3 J=1,NMAX
      IR=NSTART(J)-1
 4    IR=IR+1
      ARM(IR,1)=ARM(IR,2)
      ARM(IR,2)=ARM(IR,3)
      ARM(IR,3)=ARM(IR,4)
      IF(CAT(IR,3).LT.-.1)GO TO 4
 3    CONTINUE
      RETURN
      END
      SUBROUTINE HALF(ISO)
      COMPLEX ARM,FPOM
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/AZ/ARM(600,7)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/COEX2/NMAX,NDIM,NMAX1
      IF(ISO.EQ.0)GO TO 1
      DO 2 J=1,ISMAX
      FPOM=ARM(J,3)
      ARM(J,1)=-.0625*(ARM(J,4)+ARM(J,1))+.5625*(ARM(J,2)
     *+ARM(J,3))
      ARM(J,3)=ARM(J,3)*.75+.375*ARM(J,4)-ARM(J,2)/8.
 2    ARM(J,2)=FPOM
      RETURN
 1    DO 3 J=1,NMAX
      IR=NSTART(J)-1
 4    IR=IR+1
      FPOM=ARM(IR,3)
      ARM(IR,1)=-.0625*(ARM(IR,1)+ARM(IR,4))+.5625*
     *(ARM(IR,2)+ARM(IR,3))
      ARM(IR,3)=ARM(IR,3)*.75+.375*ARM(IR,4)-ARM(IR,2)/8.
      ARM(IR,2)=FPOM
      IF(CAT(IR,3).LT.-.1)GO TO 4
 3    CONTINUE
      RETURN
      END
      SUBROUTINE DOUBLE(ISO)
      COMPLEX ARM,FPOM
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/AZ/ARM(600,7)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/COEX2/NMAX,NDIM,NMAX1
      IF(ISO.EQ.0)GO TO 1
      DO 2 J=1,ISMAX
      FPOM=ARM(J,2)
      ARM(J,2)=-8.*ARM(J,3)+6.*ARM(J,2)+3.*ARM(J,4)
      ARM(J,1)=-16.*ARM(J,1)+9.*ARM(J,2)+9.*FPOM-ARM(J,4)
 2    ARM(J,3)=FPOM
      RETURN
 1    DO 3 J=1,NMAX
      IR=NSTART(J)-1
 4    IR=IR+1
       FPOM=ARM(IR,2)
      ARM(IR,2)=-8.*ARM(IR,3)+6.*ARM(IR,2)+3.*ARM(IR,4)
      ARM(IR,1)=-16.*ARM(IR,1)+9.*ARM(IR,2)+9.*FPOM
     *-ARM(IR,4)
      ARM(IR,3)=FPOM
      IF(CAT(IR,3).LT.-.1)GO TO 4
 3    CONTINUE
      RETURN
      END
      SUBROUTINE PATH(IRLD)
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/CLCOM8/CAT(600,3),ISMAX
      SPM=CAT(IRLD,3)
      DO 1 I=2,NMAX
      IPATH(I)=0
      IST=NSTART(I)
      IF(IST.EQ.0)GO TO 1
      ISP=NSTOP(I)
      DO 2 J=IST,ISP
      VL=CAT(J,3)
      IF(abs(vl-spm).lt.1.e-6)GO TO 3
  2   CONTINUE
      GO TO 1
  3   IPATH(I)=J
  1   CONTINUE
      IPATH(1)=IRLD
      RETURN
      END
      SUBROUTINE INTG(IEN)
      COMPLEX ARM,HOLD
      COMMON /COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/AZ/ARM(600,7)
      COMMON/RNG/IRA(8),MAXLA
      COMMON/A50/ACC50
      COMMON/CLCOM0/IFAC(75)
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/CAUX/NPT,NDIV,KDIV,LAMR(8),ISG,D2W,NSW,ISG1
      COMMON/FLA/IFLG
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/CEXC9/INTERV(50)
      INTEND=INTERV(IEN)
      D2W=.03
      NSW=1
      KAST=0
      NDIV=0
      KDIV=0
 100  IF((NPT+NSW).GT.IRA(MAXLA).AND.ISG.GT.0)RETURN
      DO 101 I=1,8
      LAMR(I)=0
      IF((NPT+NSW).LT.IRA(I))LAMR(I)=1
 101  CONTINUE
      IF(ISO.EQ.0)GO TO 290
      DO 285 IR=1,ISMAX
      ARM(IR,7)=ARM(IR,5)+D2W/24.*(55.0*ARM(IR,4)-59.0*ARM(IR,3)
     1+37.0*ARM(IR,2)-9.0*ARM(IR,1))
 285  CONTINUE
      GOTO 305
 290  DO 300 N=1,NMAX
      IR=NSTART(N)-1
 295  IR=IR+1
      ARM(IR,7)=ARM(IR,5) +D2W/24.*(55.0*ARM(IR,4)-59.0*ARM(IR,3)
     1+37.0*ARM(IR,2)-9.0*ARM(IR,1))
      MIR=CAT(IR,3)
      IR1=IR-2*MIR
      ARM(IR1,7)=IFAC(N)*ARM(IR,7)
      IF(FLOAT(MIR).LT.-0.1)GOTO 295
 300  CONTINUE
 305  CONTINUE
      NPT=NPT+NSW*ISG
      IF(NPT.LE.0)GO TO 601
      IF(NDIV.EQ.0)GO TO 325
      KDIV=KDIV+1
      IF(KDIV.GE.NDIV)GO TO 326
      GO TO 325
 326  KDIV=0
      NPT=NPT+ISG
      IF(NPT.LE.0)GO TO 601
      GO TO 325
 601  NPT=-NPT+2
      ISG=1
 325  CALL RESET(ISO)
      IFLG=1
      I57=7
      CALL AMPDER(I57)
      IF(ISO.EQ.0)GOTO 315
      DO 310 IR=1,ISMAX
      ARM(IR,5)=ARM(IR,5)+D2W/24.*(9.0*ARM(IR,4)+19.0*ARM(IR,3)
     1-5.0*ARM(IR,2)+ARM(IR,1))
 310  CONTINUE
      GOTO 330
 315  DO 311 N=1,NMAX
      IR=NSTART(N)-1
 320  IR=IR+1
      ARM(IR,5)=ARM(IR,5)+D2W/24.*(9.0*ARM(IR,4)+19.0*ARM(IR,3)
     1-5.0*ARM(IR,2)+ARM(IR,1))
      MIR=CAT(IR,3)
      IR1=IR-2*MIR
      ARM(IR1,5)=IFAC(N)*ARM(IR,5)
      IF(FLOAT(MIR).LT.-0.1)GOTO 320
 311  CONTINUE
 330  CONTINUE
      KAST=KAST+1
      IFLG=0
      I57=5
      CALL AMPDER(I57)
      IF((LAMR(2)+LAMR(3)).EQ.0)GO TO 100
      IF(KAST.LT.INTEND)GO TO 100
      KAST=0
      F=0.
      DO 401 K=1,NMAX
       IHOLD=IPATH(K)
      IF(IHOLD.EQ.0)GO TO 401
      HOLD=ARM(IHOLD,5)-ARM(IHOLD,7)
      RL=REAL(HOLD)
      RIM=AIMAG(HOLD)
      SRT=RL*RL+RIM*RIM
      F=AMAX1(F,SRT)
  401 CONTINUE
      F=SQRT(F)/14.
      IF(F.LE.ACCUR.AND.F.GE.ACC50)GO TO 100
      IF(F.LT.ACC50)GO TO 500
      CALL HALF(ISO)
      D2W=D2W/2.
      NSW=(FLOAT(NSW)+.01)/2.
      INTEND=2*INTEND
      IF(NSW.GE.1)GO TO 100
      NDIV=2*NDIV
      IF(NDIV.EQ.0)NDIV=2
      GO TO 100
 500  CALL DOUBLE(ISO)
      D2W=2.*D2W
      NSW=2*NSW
      INTEND=(FLOAT(INTEND)+.01)/2.
      IF(INTEND.EQ.0)INTEND=1
      IF(NSW.GE.1)GO TO 100
       NDIV=(FLOAT(NDIV)+.01)/2.
      IF(NDIV.GE.2)GO TO 100
      NDIV=0
      NSW=1
      GO TO 100
      END
      SUBROUTINE NEWLV(N,LD,LA)
      COMPLEX EXPO,EXPON
      COMMON /CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON /CAUX/NPT,NDIV,KDIV,LAMR(8),ISG,D2W,NSW,ISG1
      COMMON/PINT/ISSTAR(76),ISSTO(75),MSTORE(2,75)
      COMMON/ADBXI/EXPO(500)
      COMMON/FLA/IFLG
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      LD=LDNUM(LA,N)
      IF(LD.EQ.0)RETURN
      DO 1 I2=1,LD
      M=LEADF(N,I2,LA)
      ISSTAR(I2)=NSTOP(M)
      ISSTO(I2)=NSTART(M)
      MSTORE(1,I2)=M
      INDX=MEM(N,M,LA)
      MSTORE(2,I2)=INDX
      IF(IFLG.EQ.0)GO TO 1
      IF(M.EQ.N)GO TO 1
      EXPO(INDX)=EXPON(INDX,NPT,ISG,ISG1,NDIV,KDIV)
 1    CONTINUE
      RETURN
      END
      SUBROUTINE CODE7(IR,IS,N,MT,INQA,INDX)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      IAPR(INDX,1)=N
      IAPR(INDX,2)=MT
      IF(IPATH(N).EQ.0.OR.IPATH(MT).EQ.0)GO TO 6
      IDN=IR-IPATH(N)
      IDM=IS-IPATH(MT)
      ISM=IDN+IDM+3
      GO TO(1,2,3,4,5)ISM
 1    INQA=1
      RETURN
 2    INQA=2
      IF(IDN.GT.IDM)INQA=3
      RETURN
 3    INQA=4
      RETURN
 4    INQA=5
      IF(IDN.GT.IDM)INQA=6
      RETURN
 5    INQA=7
      RETURN
  6   INQA=-1
      RETURN
      END
      SUBROUTINE APRAM(IEXP,INC,INDX,IRLD,ACCA)
      COMPLEX ARM
      COMMON/AZ/ARM(600,7)
      COMMON/CXI/XI(500)
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/APRX/LERF,IDIVE(50,2)
      LERF=0
      ACCAH=ACCA
 500  I7=7
      ITM=-1
      IMG=3
      I1=1
      IF(MAGA(IEXP).NE.0)GO TO 15
      I7=4
      I1=4
      IMG=1
 15   IF(INC.EQ.0)GO TO 1
      IF(LERF.EQ.0)CALL NEWCAT(IEXP,JIDIM)
      IF(LERF.EQ.0)CALL PODZIEL(3,IEXP)
      I56=5
      DO 2 K=1,JIDIM
      ARM(K,2)=(0.,0.)
   2  ARM(K,5)=(0.,0.)
      ARM(IRLD+1,5)=(1.,0.)
 100  KTOTO=0
      LERF=0
      L1=IDIVE(IEXP,1)
      DO 300 L3=1,L1
      ACCA=ACCAH*L3/L1
      CALL POMNOZ(ACCA,1,I56,KTOTO,IMG,JIDIM)
      IF(LERF.EQ.0)GO TO 300
      CALL PODZIEL(1,IEXP)
      GO TO 500
 300  CONTINUE
      L2=IDIVE(IEXP,2)
      DO 301 L3=1,L2
      ACCA=ACCAH+ACCAH*L3/L2
      CALL POMNOZ(ACCA,2,I56,KTOTO,IMG,JIDIM)
      IF(LERF.EQ.0)GO TO 301
      CALL PODZIEL(2,IEXP)
      GO TO 500
 301  CONTINUE
      DO 3 L=1,MEMX6
      DO 3 M=I1,I7
 3    QAPR(L,1,M)=-QAPR(L,1,M)
      DO 302 L3=1,L1
      ACCA=ACCAH*2.+ACCAH*L3/L1
      CALL POMNOZ(ACCA,1,I56,KTOTO,IMG,JIDIM)
 302  CONTINUE
      ACCA=ACCAH
      DO 7 L=1,MEMX6
      DO 7 M=I1,I7
 7    QAPR(L,1,M)=-QAPR(L,1,M)
      IF(INC.EQ.0.AND.ITM.EQ.0)GO TO 1
      IF(INC.EQ.0)GO TO 200
      DO 60 JJ=2,JIDIM
 60   ARM(JJ-1,5)=ARM(JJ,5)
      RETURN
   1     ITM=ITM+1
      I56=ITM+6
      DO 10 K=1,JIDIM
 10   ARM(K,I56)=(0.,0.)
      ARM(IRLD+1,I56)=(1.,0.)
      UWA=-ITM*.0298019802+1.01
      DO 11 L=1,2
      DO 11 J=I1,I7
 11   QAPR(INDX,L,J)=QAPR(INDX,L,J)*UWA
      DO 404 J=1,JIDIM
 404  ARM(J,2)=(0.,0.)
      GO TO 100
 200  DO 20 L=1,JIDIM
      ARM(L,6)=ARM(L,6)-ARM(L,7)
 20   ARM(L,6)=50.*ARM(L,6)/ELM(INDX)
       DO 12 L=1,2
       DO 12 J=I1,I7
 12   QAPR(INDX,L,J)=QAPR(INDX,L,J)/.99
      DO 61 JJ=2,JIDIM
 61   ARM(JJ-1,6)=ARM(JJ,6)
      RETURN
      END
      SUBROUTINE NEWCAT(IEXP,JIDIM)
      COMMON/MAP/PARX(50,12,5),PARXM(50,4,10,6),XIR(6,50)
      COMMON/CXI/XI(500)
      COMMON /CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      JIDIM=NMAX+1
      IF(MAGA(IEXP).NE.0)JIDIM=3*NMAX+1
      IST=1
      DO 11 KK=1,6
      IF(MULTI(KK).EQ.0)GO TO 11
      ISTOP=MULTI(KK)-1+IST
      DO 2 K=IST,ISTOP
      XX=ABS(XI(K))
      XX=XX/XIR(KK,IEXP)
      DO 5 N=1,7,3
      IF(MAGA(IEXP).EQ.0.AND.N.NE.4)GO TO 5
      ZT=QAPR(K,1,N)
      ZT=ABS(ZT)
      XP=9.*XX
      NL=INT(XP)+1
      WG=XP-FLOAT(NL-1)
      NG=NL+1
      WL=FLOAT(NL)-XP
      A=WG*PARXM(IEXP,1,NG,KK)+WL*PARXM(IEXP,1,NL,KK)
      B=WG*PARXM(IEXP,2,NG,KK)+WL*PARXM(IEXP,2,NL,KK)
      Q1=A*ZT+B
      A=WG*PARXM(IEXP,3,NG,KK)+WL*PARXM(IEXP,3,NL,KK)
      B=WG*PARXM(IEXP,4,NG,KK)+WL*PARXM(IEXP,4,NL,KK)
      Q2=A*ZT+B
      QAPR(K,2,N)=QAPR(K,1,N)*Q2*FXIS2(K,N)
      QAPR(K,1,N)=QAPR(K,1,N)*Q1*FXIS1(K,N)
      IF(IAPR(K,1).NE.IAPR(K,2))GO TO 5
      QAPR(K,1,N)=0.
      QAPR(K,2,N)=QAPR(K,2,N)/2.
 5    CONTINUE
      IF(MAGA(IEXP).EQ.0)GO TO 2
      DO 7 N=2,6
      IF(N.EQ.4)GO TO 7
      ZT=QAPR(K,1,N)
      ZT=ABS(ZT)
      XP=4.*XX
      NL=INT(XP)+1
      WG=XP-FLOAT(NL-1)
      NG=NL+1
      WL=FLOAT(NL)-XP
      Q1=WG*PARX(IEXP,2*KK-1,NG)+WL*PARX(IEXP,2*KK-1,NL)
      Q2=WG*PARX(IEXP,2*KK,NG)+WL*PARX(IEXP,2*KK,NL)
      QAPR(K,2,N)=QAPR(K,1,N)*Q2*FXIS2(K,N)
      QAPR(K,1,N)=QAPR(K,1,N)*Q1*FXIS1(K,N)
 7    CONTINUE
  2   CONTINUE
      IST=ISTOP+1
 11   CONTINUE
      RETURN
      END
      SUBROUTINE POMNOZ(ACCA,L,IW,KTOTO,IMG,JIDIM)
      COMPLEX ARM,CI
      common/inhi/inhb
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/AZ/ARM(600,7)
      COMMON/APRX/LERF,IDIVE(50,2)
      DATA CI/(0.,-1.)/
      SIG=1.
      IF(L.NE.2)SIG=-1.
      DO 5 KK=1,JIDIM
 5    ARM(KK,1)=ARM(KK,IW)
      DO 1 K=1,100
      KTOTO=KTOTO+1
      DO 2 M=1,MEMX6
      MW1=IAPR(M,1)
      MC1=IAPR(M,2)
      IF(IPATH(MW1).EQ.0.OR.IPATH(MC1).EQ.0)GO TO 2
      MW=IPATH(MW1)+1
      MC=IPATH(MC1)+1
      IF(KTOTO.LT.ISEX(MC1))GO TO 2
      IF(IMG.EQ.1)GO TO 200
      ARM(MW,2)=ARM(MW,2)+QAPR(M,L,4)*ARM(MC,1)
      ARM(MC,2)=ARM(MC,2)+SIG*QAPR(M,L,4)*ARM(MW,1)
      ARM(MW-1,2)=ARM(MW-1,2)+QAPR(M,L,2)*ARM(MC,1)
      ARM(MC,2)=ARM(MC,2)+SIG*QAPR(M,L,2)*ARM(MW-1,1)
      ARM(MW-1,2)=ARM(MW-1,2)+QAPR(M,L,1)*ARM(MC-1,1)
      ARM(MC-1,2)=ARM(MC-1,2)+SIG*QAPR(M,L,1)*ARM(MW-1,1)
      ARM(MW,2)=ARM(MW,2)+QAPR(M,L,3)*ARM(MC-1,1)
      ARM(MC-1,2)=ARM(MC-1,2)+SIG*QAPR(M,L,3)*ARM(MW,1)
      ARM(MW,2)=ARM(MW,2)+QAPR(M,L,5)*ARM(MC+1,1)
      ARM(MC,2)=ARM(MC,2)+SIG*QAPR(M,L,6)*ARM(MW+1,1)
      ARM(MC+1,2)=ARM(MC+1,2)+SIG*QAPR(M,L,5)*ARM(MW,1)
      ARM(MW+1,2)=ARM(MW+1,2)+QAPR(M,L,6)*ARM(MC,1)
      ARM(MW+1,2)=ARM(MW+1,2)+QAPR(M,L,7)*ARM(MC+1,1)
      ARM(MC+1,2)=ARM(MC+1,2)+SIG*QAPR(M,L,7)*ARM(MW+1,1)
      GO TO 2
 200  ARM(MW,2)=ARM(MW,2)+QAPR(M,L,4)*ARM(MC,1)
      ARM(MC,2)=ARM(MC,2)+SIG*QAPR(M,L,4)*ARM(MW,1)
 2    CONTINUE
      TEST=0.
      DO 20 M=1,JIDIM
      ARM(M,1)=ARM(M,2)/K
      ARM(M,2)=(0.,0.)
      IF(L.EQ.1)GO TO 25
      ARM(M,1)=ARM(M,1)*CI
 25   ARM(M,IW)=ARM(M,IW)+ARM(M,1)
      IF(K.LE.5)GO TO 20
      U=TCABS(ARM(M,IW))
      TEST=TEST+U*U
 20   CONTINUE
      IF(ABS(TEST-1.).LT.ACCA)GO TO 30
 1    CONTINUE
      IF(inhb.eq.1)go to 30
      LERF=1
 30    CONTINUE
      RETURN
      END
      SUBROUTINE TENB(ICL,BTEN,LMAX)
      COMPLEX ARM
      DIMENSION BTEN(1200)
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/AZ/ARM(600,7)
      IHA=(-1)**IFIX(2.*SPIN(1)+.01)
      IF(ICL.NE.1)GO TO 11
      MS=16*(NMAX-1)
      DO 10 I=1,MS
  10  BTEN(I)=0.
  11  CONTINUE
      DO 1 I=2,NMAX
      MS=NSTART(I)
      IF(MS.EQ.0)GO TO 1
      MSP=NSTOP(I)
      SI=SPIN(I)
      ISI=IFIX(2.*SI+.01)
      CE=SQRT(2.*SI+1.)
      DO 2 KP=1,7,2
      K=KP-1
      KK=2*K
      IF(ISI.LT.K)GO TO 2
      ILA=-1
      DO 3 LP=1,KP
      ILA=-ILA
      L=LP-1
      LL=2*L
      IND=K*K/4+LP+(I-2)*16
      DO 4 M=MS,MSP
      MM=M
      MP=M+L
      JM=IFIX(2.01*CAT(MM,3))
      IF(MP.GT.NSTOP(I))GO TO 8
      ILG=(-1)**IFIX(SI-CAT(MP,3))
      JMP=-IFIX(2.01*CAT(MP,3))
      FC=WTHREJ(ISI,KK,ISI,JMP,LL,JM)
      ITE=1
  7   CONTINUE
      IF(ILA.EQ.1)X=REAL(ARM(MP,5))*REAL(ARM(MM,5))+
     *AIMAG(ARM(MP,5))*AIMAG(ARM(MM,5))
      IF(ILA.NE.1)X=REAL(ARM(MP,5))*AIMAG(ARM(MM,5))-
     *REAL(ARM(MM,5))*AIMAG(ARM(MP,5))
      BTEN(IND)=BTEN(IND)+X*FC*ILG
      IF(ITE.EQ.2)GO TO 4
  8   CONTINUE
      IF(IHA.EQ.1.AND.ICL.EQ.LMAX)GO TO 4
      ITE=2
      MP=MP-2*L
      IF(MP.LT.NSTART(I))GO TO 4
      JMP=IFIX(2.01*CAT(MP,3))
      JM=-JM
      FC=WTHREJ(ISI,KK,ISI,JMP,LL,JM)
      ILG=(-1)**IFIX(SI+CAT(MP,3))
      GO TO 7
  4   CONTINUE
      IF(ICL.EQ.LMAX)BTEN(IND)=BTEN(IND)*CE/(2.*SPIN(1)+1.)
  3   CONTINUE
  2   CONTINUE
  1   CONTINUE
      RETURN
      END
      SUBROUTINE TENS(BTEN)
      DIMENSION BTEN(1200)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/TCM/TETACM(50),TREP(50),DSIGS(50)
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      IX=NMAX*28
      ARG=1.570796327+TETACM(IEXP)/2.
      DO 1 I=1,IX
  1   ZETA(I)=0.
      DO 2 I=2,NMAX
      DO 3 KP=1,7,2
      K=KP-1
      K1=IFIX(FLOAT(K)/2.+.01)
      IF(K.EQ.0)GO TO 98
      DO 4 LP=1,KP
      IF(IAXS(IEXP).EQ.0.AND.LP.NE.1)GO TO 4
      INZ=(I-1)*28+K1*7+LP
      L=LP-1
      DO 5 LPP=1,KP
      IND=K*K/4+LPP+(I-2)*16
      LX=LPP-1
      LXX=LX
  13  IPH=(-1)**(L+IFIX(FLOAT(LXX)/2.))
      ZETA(INZ)=ZETA(INZ)+BTEN(IND)*IPH*DJMM(ARG,K,LX,L)
      IF(LPP.EQ.1)GO TO 5
      IF(LX.LT.0)GO TO 5
      LX=-LX
      LXX=LX-1
      GO TO 13
  5   CONTINUE
  4   CONTINUE
      GO TO 3
  98  IND=(I-2)*16+1
      INZ=(I-1)*28+1
      ZETA(INZ)=BTEN(IND)
  3   CONTINUE
  2   CONTINUE
      RETURN
      END
      FUNCTION DJMM(BETA,K,KPP,KP)
      DIMENSION DJM(525),ICZY(525)
      COMMON/IDENT/BEQ
      COMMON/CB/B(20)
      SAVE DJM
      IFZA=1
      IF(BETA.LT.0.)IFZA=(-1)**(KP+KPP)
      SK=FLOAT(K)
      UL=SK*((SK-1.)*(4.*SK+7)/6.+1.)
      LCA=INT(UL+.1)
      LOC=LCA+(2*K+1)*KP+KPP+K+1
      IF(abs(BEQ-ABS(BETA)).gt.1.e-6)GO TO 99
      IF(ICZY(LOC).EQ.1)GO TO 100
      GO TO 101
  99  BEQ=ABS(BETA)
      DO 102 ILL=1,525
  102 ICZY(ILL)=0
  101 BE=BEQ/2.
      CB=COS(BE)
      SB=SIN(BE)
      IFLA=0
      IF(BEQ.GT..01.AND.ABS(BEQ-6.2832).GT..01)IFLA=1
      IF(IFLA.EQ.1)GO TO 3
      IF(KP.EQ.KPP)GO TO 4
      DJMM=0.
      RETURN
  4   SB=1.
  3   CTB=CB*CB/SB/SB
      JA=K+KP+1
      JB=K-KP+1
      JC=K+KPP+1
      JD=K-KPP+1
      B1=B(JA)*B(JB)*B(JC)*B(JD)
      JA=KP+KPP
      JB=2*K-KP-KPP
      IF(ABS(BEQ-3.141592654).LT..01.AND.JA.LT.0)IFLA=3
      IF(IFLA.EQ.3)CB=1.
      F=(-1)**(K-KP)*(CB**JA)*(SB**JB)*SQRT(B1)
      MIS=0
      IF(JA.LT.0)MIS=-JA
      MAS=K-KPP
      IF(KPP.LT.KP)MAS=K-KP
      JA=KP+KPP+MIS+1
      JB=K-KPP-MIS+1
      JC=K-KP-MIS+1
      JD=MIS+1
      B2=B(JA)*B(JB)*B(JC)*B(JD)
      IF(IFLA.EQ.3)GO TO 10
      G=(-CTB)**MIS/B2
      DJMM=G
      JA=MIS+1
      IF(MAS.LT.JA)GO TO 2
      DO 1 J=JA,MAS
      G=-G*CTB*(K-KPP-J+1)*(K-KP-J+1)/(KP+KPP+J)/J
 1    DJMM=DJMM+G
   2  CONTINUE
      IF(IFLA.EQ.0)DJMM=G
      DJMM=DJMM*F*IFZA
      DJM(LOC)=DJMM/IFZA
      ICZY(LOC)=1
      RETURN
 10   DJMM=F*IFZA/((-SB*SB)**MIS)/B2
      DJM(LOC)=DJMM/IFZA
      ICZY(LOC)=1
      RETURN
 100  DJMM=DJM(LOC)*IFZA
      RETURN
      END
      SUBROUTINE FTBM(ICLL,CHISQ,IDR,NCALL,CHILO,BTEN)
      COMPLEX ARM
      DIMENSION JMTE(6),PROP(6),BTEN(1200)
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/TCM/TETACM(50),TREP(50),DSIGS(50)
      COMMON/CEXC0/NSTART(76),NSTOP(75)
      COMMON/CCC/eg(50),cc(50,5),NANG(200),Q(3,200,8),NICC,
     *AGELI(50,200,2)
      COMMON/ILEWY/NWR
	common/switch/jzb,ibps
      COMMON/CH1T/CHIS11
      COMMON/IGRAD/IGRD
      COMMON/CAUX0/NCM,EMMA(75)
      COMMON/CXI/XI(500)
      COMMON/MAP/PARX(50,12,5),PARXM(50,4,10,6),XIR(6,50)
      COMMON/LCZP/LFL,EMH,INM,LFL1,LFL2
      COMMON/ERRAN/KFERR
      COMMON/UWAGA/ITAK2
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/CLM/LMAX
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/CLCOM8/CAT(600,3),ISMAX
      COMMON/AZ/ARM(600,7)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/DFTB/DEVD(500),DEVU(500)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/PTH/IPATH(75),MAGA(75)
      COMMON/PRT/IPRM(20)
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/SKP/JSKIP(50)
      COMMON/LIFE/NLIFT
      COMMON/LOGY/LNY,INTR,IPS1
      ISSP=0
      CHILO=0.
      FX=2.*SPIN(1)+1.
      CHISQ=0.
      LFL=0
      CHIS1=0.
      IXX=NDIM*MEMAX+LP11
      DO 12 I1=1,IXX
 12   ZETA(I1)=0.
       DO 4372 II=1,LP6
 4372 ILE(II)=1
      ITEMP=0
      NWR=0
      IAPX=1
      DO 1 JKL=1,NEXPT
      IEXP=JKL
      IGRD=0
      LFL2=1
      IF(ITAK2.NE.-1)GO TO 605
      DO 601 LARM=1,4
      DO 601 KARM=1,LP10
 601  ARM(KARM,LARM)=(0.,0.)
 605  IFLG=0
       IF(IEXP.EQ.1)GO TO 4933
       KK=NANG(IEXP)
       DO 4934 JJJ=1,LP6
 4934 ILE(JJJ)=ILE(JJJ)+NYLDE(IEXP-1,JJJ)
 4933 LP=3
      IF(JSKIP(JKL).EQ.0)GO TO 500
      IF(MAGA(IEXP).EQ.0)LP=1
      IF(NCALL.EQ.0)GO TO 200
      IF(ICLL.EQ.4)GO TO 10
 202  LOCH=LP3*(MEMAX-1)+NMAX+LP11
      DO 6130 K=1,LOCH
 6130 ZETA(K)=0.
      CALL LOAD(IEXP,1,2,0.,JJ)
      DO 2 K=1,LMAX
      FC=2.
      IF(K.EQ.LMAX)FC=1.
      IF(FLOAT(IFIX(SPIN(1))).lt.SPIN(1))FC=2.
      LOC=0
      POLM=FLOAT(K-1)-SPIN(1)
       CALL LOAD(IEXP,3,2,POLM,JJ)
       CALL PATH(JJ)
       CALL LOAD(IEXP,2,2,POLM,JJ)
      CALL APRAM(IEXP,1,1,JJ,ACCA)
      IF(NCALL.EQ.0)GO TO 30
      IF(ICLL.EQ.3)GO TO 30
      DO 3 INDX=1,MEMX6
      CALL APRAM(IEXP,0,INDX,JJ,ACCA)
      KX=0
      DO 400 I11=1,NMAX
      IF(NSTART(I11).EQ.0)GO TO 400
      LOC=LP3*(INDX-1)+I11+LP11
      JPP=IFIX(2.*SPIN(I11)+1.)
      LPX=MIN0(LP,JPP)
      IF(ISO.NE.0)LPX=NSTOP(I11)-NSTART(I11)+1
      DO 407 LPXD=1,LPX
      KX=KX+1
 407  ZETA(LOC)=ZETA(LOC)+FC*REAL(ARM(KX,5))*REAL(ARM(KX,6))/FX+
     *FC*AIMAG(ARM(KX,5))*AIMAG(ARM(KX,6))/FX
  400 CONTINUE
 3    CONTINUE
 30   CALL TENB(K,BTEN,LMAX)
 2    CONTINUE
      IF(LOC.EQ.0)GO TO 7633
      REWIND 14
      WRITE(14,*)(ZETA(I11),I11=LP8,LOCH)
 7633 CALL TENS(BTEN)
      IF(NCALL.EQ.0)GO TO 500
      IF(ICLL.LT.2)GO TO 505
      GO TO 500
 505  LLX=28*NMAX
      DO 303 LX=1,LLX
 303  ZETA(LP9+LX)=ZETA(LX)
      IF(ICLL.EQ.1)GO TO 10
      GO TO 500
 10   IAPX=0
      ISSP=1
      CALL LOAD(IEXP,1,1,0.,JJ)
      CALL ALLOC(ACCUR)
      CALL SNAKE(IEXP,ZPOL)
      CALL SETIN
      DO 20 K=1,LMAX
      POLM=FLOAT(K-1)-SPIN(1)
      CALL LOAD(IEXP,2,1,POLM,KK)
       IF(IPRM(7).EQ.-1)WRITE(22,701)POLM,IEXP
      CALL STING(KK)
      CALL PATH(KK)
      CALL INTG(IEXP)
      CALL TENB(K,BTEN,LMAX)
      IF(IPRM(7).NE.-1)GO TO 22
      DO 21 J=1,ISMAX
      WRITE(22,702)IFIX(CAT(J,1)),CAT(J,2),
     *CAT(J,3),REAL(ARM(J,5)),AIMAG(ARM(J,5))
 21   CONTINUE
  22  CONTINUE
 20   CONTINUE
      CALL TENS(BTEN)
      IF(IPRM(7).NE.-1)GO TO 500
      DO 708 JJGG=2,NMAX
      LOCT=(JJGG-1)*28+1
 708  WRITE(22,703)JJGG,ZETA(LOCT)
      GO TO 500
 100  CONTINUE
      IF(NCALL.NE.0)GO TO 500
 200   IF(IFLG.EQ.1)GO TO 201
      IF(IFLG.EQ.2)GO TO 3111
      ITEMP=2
      IFLG=1
      GO TO 10
 201  CONTINUE
      ITEMP=1
      IFLG=2
      GO TO 202
 500  CONTINUE
      CALL CEGRY(CHISQ,ITEMP,CHILO,IDR,NWYR,ICLL,ISSP,0)
      ISSP=0
      IF(NCALL.EQ.0.AND.JSKIP(JKL).NE.0)GO TO 100
      NWR=NWR+NWYR
      IF(ICLL.GT.2.OR.JSKIP(JKL).EQ.0)GO TO 3111
      IF(IEXP.EQ.1)CHISH=CHIS11
      IF(ICLL.EQ.1)CHIS1=CHIS11
      IF(ICLL.EQ.0)CHIS1=CHISQ
      LFL2=0
      IGRD=1
      IF(ITAK2.EQ.-1)LFL=1
      REWIND 14
      READ(14,*)(ZETA(I11),I11=LP8,LOCH)
      DO 305 LARM=1,4
      DO 305 KARM=1,LP10
 305  ARM(KARM,LARM)=(0.,0.)
      CHISX=0.
      LLX=28*NMAX
      DO 600 LIX=1,LLX
 600  ZETA(LP9+LIX)=ZETA(LIX)
      CALL CEGRY(CHISX,ITEMP,CHILO,IDR,NWYR,0,0,1)
      DO 301 KNM=1,MEMAX
      INM=KNM
      CHISX=0.
      EMH=ELM(INM)
      ELM(INM)=1.05*EMH
      LCC=LP3*(INM-1)+LP11
      DO 304 LST=2,NMAX
      WZ=ZETA(LST+LCC)
      INPX=(LST-1)*28
      DO  307 JY=1,4
      INP=INPX+(JY-1)*7
      IF(JY.EQ.1)PR=ZETA(LP13+INP)+1.E-12
      JMF=2*JY-1
      IF(IAXS(IEXP).EQ.0)JMF=1
      DO 306 JM=1,JMF
      INP=INP+1
 306  ZETA(INP)=ZETA(INP+LP9)*(1.+.1*EMH*WZ/PR)
 307  CONTINUE
 304  CONTINUE
      CALL CEGRY(CHISX,ITEMP,CHILO,IDR,NWYR,0,0,0)
      ELM(INM)=EMH
 301  CONTINUE
      IF(ITAK2.NE.-1.OR.LFL1.EQ.0)GO TO 3111
      IF(IPRM(17).EQ.0)GO TO 800
      KMT=IABS(IPRM(17))
      WRITE(22,801)IEXP
      NLIN=(NMAX-2)/6+1
      NREST=NMAX-1-6*(NLIN-1)
      DO 802 ILIN=1,NLIN
      NPOZ=6
      IF(ILIN.EQ.NLIN)NPOZ=NREST
      INPO=(ILIN-1)*6+2
      INKO=INPO+NPOZ-1
      LPIT=0
      DO 888 LM=INPO,INKO
      LPIT=LPIT+1
 888  JMTE(LPIT)=LM
      WRITE(22,803)(JMTE(LM),LM=1,LPIT)
      WRITE(22,804)(ZETA(LP13+(JPZ-1)*28),JPZ=INPO,INKO)
      DO 805 JMT=1,KMT
      LPUT=0
      DO 819 LS=INPO,INKO
      LPUT=LPUT+1
      PROP(LPUT)=0.
      DO 820 LM=1,MEMX6
      INZZ=LS+LP3*(LM-1)+LP11
      INZZZ=LP13+(LS-1)*28
      IF(abs(ZETA(INZZZ)).lt.1.e-20)ZETA(INZZZ)=1.E-20
      VAL=2.*ELM(LM)*ZETA(INZZ)/ZETA(INZZZ)
      AVAL=ABS(VAL)
      IF(AVAL.LE.ABS(PROP(LPUT)))GO TO 820
      PROP(LPUT)=VAL
      LMH=LM
      JMTE(LPUT)=LM
  820 CONTINUE
      IZZZ=(LMH-1)*LP3+LP11+LS
      ZETA(IZZZ)=0.
 819  CONTINUE
      WRITE(22,807)(JMTE(LCOU),PROP(LCOU),LCOU=1,NPOZ)
 805  CONTINUE
 802  CONTINUE
      REWIND 14
      READ(14,*)(ZETA(I11),I11=LP8,LOCH)
 801  FORMAT(1X///20X,10HEXPERIMENT,11X,1I2,5X
     *,24HD(LOG(P))/D(LOG(ME)) MAP/20X,52(1H-)///)
 803  FORMAT(5X,5HLEVEL,6(8X,1I2,9X))
 807  FORMAT(10X,6(2X,1H(,1X,1I3,1X,1E8.2,1H),2X))
 804  FORMAT(1X,9HEXC.PROB.,6(5X,1E10.4,4X))
      IF(IPRM(17).LT.0)GO TO 3111
 800   LFL=0
      WRITE(22,610)IEXP
      ILE1=ILE(1)+NYLDE(IEXP,1)-1
       ILE3=ILE(1)
      LICZ=0
      DO 611 ILE2=ILE3,ILE1
      LICZ=LICZ+1
      IDEC=IY(ILE2,1)
      IF(IDEC.GT.1000)IDEC=IDEC/1000
      LUU=6*LICZ-5
      JK=(LUU-1)/LP10+1
      KK=LUU-LP10*(JK-1)
      KK6=KK+5
 611  WRITE(22,612)KSEQ(IDEC,3),KSEQ(IDEC,4),(IFIX(REAL(ARM
     *(KKX,JK))),AIMAG(ARM(KKX,JK)),KKX=KK,KK6)
 610  FORMAT(10X,10HEXPERIMENT,1X,1I2/10X,19HD(LOG(Y)/D(LOG(ME)),//)
 612  FORMAT(2X,1I2,2H--,1I2,5X,6(1H(,1I3,2X,1E8.2,1H),3X))
 3111 CONTINUE
 1    CONTINUE
      IF(ITAK2.EQ.-1.AND.ICLL.LT.2)ITAK2=0
      IF(NCALL.EQ.0)GO TO 999
      IF(ICLL.GT.2)GO TO 501
      IF(ICLL.EQ.1)CALL CEGRY(CHISQ,ITEMP,CHILO,IDR
     *,NOWR,7,ISSP,0)
 501  CALL BRANR(CHISQ,NWR,CHILO)
      CALL MIXR(NWR,0,CHISQ,CHILO)
      CALL CHMEM(NWR,CHISQ,CHILO)
      NWR=NWR+NLIFT
      CHISQ=CHISQ/NWR
      IF(INTR.EQ.0)GO TO 999
      CHX=CHISQ
      CHISQ=CHILO
      CHILO=CHX
 999  CONTINUE
      RETURN
 701  FORMAT(1X//40X,21HEXCITATION AMPLITUDES//10X,2HM=
     *,1F4.1,5X,10HEXPERIMENT,1X,1I2//5X,5HLEVEL,2X,4HSPIN
     *,2X,1HM,5X,14HREAL AMPLITUDE,2X,19HIMAGINARY AMPLITUDE
     *//)
 702  FORMAT(7X,1I2,3X,1F4.1,2X,1F4.1,2X,1E14.6,2X,1E14.6)
 703  FORMAT(2X,5HLEVEL,1X,1I2,10X,10HPOPULATION,1X,1E14.6)
      END
      SUBROUTINE MINI(CHISQ,CHIOK,NPTL,CONV,IMODE,IDR,XTEST,IPS,
     *IS,JJH,BTEN)
      DIMENSION IPM(10),BTEN(1200),GRADP(500)
      COMMON/DUMM/GRAD(500),HLMLM(500),ELMH(500)
	common/switch/jzb,ibps
      COMMON/ILEWY/NWR
      COMMON/CH1T/CHIS11
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/MINNI/IMIN,LNORM(50)
      COMMON/UWAGA/ITAK2
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      COMMON/DFTB/DEVD(500),DEVU(500)
      COMMON/PRT/IPRM(20)
      COMMON/LCZP/LFL,EMH,INM,LFL1,LFL2
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/SEL/KVAR(500)
      COMMON/FIT/LOCKF,NLOCK,IFBFL,LOCKS,DLOCK
      COMMON/ERRAN/KFERR
      COMMON/LOGY/LNY,INTR,IPS1
      COMMON/ERCAL/JENTR,ICS
      DO 1006 I=1,MEMAX
 1006 GRADP(I)=0.
      ICOUNT=0
      LNM=0
      LNY=0
      INTR=0
      METF=0
      LFL1=0
      NCALL=0
      ITAK2=0
      IF(IMODE.LT.2000)GO TO 2
      ICL1=1
      IF(IMODE.GE.2100)METF=1
      IF((IMODE-2000-100*METF).GE.10)LNM=1
      IF((IMODE-2000-100*METF-10*LNM).EQ.1)LNY=1
      ICL2=4
      IF(IPS.EQ.0)GO TO 3
      IF(IPS.EQ.1)GO TO 333
      IF(IPRM(4).LT.0)ITAK2=-2
      GO TO 334
 333  CONTINUE
      IF(IPRM(4).EQ.-1)ITAK2=-2
 334  ICL1=4
      IF(ITAK2.EQ.-2)ICL1=1
      IF(ICL1.EQ.4)GO TO 111
      GO TO 3
 2    ICL1=0
      ICL2=3
      IF(IMODE.GE.1100)METF=1
      IF((IMODE-1000-100*METF).GE.10)LNM=1
      IF((IMODE-1000-100*METF-10*LNM).EQ.1)LNY=1
      IF(JENTR.EQ.1)GO TO 111
      IF(ICS.EQ.0)GO TO 3
      REWIND 11
      DO 3711 JNM=1,LP4
 3711 read(11)(CORF(JNM,KH2),KH2=1,LP6)
      ICS=0
      GO TO 111
 3    CALL FTBM(0,CHISS,IDR,0,CHL,BTEN)
      REWIND 11
      DO 3712 JNM=1,LP4
 3712 write(11)(CORF(JNM,KH2),KH2=1,LP6)
      IF(IPS1.EQ.0)RETURN
 111  NOFLG=0
      NCALL=1
 1    CONTINUE
      SUMHT=0.
      IF(LNY.EQ.1)INTR=1
      LFL1=1
      ITAK2=ITAK2+1
      ICOUNT=ICOUNT+1
      IF(ICOUNT.GT.NPTL)GO TO 997
      IF(ITAK2.EQ.IPRM(4))ITAK2=-1
      IF(ITAK2.NE.-1)GO TO 700
      IF(KFERR.EQ.1)GO TO 700
      CALL FTBM(3,CHD,IDR,1,CHL,BTEN)
      CHIS11=CHD*NWR
      CALL FTBM(ICL1,CHISQ,IDR,NCALL,CHILO,BTEN)
 700  IF(IPS.EQ.1)RETURN
      IF(ICL1.EQ.1)CALL FTBM(4,CHISQ,IDR,NCALL,CHILO,BTEN)
      IF(IPRM(8).NE.-1.AND.IPRM(13).NE.-1)GO TO 400
      IF(IPRM(8).EQ.-1)IPRM(8)=-2
      IF(IPRM(13).EQ.-1)IPRM(13)=-2
      CALL FTBM(4,CCD,IDR,NCALL,CHL,BTEN)
      IF(IPS.EQ.2)RETURN
 400  CALL FTBM(3,CHIS12,IDR,NCALL,CHILO,BTEN)
      IF(ICL1.EQ.0)CHISQ=CHIS12
      UXA=CHISQ
      IF(INTR.EQ.1)UXA=CHILO
      IPAS=0
      IF(UXA.LT.CHIOK)CHISQ=UXA
      IF(UXA.LT.CHIOK)GO TO 998
 201  INO=1
      IF(METF.EQ.1)IPAS=IPAS+1
      IF(IFBFL.EQ.1)INO=2
      DO 207 JJJ=1,INO
      DO 702 JNM=1,MEMAX
      GRAD(JNM)=0.
      IF(IVAR(JNM).NE.1.AND.IVAR(JNM).NE.2)GO TO 702
      DO 803 JCOUP=1,MEMAX
 803  ELMH(JCOUP)=ELM(JCOUP)
      DO 800 JCOUP=1,MEMAX
      IF(JNM.EQ.JCOUP)GO TO 801
      IF(IVAR(JCOUP).LT.1000)GO TO 800
      JCP=IVAR(JCOUP)-1000
      IF(JCP.NE.JNM)GO TO 800
      IF(IVAR(JNM).EQ.0)GO TO 800
 801  FLT=1.01
      IF(JJJ.EQ.2)FLT=.99
      ELM(JCOUP)=ELMH(JCOUP)*FLT
 800  CONTINUE
      CALL FTBM(3,CHIS13,IDR,NCALL,CHX,BTEN)
      IF(JJJ.EQ.1)HLMLM(JNM)=CHIS13
      IF(IFBFL.EQ.1.AND.JJJ.EQ.1)GO TO 208
      IF(JJJ.EQ.2)CHIS12=CHIS13
      GRAD(JNM)=100.*(HLMLM(JNM)-CHIS12)/ELMH(JNM)
      IF(IFBFL.EQ.1)GRAD(JNM)=GRAD(JNM)/2.
      IF(LNM.EQ.1)GRAD(JNM)=GRAD(JNM)*ABS(ELMH(JNM))
 208  CONTINUE
      DO 804 JCOUP=1,MEMAX
 804  ELM(JCOUP)=ELMH(JCOUP)
 702  CONTINUE
 207  CONTINUE
      IF(KFERR.NE.1)GO TO 701
      GRAD(JJH)=0.
      IF(IS.NE.1.OR.ICOUNT.NE.1)GO TO 701
      WRITE(33,*)(NWR*GRAD(JNM),JNM=1,MEMAX)
 701  CONTINUE
      IF(METF.EQ.1.AND.IPAS.EQ.2)GO TO 862
      SUMG2=0.
      DO 860 JNM=1,MEMAX
      IF(IVAR(JNM).NE.1.AND.IVAR(JNM).NE.2)GO TO 860
      SUMG2=SUMG2+GRAD(JNM)*GRAD(JNM)
 860  CONTINUE
      IF(SUMG2.LT.1.E-10)GO TO 999
      SUMG2=SQRT(SUMG2)
      DO 888 JNM=1,MEMAX
 888  GRAD(JNM)=GRAD(JNM)/SUMG2
      IF(METF.EQ.0)GO TO 861
      DM=0.
      DO 852 JNM=1,MEMAX
      IF(IVAR(JNM).NE.2.AND.IVAR(JNM).NE.1)GO TO 852
      DM=DM+ELM(JNM)*ELM(JNM)*GRAD(JNM)*GRAD(JNM)
 852  CONTINUE
      DM=SQRT(DM)
      DO 851 JNM=1,MEMAX
      DEVD(JNM)=GRAD(JNM)
       DEVU(JNM)=ELM(JNM)
      SEL=DM*GRAD(JNM)/20.
      IF(LNM.EQ.1)SEL=SEL*ABS(ELM(JNM))
 851  ELM(JNM)=ELM(JNM)-SEL
      IF(IFBFL.EQ.0)CALL FTBM(3,CHIS12,IDR,NCALL,CHX,BTEN)
      GO TO 201
 862  CONTINUE
      DO 863 JNM=1,MEMAX
      ELM(JNM)=DEVU(JNM)
 863  CONTINUE
      SHL=DM/20./SUMG2
      SUMG1=0.
      DO 864 JNM=1,MEMAX
      GRAD(JNM)=(DEVD(JNM)*SUMG2-GRAD(JNM))/SHL
 864  SUMG1=SUMG1+GRAD(JNM)*GRAD(JNM)
      SUMG1=SQRT(SUMG1)
      P=0.
      DO 865 JNM=1,MEMAX
      GRAD(JNM)=GRAD(JNM)/SUMG1
      DEVU(JNM)=ELM(JNM)
      SEL=DM*GRAD(JNM)/100.
      IF(LNM.EQ.1)SEL=SEL*ABS(DEVU(JNM))
      P=P+DEVD(JNM)*GRAD(JNM)
 865  ELM(JNM)=ELM(JNM)+SEL
      CALL FTBM(3,CHIS13,IDR,NCALL,CHX,BTEN)
      SHL=DM/100.
      DO 889 JNM=1,MEMAX
      SEL=DM*GRAD(JNM)/50.
      IF(LNM.EQ.1)SEL=SEL*ABS(DEVU(JNM))
 889  ELM(JNM)=ELM(JNM)-SEL
      CALL FTBM(3,CHIS12,IDR,NCALL,CHX,BTEN)
      Q=(CHIS12+CHIS13-2.*CHISQ)/SHL/SHL
      A0=Q*SUMG2/SUMG1-P
      A1=P*P-1.
      SUMG1=SQRT(A0*A0+A1*A1+2.*A0*A1*P)
      DO 866 JNM=1,MEMAX
      ELM(JNM)=DEVU(JNM)
      GRAD(JNM)=(GRAD(JNM)*A1+DEVD(JNM)*A0)/SUMG1
  866 CONTINUE
 861  LFL1=0
      IF(LNM.EQ.0)GO TO 894
      DO 893 JNM=1,MEMAX
 893  GRAD(JNM)=GRAD(JNM)*ABS(ELM(JNM))
 894  CONTINUE
      SUMG1=0.
      DO 895 JNM=1,MEMAX
 895  SUMG1=SUMG1+GRAD(JNM)*GRAD(JNM)
      SUMG1=SQRT(SUMG1)
      DO 896 JNM=1,MEMAX
 896  GRAD(JNM)=GRAD(JNM)/SUMG1
      IF(LNY.EQ.1)CHISQ=CHILO
      IF(NOFLG.EQ.0)CHIRF=CHISQ
      NOFLG=1
      CHIL=CHISQ
      IF(KFERR.EQ.1)GOTO 666
      IF(MOD(ICOUNT,IPRM(5)).EQ.0.OR.ICOUNT.EQ.1)WRITE(22,54)CHISQ
c      WRITE(*,54)CHISQ
      IF(MOD(ICOUNT,IPRM(6)).NE.0)GO TO 666
      WRITE(22,601)
      NLINN=MEMAX/10+1
      DO 602 JLIN=1,NLINN
      JSA=(JLIN-1)*10+1
      DO 608 JIN=1,10
 608  IPM(JIN)=JSA+JIN-1
      JST=MIN0(JSA+9,MEMAX)
      JPR=MIN0(10,MEMAX-JSA+1)
      WRITE(22,603)(IPM(JIN),JIN=1,JPR)
      WRITE(22,604)(GRAD(JIN),JIN=JSA,JST)
 602  CONTINUE
 603  FORMAT(5X,10(5X,1I3,4X))
 604  FORMAT(5X,10(1X,1E10.4,1X)/)
 601  FORMAT(20X,8HGRADIENT//)
 666  CONTINUE
      IF(CHIL.LT.CHIOK)GO TO 998
      DO 5 L=1,MEMAX
      HLMLM(L)=ELM(L)
 5    CONTINUE
      DO 22 L=1,MEMAX
      IF(ABS(GRAD(L)).GT.DLOCK.OR.LOCKS.NE.1.OR.ICOUNT
     *.NE.1.OR.IVAR(L).GT.999.OR.IVAR(L).EQ.0)GO TO 22
      IF(KFERR.NE.1)
     *KVAR(L)=0
      IF(KFERR.NE.1)
     *WRITE(22,227)L,GRAD(L)
      IVAR(L)=0
 227  FORMAT(1X,14HMATRIX ELEMENT,1X,1I3,1X,6HLOCKED,3X,
     *11HDERIVATIVE=,1E14.6)
 22   CONTINUE
      ISTEC=0
 500  DO 66 J=1,MEMAX
 66   ELMH(J)=ELM(J)
      ISTEC=ISTEC+1
      CMAX=0.
      INTR=0
      INMX=1
      DO 777 IHT=1,MEMAX
      IF(ABS(GRAD(IHT)).LE.CMAX)GO TO 777
      CMAX=ABS(GRAD(IHT))
      INMX=IHT
 777  CONTINUE
      HT=.01*ABS(ELM(INMX))/CMAX
      MVFL=0
      IF(ICOUNT.EQ.1.OR.ISTEC.NE.1)GO TO 1000
      XKAT=0.
      DO 1001 J=1,MEMAX
 1001 XKAT=XKAT+GRAD(J)*GRADP(J)
      DO 1002 J=1,MEMAX
 1002 GRADP(J)=GRAD(J)
      IF(XKAT.LT..8)GO TO 1000
      A=0.
      DO 1003 J=1,MEMAX
      IF(IVAR(J).EQ.0.OR.IVAR(J).GT.999)GO TO 1003
      A=AMAX1(A,ABS(GRAD(J)))
      IF(abs(A-ABS(GRAD(J))).lt.1.e-9)IIN=J
 1003 CONTINUE
      WRITE(22,215)IIN
      IVAR(IIN)=0
      GRAD(IIN)=0.
      GRADP(IIN)=0.
 1000 CONTINUE
 300  DO 6 J=1,MEMAX
 6    ELM(J)=ELMH(J)-HT*GRAD(J)
      DO 7 J=1,MEMAX
      IF(IVAR(J).LT.1000)GO TO 7
      INDX1=IVAR(J)-1000
      ELM(J)=ELM(INDX1)*SA(J)
 7    CONTINUE
      IF(MVFL.EQ.0)GO TO 13
      CALL LIMITS
      GO TO 100
 13   CALL FTBM(ICL2,CHISP,IDR,NCALL,CHILO,BTEN)
      DO 8 J=1,MEMAX
   8  ELM(J)=2.*ELMH(J)-ELM(J)
      CALL FTBM(ICL2,CHISF,IDR,NCALL,CHILO,BTEN)
      C=(CHISP+CHISF-2.*CHIL)/HT/HT
      B=(CHISP-CHISF)/HT/2.
      DL=B*B-2.*C*CHIL
      IF(DL.GT.0.)GO TO 11
      F1=B
      F2=C
      GO TO 12
   11 F1=CHIL
      F2=B
 12   MVFL=1
      IF(ABS(F2).LT.1.E-10)GO TO 14
      HT=-F1/F2
      GO TO 300
  14  HT=1.
      GO TO 300
 100  CALL FTBM(ICL2,CHISQ,IDR,NCALL,CHILO,BTEN)
      IF(CHISQ.GE.CHIL)GO TO 16
      CHIL=CHISQ
      SUMHT=SUMHT+HT
      IF(ABS(HT/SUMHT).LT..01)GO TO 17
      GO TO 500
 16   HT=HT/2.
      IF(ABS(HT).LT.CONV)GO TO 17
      GO TO 300
 17   CRIT=0.
      DO 18 JJJ=1,MEMAX
 18   CRIT=CRIT+(ELM(JJJ)-HLMLM(JJJ))**2
      CRIT=SQRT(CRIT)
      IF(CRIT.LT.CONV)GO TO 999
      if(chisq.lt.chiok)go to 998
      RFK=CHIRF/CHISQ
      IF(RFK.GT.XTEST.AND.ICOUNT.LT.NPTL)GO TO 3
      GO TO 1
 997  CONTINUE
        IF(KFERR.EQ.1)RETURN
      IF(IPS.EQ.0)WRITE(22,50)NPTL
      IF(IPS.EQ.0)WRITE(22,54)CHIL
      INTR=0
      RETURN
  998 CHIL=CHISQ
      IF(IPS.EQ.0)WRITE(22,51)ICOUNT
      IF(IPS.EQ.0)WRITE(22,54)CHIL
      RETURN
 999  IF(LOCKF.EQ.0)GO TO 210
      DO 212 KKK=1,NLOCK
      A=0.
      IIN=1
      DO 211 JJJ=1,MEMAX
      IF(IVAR(JJJ).EQ.0.OR.IVAR(JJJ).GT.999)GO TO 211
      A=AMAX1(A,ABS(GRAD(JJJ)))
      IF(abs(A-ABS(GRAD(JJJ))).lt.1.e-9)IIN=JJJ
 211  CONTINUE
      IVAR(IIN)=0
      WRITE(22,215)IIN
 212  CONTINUE
      ITF=0
      DO 213 JJJ=1,MEMAX
      IF(IVAR(JJJ).GT.999)GO TO 213
      IF(IVAR(JJJ).NE.0)ITF=ITF+1
 213  CONTINUE
      if(itf.ne.1)go to 3233
      metf=0
      write(22,3234)
 3234 format(2x,'Warning - only one matrix element free',//2x,
     *'Mode reset to single gradient, execution continues',/)
 3233 continue
      IF(ITF.NE.0)GO TO 1
      WRITE(22,214)
 214  FORMAT(1X/////5X,5H*****,2X,
     *27HALL MATRIX ELEMENTS LOCKED!,2X,5H*****/////)
      INTR=0
      RETURN
 215  FORMAT(1X/5X,14HMATRIX ELEMENT,1X,1I3,1X,7HLOCKED!)
 210  CONTINUE
      IF(CHISQ.LT.CHIL)GO TO 698
      DO 699 JJJ=1,MEMAX
 699  ELM(JJJ)=ELMH(JJJ)
 698  CONTINUE
      IF(KFERR.EQ.1)RETURN
      IF(IPS.EQ.0)WRITE(22,53)ICOUNT,CRIT
      IF(IPS.EQ.0)WRITE(22,54)AMIN1(CHIL,CHISQ)
      INTR=0
      RETURN
 50   FORMAT(5X,42HMINIMIZATION STOPPED-NUMBER OF STEPS NPTL=,
     *1I5,1X,8HEXCEEDED)
 51   FORMAT(5X,7HAT STEP,1X,1I5,1X,
     *25HCHISQ CRITERION FULFILLED)
 53   FORMAT(5X,7HAT STEP,1X,1I5
     *,21HCONVERGENCE ACHIEVED(,1E14.6,1H))
 54   FORMAT(5X,10H*** CHISQ=,1E14.6,1X,3H***)
      END
      SUBROUTINE CEGRY(CHISQ,ITEMP,CHILO,IDR,NWYR,ICALL,ISSP,IREDV)
      CHARACTER*4 WUPL,WAR,OP2
      DIMENSION LIC(32),LTH(500)
     *,CNR(32,50)
      common/clust/iclust(50,200),lastcl(50,20),sumcl(20,500)
     *,irawex(50)
      COMMON/ODCH/DEV(500)
	common/tutaj/part(32,50,2),partl(32,50,2)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/TRA/DELTA(500,3),ENDEC(500),ITMA(50,200),ENZ(200)
      COMMON/BREC/BETAR(50)
      COMMON/DIMX/DIX(4),ODL(200)
      COMMON/VAC/VACDP(3,75),QCEN,DQ,XNOR,AKS(6,75),IBYP
      COMMON/CINIT/cnor(32,75),innr,mcfix,cnor1(32,75),cnor2(32,75)
      COMMON/PRT/IPRM(20)
      COMMON/LIFE/NLIFT
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/IGRAD/IGRD
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON/MINNI/IMIN,LNORM(50)
      COMMON/LCZP/LFL,EMH,INM,LFL1,LFL2
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/YTEOR/YGN(500),YGP(500),IFMO
      COMMON/SEL/KVAR(500)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/CCC/eg(50),cc(50,5),NANG(200),Q(3,200,8),NICC,
     *AGELI(50,200,2)
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/WARN/SGW,SUBCH1,SUBCH2,IWF
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/SKP/JSKIP(50)
      COMMON/TRB/ITS
      COMMON/TCM/TETACM(50),TREP(50),DSIGS(50)
      common/cccds/ndst(50)
	OP2='    '
      IFXD=0
      TETRC=TREP(IEXP)
      IF(ICALL.NE.4.OR.IPRM(13).NE.-2)GO TO 400
      IPRM(13)=0
      WRITE(22,402)
      DO 401 JPC=1,NEXPT
      K=ndst(JPC)
 401  WRITE(22,403)JPC,(CNOR(L,JPC),L=1,K)
 402  FORMAT(1X//20X,23HNORMALIZATION CONSTANTS//
     *2X,10HEXPERIMENT,5X,15HDETECTORS(1-32))
 403  FORMAT(1X,1I2,2X,32(1E10.4,1X))
      WRITE(22,602)
      DO 601 JPC=1,NEXPT
      IF(abs(CNR(1,JPC)).lt.1.e-9)CNR(1,JPC)=1.
      K=ndst(JPC)
      WRITE(22,403)JPC,(CNR(L,JPC)/CNR(1,JPC),L=1,K)
 601  CONTINUE
 602  FORMAT(1X//20X,40HRECOMMENDED RELATIVE GE(LI) EFFICIENCIES//
     *2X,10HEXPERIMENT)
 400  CONTINUE
      DO 824 JPC=1,LP6
 824  LIC(JPC)=0
      IF(ICALL.EQ.7)GO TO 201
      IF(ITEMP.EQ.0)GO TO 100
      IFXD=1
      IF(ITEMP.NE.2)IFXD=0
      NWYR=1
      CALL DECAY(CCD,0,CCC)
      FI0=FIEX(IEXP,1)
      FI1=FIEX(IEXP,2)
      NA=NANG(IEXP)
      do 6257 k=1,lp2
      do 6257 kj=1,20
 6257 sumcl(kj,k)=0
      k9=0
      DO 105 K=1,NA
      GTH=AGELI(IEXP,K,1)
      FIGL=AGELI(IEXP,K,2)
      fm=(fi0+fi1)/2.
      CALL ANGULA(YGN,IDR,IFXD,FI0,FI1,TETRC,GTH,FIGL,K,OP2)
      IF(IFMO.EQ.0)GO TO 7232
      ID=ITMA(IEXP,K)
      D=ODL(ID)
      RX=D*SIN(GTH)*COS(FIGL-FM)-.25*SIN(TETRC)*COS(FM)
      RY=D*SIN(GTH)*SIN(FIGL-FM)-.25*SIN(TETRC)*SIN(FM)
      RZ=D*COS(GTH)-.25*COS(TETRC)
      RL=SQRT(RX*RX+RY*RY+RZ*RZ)
      SF=D*D/RL/RL
      THC=TACOS(RZ/RL)
      FIC=ATAN2(RY,RX)
      CALL ANGULA(YGP,IDR,IFXD,FI0,FI1,TETRC,THC,FIC,K,OP2)
      DO 7233 IXL=1,IDR
      IXM=KSEQ(IXL,3)
      TFAC=TAU(IXM)
      IF(TFAC.GT.1.E+4)GO TO 7235
 7233 YGN(IXL)=YGN(IXL)+.01199182*TFAC*BETAR(IEXP)*
     *(SF*YGP(IXL)-YGN(IXL))
 7235 IFMO=0
      WRITE(22,7236)
 7236 FORMAT(1X,/,2X,23HDURING THE MINIMIZATION,1X,
     *60HIT WAS NECESSARY TO SWITCH OFF THE TIME-OF-FLIGHT CORRECTION)
 7232 CONTINUE
      if(irawex(iexp).eq.0)go to 6252
      ipd=itma(iexp,k)
      do 6253 l=1,idr
      decen=endec(l)
      cocos=sin(tetrc)*sin(gth)*cos(fm-figl)+cos(tetrc)*
     *cos(gth)
      decen=decen*(1.+betar(iexp)*cocos)
      call effix(ipd,decen,effi)
      ygn(l)=ygn(l)*effi
 6253 continue
      inclus=iclust(iexp,k)
      if(inclus.eq.0)go to 6252
      do 6254 l=1,idr
 6254 sumcl(inclus,l)=sumcl(inclus,l)+ygn(l)
      if(k.ne.lastcl(iexp,inclus))go to 105
      do 6255 l=1,idr
 6255 ygn(l)=sumcl(inclus,l)
 6252 continue
      k9=k9+1
      IYEX=NYLDE(IEXP,K9)+ILE(K9)-1
       ILE2=ILE(K9)
      DO 102 L=ILE2,IYEX
      IF(JSKIP(IEXP).EQ.0)GO TO 102
      IDC=IY(L,K9)
      IF(IDC.LT.1000)GO TO 612
      IDC=IDC/1000
      LL1=IY(L,K9)-IDC*1000
      YGN(IDC)=YGN(IDC)+YGN(LL1)
 612  IF(ITEMP.EQ.1)GO TO 101
      CORF(L,K9)=YGN(IDC)
      IF(IMIN.GT.1.OR.L.NE.IYEX)GO TO 103
      CNOR(K9,IEXP)=YEXP(K9,L)/YGN(IDC)
 103  GO TO 102
  101 CONTINUE
      CORF(L,K9)=CORF(L,K9)/(YGN(IDC)+1.E-24)
 102  CONTINUE
 105   CONTINUE
      RETURN
 100  CONTINUE
      NWYR=0
      IF(IGRD.EQ.1)GO TO 1777
      IF(IEXP.EQ.1)SUMPR=0.
      IF(IEXP.EQ.1)SUM3=0.
      DO 107 JJ=1,LP6
      DO 107 JK=1,2
      PARTL(JJ,IEXP,JK)=0.
 107  PART(JJ,IEXP,JK)=0.
 1777 CALL DECAY(CHISQ,NLIFT,CHILO)
      IF(ICALL.NE.4.OR.IPRM(14).NE.-1)GO TO 414
      IF(IEXP.EQ.NEXPT)IPRM(14)=0
      WRITE(22,462)
      WRITE(22,467)IEXP
      DO 463 IVA=2,NMAX
  463 WRITE(22,464)IVA,(VACDP(II,IVA),II=1,3)
  462 FORMAT(1X//20X,35HVACUUM DEPOLARIZATION COEFFICIENTS //)
  467 FORMAT(5X,10HEXPERIMENT,1X,1I2/5X,5HLEVEL,10X,2HG2,10X,
     *2HG4,10X,2HG6/)
  464 FORMAT(7X,1I2,9X,3(1F6.4,6X))
  414 CONTINUE
      FI0=FIEX(IEXP,1)
      FI1=FIEX(IEXP,2)
      NA=NANG(IEXP)
      do 8257 k=1,lp2
      do 8257 k9=1,20
 8257 sumcl(k9,k)=0.
      k9=0
      DO 2 K=1,NA
      GTH=AGELI(IEXP,K,1)
      FIGL=AGELI(IEXP,K,2)
      IFXD=0
      fm=(fi0+fi1)/2.
      IF(ICALL.EQ.4)IFXD=1
      CALL ANGULA(YGN,IDR,IFXD,FI0,FI1,TETRC,GTH,FIGL,K,OP2)
      IF(IFMO.EQ.0)GO TO 8232
      ID=ITMA(IEXP,K)
      D=ODL(ID)
      RX=D*SIN(GTH)*COS(FIGL-FM)-.25*SIN(TETRC)*COS(FM)
      RY=D*SIN(GTH)*SIN(FIGL-FM)-.25*SIN(TETRC)*SIN(FM)
      RZ=D*COS(GTH)-.25*COS(TETRC)
      RL=SQRT(RX*RX+RY*RY+RZ*RZ)
      SF=D*D/RL/RL
      THC=TACOS(RZ/RL)
      FIC=ATAN2(RY,RX)
      CALL ANGULA(YGP,IDR,IFXD,FI0,FI1,TETRC,THC,FIC,K,OP2)
      DO 8233 IXL=1,IDR
      IXM=KSEQ(IXL,3)
      TFAC=TAU(IXM)
 8233 YGN(IXL)=YGN(IXL)+.01199182*TFAC*BETAR(IEXP)*
     *(SF*YGP(IXL)-YGN(IXL))
 8232 CONTINUE
      if(irawex(iexp).eq.0)go to 8252
      ipd=itma(iexp,k)
      do 8253 l=1,idr
      decen=endec(l)
      cocos=sin(tetrc)*sin(gth)*cos(fm-figl)+cos(tetrc)*
     *cos(gth)
      decen=decen*(1.+betar(iexp)*cocos)
      call effix(ipd,decen,effi)
      ygn(l)=ygn(l)*effi
 8253 continue
      inclus=iclust(iexp,k)
      if(inclus.eq.0)go to 8252
      do 8254 l=1,idr
 8254 sumcl(inclus,l)=sumcl(inclus,l)+ygn(l)
      if(k.ne.lastcl(iexp,inclus))go to 2
      do 8255 l=1,idr
 8255 ygn(l)=sumcl(inclus,l)
 8252 continue
      k9=k9+1
      IF(ICALL.NE.4.OR.IPRM(8).NE.-2)GO TO 410
      WRITE(22,404)IEXP,K9
 404  FORMAT(1X//5X,47HCALCULATED AND EXPERIMENTAL YIELDS   EXPERIMENT
     *,1X,1I2,1X,8HDETECTOR,1X,1I2//6X,
     *2HNI,5X,2HNF,7X,2HII,8X,2HIF,9X,11HENERGY(MEV),
     *6X,4HYCAL,8X,4HYEXP,7X,9HPC. DIFF.,
     *2X,13H(YE-YC)/SIGMA)
 410  LU=ILE(K9)
      DO 460 IABC=1,LP2
 460  LTH(IABC)=0
      DO 3 L=1,IDR
      NI=KSEQ(L,3)
      NF=KSEQ(L,4)
      IF(L.EQ.IY(LU,K9).OR.L.EQ.(IY(LU,K9)/1000))GO TO 300
      IF(JSKIP(IEXP).EQ.0)YGN(IDRN)=1.E+10
      RY=YGN(L)/YGN(IDRN)
      IF(ICALL.NE.4.OR.IPRM(8).NE.-2)GO TO 700
      WUPL='    '
      IF(RY.GT.UPL(K9,IEXP).AND.LTH(L).EQ.0)WUPL='UPL!'
      IF(IPRM(16).EQ.0.AND.WUPL.EQ.'    ')GO TO 700
      IF(WUPL.EQ.'    ')
     *WRITE(22,405)NI,NF,SPIN(NI),SPIN(NF),ENDEC(L),
     *YGN(L)*CNOR(K9,IEXP),WUPL
      IF(WUPL.EQ.'    ')GO TO 700
      SGM=(RY-UPL(K9,IEXP))/UPL(K9,IEXP)
      WRITE(22,406)NI,NF,SPIN(NI),
     *SPIN(NF),ENDEC(L),YGN(L)*CNOR(K9,IEXP),UPL(K9,IEXP)
     **CNOR(K9,IEXP)
     **YGN(IDRN),100.*(1.-YGN(L)/UPL(K9,IEXP)/YGN(IDRN)),
     *SGM,WUPL
      SUBCH1=SUBCH1+SGM*SGM
 700  CONTINUE
      IF(RY.LT.UPL(K9,IEXP).OR.LTH(L).EQ.1)GO TO 3
      CHISQ=CHISQ+(RY-UPL(K9,IEXP))*(RY-UPL(K9,IEXP))
     */UPL(K9,IEXP)/UPL(K9,IEXP)
      CHILO=CHILO+ALOG(RY/UPL(K9,IEXP))**2
      IF(IWF.EQ.0)GO TO 504
      WRITE(22,503)IEXP,NI,NF,RY/UPL(K9,IEXP)
 503  FORMAT(5X,13HWARNINIG-EXP. ,1I2,2X,7HTRANS. ,1I2,2H--,
     *1I2,5X,27HEXCEEDS UPPER LIMIT (RATIO=,1E14.6,1H))
 504  CONTINUE
      GO TO 3
 300  CONTINUE
      IFDU=0
      LIC(K9)=LIC(K9)+1
      LICZ=LIC(K9)
      NWYR=NWYR+1
      WF=CORF(LU,K9)
       IF(ICALL.EQ.4)WF=1.
       IF(ICALL.EQ.1.AND.ISSP.EQ.1)WF=1.
      IF(IY(LU,K9).LT.1000)GO TO 1527
      IFDU=1
      L1=IY(LU,K9)/1000
      L1=IY(LU,K9)-1000*L1
      YGN(L)=YGN(L)+YGN(L1)
      LTH(L1)=1
      IF(ICALL.NE.4.OR.IPRM(8).NE.-2)GO TO 701
      WAR='    '
      SGM=(YEXP(K9,LU)-YGN(L)*CNOR(K9,IEXP))/DYEX(K9,LU)
      IF(ABS(SGM).GE.SGW)WAR='*?!*'
      NI1=KSEQ(L1,3)
      NF1=KSEQ(L1,4)
      WRITE(22,408)NI,NI1,NF,NF1,SPIN(NI),SPIN(NI1),SPIN(NF),
     *SPIN(NF1),ENDEC(L),ENDEC(L1),YGN(L)*CNOR(K9,IEXP),
     *YEXP(K9,LU),100.*(YEXP(K9,LU)-YGN(L)*
     *CNOR(K9,IEXP))/YEXP(K9,LU),
     *SGM,WAR
      SUBCH1=SUBCH1+SGM*SGM
 701  CONTINUE
 1527 RY=YGN(L)*WF*CNOR(K9,IEXP)-YEXP(K9,LU)
      IF(IFDU.EQ.1)GO TO 702
      IF(ICALL.NE.4.OR.IPRM(8).NE.-2)GO TO 702
      WAR='    '
      SGM=(YEXP(K9,LU)-YGN(L)*CNOR(K9,IEXP))/DYEX(K9,LU)
      IF(ABS(SGM).GE.SGW)WAR='*?!*'
      WRITE(22,406)NI,NF,SPIN(NI),SPIN(NF),ENDEC(L),
     *YGN(L)*CNOR(K9,IEXP),YEXP(K9,LU),
     *100.*(YEXP(K9,LU)-YGN(L)*CNOR(K9,IEXP))/YEXP(K9,LU),
     *SGM,WAR
      SUBCH1=SUBCH1+SGM*SGM
 702  CONTINUE
      RYS=RY*RY
      IF(IGRD.EQ.1)CHISQ=CHISQ+RYS/DYEX(K9,LU)/DYEX(K9,LU)
      IF(K9.EQ.1.AND.IREDV.EQ.1)DEV(LICZ)=RY
      IF(IREDV.EQ.1)GO TO 5
      IF(LFL.NE.1)GO TO 5
      IF(K9.NE.1)GO TO 5
      LUU=6*LICZ-5
      JK=(LUU-1)/LP10+1
      KK=LUU-LP10*(JK-1)
      RIK=DEV(LICZ)+YEXP(K9,LU)
      SGM=-DEV(LICZ)/DYEX(K9,LU)
      IF(ITS.EQ.1.AND.KVAR(INM).NE.0)WRITE(17,*)NI,NF,SGM,
     *YGN(L)*CNOR(K9,IEXP)/DYEX(K9,LU)
      IF(ITS.EQ.1.AND.INM.EQ.1)WRITE(15,*)IEXP,RIK/CNOR(1,IEXP),
     *CNOR(1,IEXP),DYEX(K9,LU),YEXP(K9,LU)
      CALL SIXEL(RIK,RY,EMH,JK,KK,INM,LICZ)
 5    CONTINUE
      IF(IGRD.EQ.1)GO TO 1003
      IF(JSKIP(IEXP).EQ.0)GO TO 1003
      DL=DYEX(K9,LU)*DYEX(K9,LU)
      PART(K9,IEXP,1)=PART(K9,IEXP,1)+YGN(L)*YGN(L)
     **WF*WF/DL
      PART(K9,IEXP,2)=PART(K9,IEXP,2)-2.*YGN(L)*WF
     **YEXP(K9,LU)/DL
      SUMPR=SUMPR+YEXP(K9,LU)*YEXP(K9,LU)/DL
      PARTL(K9,IEXP,1)=PARTL(K9,IEXP,1)+YEXP(K9,LU)*YEXP(K9,LU)/DL
      PARTL(K9,IEXP,2)=PARTL(K9,IEXP,2)+ALOG(WF*YGN(L)/YEXP(K9,LU))
     **YEXP(K9,LU)*YEXP(K9,LU)/DL
      SUM3=SUM3+YEXP(K9,LU)*YEXP(K9,LU)*ALOG(WF*YGN(L)/YEXP(K9,LU))
     ***2/DL
 1003 LU=LU+1
 3    CONTINUE
      IF(IEXP.EQ.NEXPT)IWF=0
      IF(ICALL.NE.4.OR.IPRM(8).NE.-2) GO TO 2
      WRITE(22,620) SUBCH1-SUBCH2
      SUBCH2=SUBCH1
 620  FORMAT(1X/50X,17HCHISQ SUBTOTAL = ,E14.6)
 2    CONTINUE
      IF(IGRD.EQ.1)RETURN
      IF(IEXP.NE.NEXPT)RETURN
      IF(ICALL.EQ.1)RETURN
 201  DO 202 JJ=1,NEXPT
      IF(JSKIP(JJ).EQ.0)GO TO 202
      KC=ndst(JJ)
      DO 600 JK=1,KC
      CNR(JK,JJ)=-.5*PART(JK,JJ,2)/PART(JK,JJ,1)
      IF(INNR.EQ.0)GO TO 600
      CNOR(JK,JJ)=CNR(JK,JJ)
 600  CONTINUE
      IF(INNR.EQ.1)GO TO 202
      D=0.
      G=0.
      DO 203 JJ1=JJ,NEXPT
      IF(LNORM(JJ1).NE.JJ)GO TO 203
      K=ndst(JJ1)
      DO 204 JK=1,K
      D=D+YNRM(JK,JJ1)*PART(JK,JJ1,1)*YNRM(JK,JJ1)
 204  G=G-.5*YNRM(JK,JJ1)*PART(JK,JJ1,2)
 203  CONTINUE
      IF(LNORM(JJ).NE.JJ)GO TO 202
      CNOR(1,JJ)=G*YNRM(1,JJ)/D
      K=ndst(JJ)
      IF(K.EQ.1)GO TO 202
      DO 205 JK=2,K
 205  CNOR(JK,JJ)=YNRM(JK,JJ)*CNOR(1,JJ)/YNRM(1,JJ)
 202  CONTINUE
      IF(INNR.EQ.1)GO TO 209
      DO 206 JJ=1,NEXPT
      IF(LNORM(JJ).EQ.JJ)GO TO 206
      IW=LNORM(JJ)
      K=ndst(JJ)
      DO 207 JK=1,K
 207  CNOR(JK,JJ)=CNOR(1,IW)*YNRM(JK,JJ)/YNRM(1,IW)
 206  CONTINUE
 209  CONTINUE
      IF(ICALL.EQ.7)CHISQ=0.
      DO 501 JJ=1,NEXPT
      K=ndst(JJ)
      DO 502 JK=1,K
	if(mcfix.eq.0)cnor(jk,jj)=cnor1(jk,jj)
      CHILO=CHILO+PARTL(JK,JJ,1)*ALOG(CNOR(JK,JJ))**2+
     *PARTL(JK,JJ,2)*2.*ALOG(CNOR(JK,JJ))
 502  CHISQ=CHISQ+CNOR(JK,JJ)*CNOR(JK,JJ)*PART(JK,JJ,1)
     *+CNOR(JK,JJ)*PART(JK,JJ,2)
 501  CONTINUE
      CHISQ=CHISQ+SUMPR
      CHILO=CHILO+SUM3
      RETURN
 408  FORMAT(4X,1I2,1H+,1I2,2H--,1I2,1H+,1I2,3X,1F4.1,
     *1H+,1F4.1,2H--,1F4.1,1H+,1F4.1,3X,1F6.4,1H+,1F6.4,
     *2X,1E9.4,6X,1E9.4,3X,1F6.1,5X,1F4.1,10X,1A4)
 406  FORMAT(6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,1F6.4,
     *6X,1E9.4,6X,1E9.4,3X,1F6.1,5X,1F4.1,10X,1A4)
 405  FORMAT(6X,1I2,5X,1I2,7X,1F4.1,6X,1F4.1,9X,1F6.4,
     *6X,1E9.4,10X,1A4)
      END
      SUBROUTINE FAKP
      COMMON/FAKUL/IP(26),IPI(26),KF(101,26),PILOG(26)
      DO 99 I=1,26
      X=FLOAT(IP(I))
      PILOG(I)=ALOG(X)
   99 CONTINUE
      DO 1 L=1,26
      KF(1,L)=0
    1 KF(2,L)=0
      DO 3 K=3,101
      CALL PRIM(K-1)
      DO 2 I=1,26
    2 KF(K,I)=KF(K-1,I)+IPI(I)
    3 CONTINUE
      RETURN
      END
      SUBROUTINE PRIM(N)
      COMMON/FAKUL/IP(26),IPI(26),KF(101,26),PILOG(26)
      NNK=N
      DO 3 I=1,26
      NNI=NNK
      IPI(I)=0
    1 NNI=NNI/IP(I)
      IF(IP(I)*NNI.NE.NNK) GOTO 3
      IPI(I)=IPI(I)+1
      NNK=NNI
      GOTO1
    3 CONTINUE
      RETURN
      END
        SUBROUTINE SEQ(IDR)
        COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/TRA/DELTA(500,3),ENDEC(500),ITMA(50,200),ENZ(200)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
        COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/CATLF/FP(4,500,3),GKP(4,500,2),KLEC(75)
      M6=0
      DO 5 L=1,6
  5   M6=M6+MULTI(L)
      IDECAY=0
      IDR=0
      DO 515 L=1,LP3
 515  KLEC(L)=0
      DO 913 K=1,LP2
      DO 914 J=1,3
      DO 915 L=1,4
      FP(L,K,J)=0.
      IF(J.EQ.3)GO TO 915
      GKP(L,K,J)=0.
 915  CONTINUE
 914  DELTA(K,J)=0.
 913   CONTINUE
        DO 1 N=1,NMAX
        TAU(N)=EN(N)
  1     CONTINUE
        DO 2 N=1,NMAX
        EMAX=0.
         DO 3 J=1,NMAX
        IF(TAU(J).LT.EMAX)  GOTO 3
        EMAX=TAU(J)
        JSAVE=J
  3     CONTINUE
      DO 100 IS=1,NMAX
      DO 101 LA=1,8
      IF(LA.GT.3.AND.LA.NE.7.AND.LA.NE.8)GO TO 101
      LD=LDNUM(LA,IS)
      IF(LD.EQ.0)GO TO 101
      DO 102 IR=1,LD
      M=LEADF(IS,IR,LA)
      IF(M.NE.JSAVE.AND.IS.NE.JSAVE)GO TO 102
      IF(IS.EQ.JSAVE.AND.EN(M).GE.EN(IS))GO TO 102
      IF(M.EQ.JSAVE.AND.EN(IS).GE.EN(M))GO TO 102
      INDX=MEM(IS,M,LA)
      IDECAY=IDECAY+1
      KSEQ(IDECAY,1)=M
      KSEQ(IDECAY,2)=IS
      KSEQ(IDECAY,3)=INDX
      KSEQ(IDECAY,4)=LA+10
      IF(EN(M).GT.EN(IS))GO TO 102
      KSEQ(IDECAY,1)=IS
      KSEQ(IDECAY,2)=M
 102  CONTINUE
 101  CONTINUE
 100  CONTINUE
        TAU(JSAVE)=-1.
  2     CONTINUE
      DO 300 L=1,IDECAY
      ISTR1=0
      IF(KSEQ(L,4).LT.10)GO TO 300
      ISTR2=0
      N=KSEQ(L,1)
      M=KSEQ(L,2)
      INX=KSEQ(L,3)
      LA=KSEQ(L,4)-10
      EGA=EN(N)-EN(M)
      TWOI=1./SQRT(2.*SPIN(N)+1.)
      SPINI=SPIN(N)+.001
      SPINF=SPIN(M)+.001
      EGS=SQRT(EGA)*TWOI
      JS=L+1
      LA1=0
      INX1=0
      DO 301 J=JS,IDECAY
      IF(KSEQ(J,4).LT.10)GO TO 301
      N1=KSEQ(J,1)
        M1=KSEQ(J,2)
      IF(N1.NE.N.OR.M1.NE.M)GO TO 301
      INX1=KSEQ(J,3)
      LA1=KSEQ(J,4)-10
      KSEQ(J,4)=KSEQ(J,4)-10
 301  CONTINUE
      KSEQ(L,4)=KSEQ(L,4)-10
      IDR=IDR+1
      MULE=0
      MULM=0
      NOB=1
 700  IF(LA.GT.3)GO TO 400
      GO TO(500,501,502)LA
 400  LA=LA-6
      IF(LA.EQ.2)GO TO 401
      DELTA(IDR,2)=4.1952*EGA*EGS
      MULM=1
      ISTR2=4
      GO TO 600
 401  DELTA(IDR,2)=.0368*EGA*EGA*EGS
      MULM=2
      ISTR2=5
      GO TO 600
 500  DELTA(IDR,1)=399.05*EGA*EGS
      MULE=1
      ISTR1=1
      GO TO 600
 501  DELTA(IDR,1)=3.4928*EGS*EGA*EGA
      MULE=2
      ISTR1=2
      GO TO 600
 502  DELTA(IDR,1)=.02391*EGA*EGA*EGA*EGS
      MULE=3
      ISTR1=3
 600  CONTINUE
      IF(NOB.EQ.2)GO TO 303
      IF(MULE.EQ.1)GO TO 304
      NOB=NOB+1
      IF(LA.GT.3)INX1=INX
      IF(LA1.EQ.0)GO TO 304
      LA=LA1
      GO TO 700
 304  INX1=0
 303  DELTA(IDR,3)=DELTA(IDR,1)*DELTA(IDR,2)
      DELTA(IDR,1)=DELTA(IDR,1)*DELTA(IDR,1)
      DELTA(IDR,2)=DELTA(IDR,2)*DELTA(IDR,2)
      KSEQ(IDR,1)=INX
      KSEQ(IDR,2)=INX1
      KSEQ(IDR,3)=N
      KSEQ(IDR,4)=M
      IF(INX.LE.M6)GO TO 6
      KSEQ(IDR,2)=INX
      KSEQ(IDR,1)=0
  6   ENDEC(IDR)=EN(N)-EN(M)
      DO 8 MK=1,7,2
      KPA=MK/2+1
      K=MK-1
      IF(MULE.LT.3.AND.K.EQ.6)GO TO 8
      GKP(KPA,IDR,1)=GF(K,SPINI,SPINF,MULE)*DELTA
     *(IDR,1)*(1.+CONV(EGA,ISTR1))
      GKP(KPA,IDR,2)=GF(K,SPINI,SPINF,MULM)*DELTA
     *(IDR,2)*(1.+CONV(EGA,ISTR2))
      FP(KPA,IDR,1)=F(K,SPINI,SPINF,MULE,MULE)
     **DELTA(IDR,1)
      FP(KPA,IDR,3)=F(K,SPINI,SPINF,MULM,MULE)*
     *DELTA(IDR,3)
      FP(KPA,IDR,2)=F(K,SPINI,SPINF,MULM,MULM)*
     *DELTA(IDR,2)
 8    CONTINUE
      DELTA(IDR,1)=DELTA(IDR,1)*( 1.+CONV(EGA,ISTR1))
      DELTA(IDR,2)=DELTA(IDR,2)*(CONV(EGA,ISTR2)+1.)
      KLEC(N)=KLEC(N)+1
 300  CONTINUE
      NMAX1=0
      DO 71 N=1,NMAX
      IF(KLEC(N).NE.0)NMAX1=NMAX1+1
 71   CONTINUE
        RETURN
         END
        FUNCTION GF(K,SJI,SJF,L)
      GF=0.
      IF(L.EQ.0)RETURN
      IX=IFIX(SJI+SJF+.0001)
        I=IX+L+K
        PHASE=1.
        IF(I/2*2.NE.I)  PHASE=-1.
        KZ=K*2
        JIZ=SJI*2
        JFZ=SJF*2
        LZ=L*2
        GF=PHASE*SQRT((JIZ+1.)*(JFZ+1.))
     1      *WSIXJ(JIZ,JIZ,KZ,JFZ,JFZ,LZ)
        RETURN
        END
        FUNCTION F(K,SJI,SJF,L1,L2)
      F=0.
      IF((L1*L2).EQ.0)RETURN
      IX=IFIX(SJI+SJF+.0001)
        L=IX-1
        PHASE=1.
        IF(L/2*2.NE.L)  PHASE=-1.
        KZ=K*2
        JIZ=SJI*2
        JFZ=SJF*2
        L1Z=L1*2
        L2Z=L2*2
        F=PHASE*SQRT((L1Z+1.)*(L2Z+1.)*(JIZ+1.)*(KZ+1.))
     1          *WTHREJ(L1Z,L2Z,KZ,2,-2,0)
     2          *WSIXJ(JIZ,JIZ,KZ,L2Z,L1Z,JFZ)
        RETURN
        END
        FUNCTION CONV(EGA,N)
      DIMENSION CPO(51),CPO1(51)
      COMMON/CCC/eg(50),cc(50,5),NANG(200),Q(3,200,8),NICC,
     *AGELI(50,200,2)
      IF(N.EQ.0)GO TO 10
        IF(abs(CC(1,N)).lt.1.e-9) GO TO 10
      NEN=4
      DO 20 J=1,NICC
      IF(EGA.LE.EG(J))GO TO 11
 20   CONTINUE
 11   N1=J-2
      IF(N1.LT.1)N1=1
      IF((J+1).GT.NICC)N1=N1-1
      IF(NICC.GT.4)GO TO 12
      N1=1
      NEN=NICC
 12    CONTINUE
      DO 1 J=1,NEN
      CPO(J)=CC(N1+J-1,N)
      CPO1(J)=EG(N1+J-1)
 1    CONTINUE
        CALL LAGRAN(CPO1,CPO,4,1,EGA,CV,2,1)
      CONV=CV
        RETURN
 10   CONV=0.0
      RETURN
        END
      FUNCTION WTHREJ(J1,J2,J3,M1,M2,M3)
      DIMENSION JVORA(26)
      COMMON/FAKUL/IP(26),IPI(26),KF(101,26),PILOG(26)
      WTHREP=0.E+00
      JJHA=(J1+J2-J3)/2+1
      JJHB=(J1-J2+J3)/2+1
      JJHC=(-J1+J2+J3)/2+1
      IF((JJHA.LT.1).OR.(JJHB.LT.1).OR.(JJHC.LT.1).OR.
     *((M1+M2+M3).NE.0))GO TO 8
      JJHD=(J1+J2+J3+4)/2
      JMAX=MAX0(J1,J2,J3)
      IF(JMAX.EQ.J1) GOTO1
      IF(JMAX.EQ.J2) GOTO2
      IF(JMAX.EQ.J3) GOTO3
    1 JJ1=J2
      JJ2=J3
      JJ3=J1
      MM1=M2
      MM2=M3
      MM3=M1
      GOTO4
    2 JJ1=J3
      JJ2=J1
      JJ3=J2
      MM1=M3
      MM2=M1
      MM3=M2
      GOTO4
    3 JJ1=J1
      JJ2=J2
      JJ3=J3
      MM1=M1
      MM2=M2
      MM3=M3
    4 JMA=(JJ1+MM1)/2
      JMB=(JJ1-MM1)/2
      JMC=(JJ2+MM2)/2
      JMD=(JJ2-MM2)/2
      JME=(JJ3+MM3)/2
      JMF=(JJ3-MM3)/2
      JABC=(JJ1+JJ2-JJ3)/2
      JABM=(JJ2-JJ3-MM1)/2
      JBMA=(JJ1+MM2-JJ3)/2
      IZMIN=MAX0(JABM,JBMA,0)
      IZMAX=MIN0(JABC,JMB,JMC)
      NMAX=MAX0(JJHD,IZMAX+1)
      DO 9 N=1,26
      IF(IP(N).GE.NMAX) GOTO10
    9 CONTINUE
      GOTO8
   10 DO 5 JLP=1,N
      JTA=KF(JJHA,JLP)+KF(JJHB,JLP)+KF(JJHC,JLP)-KF(JJHD,JLP)
      JTB=KF(JMA+1,JLP)+KF(JMB+1,JLP)+KF(JMC+1,JLP)
      JTC=KF(JMD+1,JLP)+KF(JME+1,JLP)+KF(JMF+1,JLP)
      JVORA(JLP)=JTA+JTB+JTC
    5 CONTINUE
      VORZ=-1.E+00
      IF(2*(IZMIN/2).EQ.IZMIN) VORZ=+1.E+00
      IF(IZMIN.GT.IZMAX) GOTO8
      DO7 IZ=IZMIN,IZMAX
      QSUMLO=0.E+00
      IZA=IZ+1
      IZB=JABC+1-IZ
      IZC=JMB+1-IZ
      IZD=JMC+1-IZ
      IZE=IZ-JABM+1
      IZF=IZ-JBMA+1
      DO6 JLP=1,N
      IZEXP=JVORA(JLP)-2*KF(IZA,JLP)-2*KF(IZB,JLP)-2*KF(IZC,JLP)
     *-2*KF(IZD,JLP)-2*KF(IZE,JLP)-2*KF(IZF,JLP)
      SUMLO=IZEXP
      QSUMLO=QSUMLO+SUMLO*PILOG(JLP)*(.5E+00)
    6 CONTINUE
      ZUTHRE=VORZ* EXP(QSUMLO)
      WTHREP=WTHREP+ZUTHRE
    7 VORZ=-VORZ
      JVO=JJ1-JJ2-MM3
      IF(4*(JVO/4).EQ.JVO) GOTO8
      WTHREP=-WTHREP
    8 WTHREJ=WTHREP
      RETURN
      END
      FUNCTION WSIXJ(J1,J2,J3,L1,L2,L3)
      DIMENSION ISUMFA(26),IVORFA(26)
      COMMON/FAKUL/IP(26),IPI(26),KF(101,26),PILOG(26)
      WSIXP=0.E+00
      IF(((J1+J2-J3).LT.0).OR.((J1-J2+J3).LT.0).OR.((-J1+J2+J3).LT.0))
     1 GOTO8
      IF(((J1+L2-L3).LT.0).OR.((J1-L2+L3).LT.0).OR.((-J1+L2+L3).LT.0))
     1 GOTO8
      IF(((L1+J2-L3).LT.0).OR.((L1-J2+L3).LT.0).OR.((-L1+J2+L3).LT.0))
     1 GOTO8
      IF(((L1+L2-J3).LT.0).OR.((L1-L2+J3).LT.0).OR.((-L1+L2+J3).LT.0))
     1 GOTO8
      KQA=(J1+J2-J3)/2
      KQB=(J1-J2+J3)/2
      KQC=(J2+J3-J1)/2
      KQD=(J1+J2+J3)/2
      KRA=(J1+L2-L3)/2
      KRB=(J1-L2+L3)/2
      KRC=(L2+L3-J1)/2
      KRD=(J1+L2+L3)/2
      KSA=(L1+J2-L3)/2
      KSB=(L1-J2+L3)/2
      KSC=(J2+L3-L1)/2
      KSD=(L1+J2+L3)/2
      KTA=(L1+L2-J3)/2
      KTB=(L1-L2+J3)/2
      KTC=(L2+J3-L1)/2
      KTD=(L1+L2+J3)/2
      IZMIN=MAX0(KQD,KRD,KSD,KTD)
      KUA=KQA+KTA+J3
      KUB=KSC+KTC+L1
      KUC=KRB+KTB+L2
      IZMAX=MIN0(KUA,KUB,KUC)
      IF(IZMIN.GT.IZMAX) GOTO8
      NMAX=MAX0(IZMAX+2,KQD+2,KRD+2,KSD+2,KTD+2)
      DO 3 N=1,26
      IF(IP(N).GE.NMAX) GOTO4
    3 CONTINUE
      GOTO8
    4 VORZ=-1.E+00
      IF(2*(IZMIN/2).EQ.IZMIN) VORZ=+1.E+00
      DO1 IRL=1,N
      IVA=KF(KQA+1,IRL)+KF(KQB+1,IRL)+KF(KQC+1,IRL)-KF(KQD+2,IRL)
      IVB=KF(KRA+1,IRL)+KF(KRB+1,IRL)+KF(KRC+1,IRL)-KF(KRD+2,IRL)
      IVC=KF(KSA+1,IRL)+KF(KSB+1,IRL)+KF(KSC+1,IRL)-KF(KSD+2,IRL)
      IVD=KF(KTA+1,IRL)+KF(KTB+1,IRL)+KF(KTC+1,IRL)-KF(KTD+2,IRL)
      IVORFA(IRL)=IVA+IVB+IVC+IVD
    1 CONTINUE
      DO7 IZ=IZMIN,IZMAX
      SUMLO=0.E+00
      IZA=IZ+2
      IZB=IZ-KQD+1
      IZC=IZ-KRD+1
      IZD=IZ-KSD+1
      IZE=IZ-KTD+1
      IZF=KUA-IZ+1
      IZG=KUB-IZ+1
      IZH=KUC-IZ+1
      DO5 IRJ=1,N
      ISA=2*KF(IZA,IRJ)-2*KF(IZB,IRJ)-2*KF(IZC,IRJ)
      ISB=-2*KF(IZD,IRJ)-2*KF(IZE,IRJ)-2*KF(IZF,IRJ)
      ISC=IVORFA(IRJ)-2*KF(IZG,IRJ)-2*KF(IZH,IRJ)
      ISUMFA(IRJ)=ISA+ISB+ISC
      QSUMFA=ISUMFA(IRJ)
      SUMLO=SUMLO+QSUMFA*PILOG(IRJ)
    5 CONTINUE
      QSUMLO=(.5E+00)*SUMLO
      ZUSIX=EXP(QSUMLO)*VORZ
      WSIXP=WSIXP+ZUSIX
    7 VORZ=-VORZ
    8 WSIXJ=WSIXP
      RETURN
       END
      SUBROUTINE LAGRAN(X,Y,NDATA,IPC,XX,YY,ISCAL,IRC)
      DIMENSION X(51),Y(51),W(51),ARH(51,51)
      GO TO(50,60,70,80)IRC
 50   DO 1 I=1,NDATA
      W(I)=1.
      DO 2 J=1,NDATA
      IF(I.EQ.J)GO TO 2
      W(I)=W(I)*(XX-X(J))/(X(I)-X(J))
 2    CONTINUE
 1    CONTINUE
 60   YY=0.
      DO 11 J=1,NDATA
      Y1=Y(J)
 11   YY=YY+W(J)*FUNC(Y1,ISCAL)
      YY=FUNC1(YY,ISCAL)
      RETURN
 70   DO 21 I=1,NDATA
      T=1.
      DO 22 J=1,NDATA
      IF(I.EQ.J)GO TO 22
      T=T*(XX-X(J))/(X(I)-X(J))
 22   CONTINUE
 21   ARH(IPC,I)=T
 80   YY=0.
      DO 31 J=1,NDATA
      Y1=Y(J)
 31   YY=YY+ARH(IPC,J)*FUNC(Y1,ISCAL)
      YY=FUNC1(YY,ISCAL)
      RETURN
      END
      FUNCTION FUNC(Y,I)
      GOTO (1,2,3),I
  1   FUNC=Y
      RETURN
 2    CONTINUE
      IF(Y.LT.1.E-12)Y=1.E-12
      FUNC=ALOG(Y)
      RETURN
  3   FUNC=SQRT(Y)
      RETURN
      END
      FUNCTION FUNC1(Y,I)
      GOTO (1,2,3),I
  1   FUNC1=Y
      RETURN
  2   FUNC1=EXP(Y)
      RETURN
  3   FUNC1=Y*Y
      RETURN
      END
        SUBROUTINE GKVAC(IL)
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/BREC/BETAR(50)
      COMMON/GGG/AVJI,GAMMA,XLAMB,TIMEC,GFAC,FIEL,POWER
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON/GVAC/GKI(3),SUM(3)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/VAC/VACDP(3,75),QCEN,DQ,XNOR,AKS(6,75),IBYP
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
       COMMON/THTAR/ITTE(50)
       IF(abs(XLAMB).lt.1.e-9)GO TO 5
       IF(ITTE(IEXP).EQ.0)GO TO 1
  5    CONTINUE
       DO 2 I=1,3
  2    VACDP(I,IL)=1.
       RETURN
  1    CONTINUE
       SP=SPIN(IL)
      BETA=BETAR(IEXP)
       TIME=TAU(IL)
       CALL GKK(IZ,BETA,SP,TIME,IL)
       VACDP(1,IL)=GKI(1)
       VACDP(2,IL)=GKI(2)
       VACDP(3,IL)=GKI(3)
        RETURN
        END
       SUBROUTINE GKK(IZ,BETA,SPIN,TIME,IL)
        COMMON/GVAC/GKI(3),SUM(3)
      COMMON/VAC/VACDP(3,75),QCEN,DQ,XNOR,AKS(6,75),IBYP
      COMMON/GGG/AVJI,GAMMA,XLAMB,TIMEC,GFAC,FIEL,POWER
      IF(IBYP.EQ.1)GO TO 100 
      IMEAN=0
      CALL XSTATIC(IZ,INQ,IFQ,BETA)
      L=0
      DO 77 I=1,6
  77  AKS(I,IL)=0.
 10   CONTINUE
      IF(IMEAN.EQ.1)INQ=1
      IF(IMEAN.EQ.1)IFQ=1
      DO 1 J=INQ,IFQ
      L=L+1
      NZ=IZ-J
      XJI=ATS(NZ)
      SM=SPIN
      IF(IMEAN.EQ.1)XJI=AVJI
      IF(SPIN.GT.XJI)SM=XJI
      NCOUP=INT(2.*SM+.5)+1
      SUM(1)=0.
      SUM(2)=0.
      SUM(3)=0.
      VALMI=SPIN-XJI
      IF(VALMI.LT.0.)VALMI=-VALMI
      DO 2 M=1,NCOUP
      F=VALMI+FLOAT(M)-1.
      DO 3 K=1,3
      RK=2.*FLOAT(K)
      IF2=F*2.+0.0001
      IRK2=RK*2.+0.0001
      ISPIN2=SPIN*2.+0.0001
      IXJI2=XJI*2.+0.0001
 3    SUM(K)=SUM(K)+((2.*F+1.)*WSIXJ(IF2,IF2,IRK2,ISPIN2,ISPIN2,IXJI2)
     *)**2/(2.*XJI+1.)
 2    CONTINUE
      IF(IMEAN.EQ.1)GO TO 11
      DO 4 K=1,3
      K1=2*K-1
 4    AKS(K1,IL)=AKS(K1,IL)
     *+SUM(K)*EXP(-((QCEN-FLOAT(J))/DQ)**2/2.)/XNOR
      IF(IMEAN.EQ.0)GO TO 1
 11   DO 12 K=1,3
      K1=2*K
 12   AKS(K1,IL)=AKS(K1,IL)+SUM(K)
 1    CONTINUE
      IMEAN=IMEAN+1
      IF(IMEAN.EQ.1)GO TO 10
 100  CONTINUE
      HMEAN=FIEL*IZ*(BETA**POWER)
      WSP=4789.*GFAC*HMEAN/AVJI
      WSP=WSP*TIMEC
      WSP=WSP*WSP*AVJI*(AVJI+1.)/3.
      DO 5 K=1,3
      K2=2*K
      K1=2*K-1
      WRT=WSP*K2*(K2+1)
      W2=WRT
      WRT=-WRT/(1.-AKS(K2,IL))
      XLAM=(1.-AKS(K2,IL))*(1.-EXP(WRT))/TIMEC
      UP=(GAMMA*TIME*AKS(K1,IL)+1.)/(TIME*GAMMA+1.)
      UP=UP*XLAMB*TIME+1.
      DOWN=TIME*(XLAM+XLAMB)+1.
      GKI(K)=UP/DOWN
      ALP=9.*XLAM*XLAM+8.*XLAM*TIMEC*(W2-XLAM*XLAM)
      ALP=SQRT(ALP)-3.*XLAM
      ALP=ALP/4./XLAM/TIMEC
      UPC=XLAM*TIME*(DOWN-2.*ALP*ALP*TIME*TIMEC)
      DWC=(DOWN+ALP*TIME)*(DOWN+2.*ALP*TIME)
      CCF=1.+UPC/DWC
      GKI(K)=GKI(K)*CCF
  5   CONTINUE
      RETURN
      END
      SUBROUTINE XSTATIC(IZ,IDO,IUP,BETA)
      COMMON/VAC/VACDP(3,75),QCEN,DQ,XNOR,AKS(6,75),IBYP
      H=1./(1.+(IZ**.45*.012008/BETA)**1.666667)
      QCEN=IZ*H**.6
      DQ=SQRT(QCEN*(1.-H))/2.
      IUP=INT(QCEN+3.*DQ+.5)
      IDO=INT(QCEN-3.*DQ-.5)
      IF(IUP.GT.IZ)IUP=IZ
      IF(IDO.LT.1)IDO=1
      XNOR=0.
      DO 1 LQ=IDO,IUP
  1   XNOR=XNOR+EXP(-((QCEN-FLOAT(LQ))/DQ)**2/2.)
      RETURN
      END
      FUNCTION ATS(N)
      IF(N.LE.0.OR.N.GT.96)GO TO 15
      X=N/2.+1
      M=N/2+1
      XM=FLOAT(M)
      IF(abs(x-xm).lt.1.e-9)GO TO 100
      GO TO(2,2,2,3,3,2,2,3,3,2,3,5,9,5,2,2,3,3,2,3,5,9,5
     *,2,2,3,3,2,7,7,2,5,13,17,2,3,2,9,5,2,2,3,3,2,
     *3,11,11,7),M
 100  M=M-1
      GO TO(1,1,1,4,1,1,1,4,1,16,1,6,8,1,1,1,4,1,1,1,6,8,1,
     *1,1,4,1,4,8,4,8,4,14,8,1,16,6,8,1,1,1,4,1,1,4,12,8,4),M
 1    ATS=0.
      RETURN
 2    ATS=.5
      RETURN
 3    ATS=1.5
      RETURN
 4    ATS=2.
      RETURN
 5    ATS=2.5
      RETURN
 6    ATS=3.
      RETURN
 7    ATS=3.5
      RETURN
   8  ATS=4.
      RETURN
 9    ATS=4.5
      RETURN
 11   ATS=5.5
      RETURN
 12   ATS=6.
      RETURN
 13   ATS=7.5
      RETURN
 14   ATS=8.
      RETURN
 15   ATS=0.
      RETURN
  16  ATS=1.
      RETURN
  17  ATS=6.5
      RETURN
      END
      SUBROUTINE YLM(THETA,YLMR)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      DIMENSION YLMR(9,9),ST(7)
      CT=COS(THETA)
      CTSQ=CT*CT
      IF(IAXS(IEXP).EQ.0)GO TO 120
      ST(1)=SIN(THETA)
      DO 10  I=2,7
      J=I-1
      ST(I)=ST(J)*ST(1)
  10  CONTINUE
      YLMR(1,3)=.1089659406
      YLMR(1,2)=-.2179318812*CT
      YLMR(1,1)=.0889703179*(3.*CTSQ-1.)
      YLMR(2,5)=.1248361677
      YLMR(2,4)=-.3530900028*CT
      YLMR(2,3)=.0943672726*(7.*CTSQ-1.)
      YLMR(2,2)=-.1334554768*CT*(7.*CTSQ-3.)
      YLMR(2,1)=.0298415518*((35.*CTSQ-30.)*CTSQ+3.)
      YLMR(3,7)=.1362755124
      YLMR(3,6)=-.4720722226*CT
      YLMR(3,5)=.100646136*(11.*CTSQ-1.)
      YLMR(3,4)=-.1837538634*CT*(11.*CTSQ-3.)
      YLMR(3,3)=.0918769316*((33.*CTSQ-18.)*CTSQ+1.)
      YLMR(3,2)=-.1162161475*CT*((33.*CTSQ-30.)*CTSQ+5.)
      YLMR(3,1)=.0179325408*(((231.*CTSQ-315.)*CTSQ
     1+105.)*CTSQ-5.)
      DO 100 L=1,3
      LF=2*L+1
      DO 99 M=2,LF
      YLMR(L,M)=YLMR(L,M)*ST(M-1)
  99  CONTINUE
 100  CONTINUE
      RETURN
 120  YLMR(1,1)=.0889703179*(3.*CTSQ-1.)
      YLMR(2,1)=.0298415518*((35.*CTSQ-30.)*CTSQ+3.)
      YLMR(3,1)=.0179325408*(((231.*CTSQ-315.)*CTSQ
     *+105.)*CTSQ-5.)
      RETURN
      END
      SUBROUTINE DECAY(CHISQ,NLIFT,CHILO)
      COMMON/TRA/DELTA(500,3),ENDEC(500),ITMA(50,200),ENZ(200)
      COMMON/GGG/AVJI,GAMMA,XLAMB,TIMEC,GFAC,FIEL,POWER
      COMMON/LIFE1/LIFCT(50),TIMEL(2,50)
      COMMON/VAC/VACDP(3,75),QCEN,DQ,XNOR,AKS(6,75),IBYP
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/CATLF/FP(4,500,3),GKP(4,500,2),KLEC(75)
      COMMON/LCDL/DELLA(500,3)
	common/switch/jzb,ibps
      DIMENSION GK(4)
      IDR=1
      DO 1 IL=1,NMAX1
      L=KSEQ(IDR,3)
      N1=28*(L-1)
      IBRA=KLEC(L)
      BSUM=0.
      IDRH=IDR
      DO 2 J=1,IBRA
      INX=KSEQ(IDR,1)
      INX1=KSEQ(IDR,2)
      EL1=0.
      IF(INX.NE.0)EL1=ELM(INX)
      EMT=EL1*EL1
      DELLA(IDR,1)=EMT
      IF(INX1.NE.0)EMT1=ELM(INX1)*ELM(INX1)
      BSUM=BSUM+DELTA(IDR,1)*EMT
      IF(INX1.EQ.0)GO TO 2
      DELLA(IDR,3)=EL1*ELM(INX1)
      DELLA(IDR,2)=EMT1
      BSUM=BSUM+DELTA(IDR,2)*EMT1
  2   IDR=IDR+1
      IDR=IDRH
      TAU(L)=1./BSUM
      CALL GKVAC(L)
      DO 20 J=1,IBRA
      L1=KSEQ(IDR,4)
      N2=28*(L1-1)
      INX1=KSEQ(IDR,2)
      DO 10 I=1,4
 10   GK(I)=GKP(I,IDR,1)*DELLA(IDR,1)
      IF(INX1.EQ.0)GO TO 11
      DO 12 I=1,4
 12   GK(I)=GK(I)+GKP(I,IDR,2)*DELLA(IDR,2)
 11   DO 13 I=1,4
      VCD=1.
      IF(I.NE.1)VCD=VACDP(I-1,L)
      GK(I)=GK(I)*TAU(L)
      IFN=2*I-1
      IU=(I-1)*7
      IF(IAXS(IEXP).EQ.0)IFN=1
      DO 22 KQ=1,IFN
      LC1=N1+IU+KQ
      LC2=N2+IU+KQ
 22   ZETA(LC2)=ZETA(LC2)+GK(I)*VCD*ZETA(LC1)
 13   CONTINUE
      IDR=IDR+1
 20   CONTINUE
 1    CONTINUE
      IBYP=1
      IF(NLIFT.EQ.0.OR.IEXP.NE.1)GO TO 50
      DO 51 JLT=1,NLIFT
      KL=LIFCT(JLT)
      DF=(TAU(KL)-TIMEL(1,JLT))/TIMEL(2,JLT)
      CHILO=CHILO+(ALOG(TAU(KL)/TIMEL(1,JLT))*TIMEL(1,JLT)/
     *TIMEL(2,JLT))**2
 51   CHISQ=CHISQ+DF*DF
 50   DO 30 L=2,NMAX
      IF(KLEC(L).EQ.0)GO TO 30
      N1=28*(L-1)
      DO 31 J=1,4
      VCD=1.
      IF(J.NE.1)VCD=VACDP(J-1,L)
      IFN=2*J-1
      IU=(J-1)*7
      DO 32 K=1,IFN
      LC1=N1+IU+K
 32   ZETA(LC1)=ZETA(LC1)*VCD
 31   CONTINUE
 30   CONTINUE
      RETURN
      END
        SUBROUTINE ANGULA(YGN,IDR,IFUL,FI0,FI1,TREC,GTH,FIGL,NGL,OP2)
	character*4 OP2  
      DIMENSION F(4),YLMR(9,9),AT(28),ALAB(9,9),ATTL(9,9)
     *,YGN(500)
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/TRA/DELTA(500,3),ENDEC(500),ITMA(50,200),ENZ(200)
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/CCC/eg(50),cc(50,5),NANG(200),Q(3,200,8),NICC,
     *AGELI(50,200,2)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/LCDL/DELLA(500,3)
      COMMON/CATLF/FP(4,500,3),GKP(4,500,2),KLEC(75)
      COMMON/BREC/BETAR(50)
      COMMON/THTAR/ITTE(50)
      COMMON/TCM/TETACM(50),TREP(50),DSIGS(50)
	COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      DO 1 L=1,IDR
      NLV=KSEQ(L,3)
      IL=(NLV-1)*28
      INX1=KSEQ(L,2)
      DO 2 J=1,4
      F(J)=FP(J,L,1)*DELLA(L,1)
 2    CONTINUE
      IF(INX1.EQ.0)GO TO 10
      DO 3 J=1,4
      F(J)=F(J)+2.*FP(J,L,3)*DELLA(L,3)
     *+FP(J,L,2)*DELLA(L,2)
 3    CONTINUE
 10   DO 11 J=1,4
      F(J)=F(J)*TAU(NLV)
      IU=(J-1)*7
      IFN=2*J-1
      IF(IAXS(IEXP).EQ.0)IFN=1
      DO 12 KQ=1,IFN
      IS=IU+KQ
      IG=IS+IL
 12   AT(IS)=ZETA(IG)*F(J)
 11   CONTINUE
      IF(IFUL.EQ.1)GO TO 100
      IXS=IAXS(IEXP)
      FI01=FI0-FIGL
      FI11=FI1-FIGL
      CALL FIINT(FI01,FI11,AT,IXS)
      IF(L.EQ.1)CALL YLM(GTH,YLMR)
      YGN(L)=AT(1)*.0795774715
      DO 14 JJ=1,3
      JI=JJ*7+1
      SM=YLMR(JJ,1)*AT(JI)
      IF(IAXS(IEXP).EQ.0)GO TO 16
      MIND=2*JJ+1
      DO 13 JM=2,MIND
      JI=JI+1
      SM=YLMR(JJ,JM)*AT(JI)*2.+SM
 13   CONTINUE
  16  IPD=ITMA(IEXP,NGL)
      ARG=(ENDEC(L)-ENZ(IPD))**2
      QV=(Q(3,IPD,2*JJ)*Q(2,IPD,2*JJ)+Q(1,IPD,2*JJ)*ARG)
     */(Q(2,IPD,2*JJ)+ARG)
      YGN(L)=YGN(L)+SM*QV
 14   CONTINUE
      GO TO 1
 100  DO 101 J=1,9
      DO 101 K=1,9
      ALAB(J,K)=0.
 101  ATTL(J,K)=0.
      DO 102 J=1,4
      LF=2*J-1
      LF1=LF
      IF(IAXS(IEXP).EQ.0)LF1=1
      DO 103 K=1,LF1
      INAT=(J-1)*7+K
 103  ALAB(LF,K)=AT(INAT)
 102  CONTINUE
      BT=BETAR(IEXP)
      IF(ITTE(IEXP).EQ.1)GO TO 200
      CALL RECOIL(ALAB,ATTL,BT,TREC)
 200  CONTINUE
      IF(L.EQ.1)CALL YLM1(GTH,YLMR)
      IXS=IAXS(IEXP)
      FI01=FI0-FIGL
      FI11=FI1-FIGL
      CALL FIINT1(FI01,FI11,ALAB,IXS)
      YGN(L)=ALAB(1,1)*.0795774715
      DO 104 J=2,9
      SM=YLMR(J,1)*ALAB(J,1)
      IF(IAXS(IEXP).EQ.0)GO TO 105
      DO 106 K=2,J
 106  SM=SM+2.*YLMR(J,K)*ALAB(J,K)
  105 IPD=ITMA(IEXP,NGL)
      ARG=(ENDEC(L)-ENZ(IPD))**2
      QV=(Q(3,IPD,J-1)*Q(2,IPD,J-1)+Q(1,IPD,J-1)*ARG)
     */(Q(2,IPD,J-1)+ARG)
      YGN(L)=YGN(L)+SM*QV
 104  CONTINUE
 1    CONTINUE
      if(OP2.eq.'INTG')return
      dsig=dsigs(iexp)
	ttx=tlbdg(iexp)/57.2957795
	do 300 j=1,idr
 300  ygn(j)=ygn(j)*dsig*sin(ttx)	
      RETURN
      END
      SUBROUTINE READY(IDR,NTAP,IPRI)
      DIMENSION IYTOT(32)
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/CCCDS/ndst(50)
      REWIND NTAP
      DO 33 K=1,LP6
 33   IYTOT(K)=0
      IF(IPRI.EQ.1)WRITE(22,99)
      DO 200 LXP=1,NEXPT
      DO 501 KKL=1,LP6
 501  NYLDE(LXP,KKL)=0
      II=ndst(LXP)
      DO 500 KK=1,II
      READ(NTAP,*)NE,NANX,ZP,AP,XEP,NVAL,WAGA
      IF(IPRI.EQ.1)WRITE(22,1)NE,ZP,AP,XEP,ndst(NE),WAGA
 1    FORMAT(1X,///5X,10HEXPERIMENT,1X,1I2/
     *2X,10HPROJECTILE,1X,1H(,1F4.0,1H,,1F4.0,1H),1X,
     *1F7.3,1X,3HMEV,1X,3H---,1I1,1X,
     *18HGE(LI) DETECTOR(S),2X,7HWEIGHT=,1E8.2/20X,11H** EXPERIME
     *,14HNTAL YIELDS **)
      IF(IPRI.EQ.1)WRITE(22,888)
      DO 3 J=1,NVAL
      READ(NTAP,*)NS1,NS2,U,W
      NSXH=NS1
      NSYH=NS2
      IF(NS1.LT.100)GO TO 450
      NS1=NS1/100
      NS2=NS2/100
 450  DO 15 NDE=1,IDR
      IF(NS1.EQ.KSEQ(NDE,3).AND.NS2.EQ.KSEQ(NDE,4))GO TO 16
 15   CONTINUE
      IF(IPRI.EQ.1)WRITE(22,17)NS1,NS2
 17    FORMAT(1X///5X,38HERROR-NO MATRIX ELEMENT BETWEEN STATES
     *,1X,1I2,5H AND ,1I2,/10X,23HTHIS TRANSITION IGNORED,//)
      GO TO 3
 16   IDC=NDE
      IYTOT(KK)=IYTOT(KK)+1
      IDC1=0
      IF(NSXH.LT.100)GO TO 460
      NS3=NSXH-100*NS1
      NS4=NSYH-100*NS2
      DO 461 NDE1=1,IDR
      IF(NS3.EQ.KSEQ(NDE1,3).AND.NS4.EQ.KSEQ(NDE1,4))GO TO 462
 461  CONTINUE
      IF(IPRI.EQ.1)WRITE(22,17)NS3,NS4
      GO TO 463
 462  IDCX=IDC*1000+NDE1
      IF(IDC.GT.NDE1)IDCX=NDE1*1000+IDC
      IDC=IDCX
 463  CONTINUE
 460  CONTINUE
      IDC1=IDC
      IF(IDC1.GT.1000)IDC1=IDC/1000
      IF(IPRI.EQ.1)WRITE(22,4)IDC,KSEQ(IDC1,3),KSEQ(IDC1,4),
     *U,W
 4    FORMAT(2X,1I6,2X,1I2,2X,1I2,1(1E14.6,3X,1E14.6))
 888  FORMAT(4X,5HDECAY,1X,2HIS,2X,2HIF,1(9X,13HYIELD+/-ERROR,9X)/)
      IYTT=IYTOT(KK)
      YEXP(KK,IYTT)=U
      DYEX(KK,IYTT)=W/(SQRT(WAGA)+1.E-4)
      IY(IYTT,KK)=IDC
 3    CONTINUE
      IYTT=IYTOT(KK)
      LBG=IYTT-NVAL+1
      CALL SZEREG(LBG,IYTT,KK)
 500  NYLDE(LXP,KK)=NVAL
 200  CONTINUE
 99   FORMAT(5X/47X,
     *41HREPRINT OF EXPERIMENTAL DATA TO BE FITTED//)
      RETURN
      END
      SUBROUTINE BRANR(CHISQ,NWYR,CHILO)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/BRNCH/BRAT(50,2),IBRC(2,50),NBRA
      COMMON/TRA/DELTA(500,3),ENDEC(500),ITMA(50,200),ENZ(200)
      COMMON/PRT/IPRM(20)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/LEV/TAU(75),KSEQ(500,4)
      IF(NBRA.EQ.0)RETURN
      IF(IPRM(3).EQ.-1)WRITE(22,11)
      NWYR=NWYR+NBRA
      MUL2=MULTI(1)+MULTI(2)
      DO 1 K=1,NBRA
      CH1=0.
      CH2=0.
      IFLG=1
      ITT=1
      IOUT=0
      N1=IBRC(1,K)
      N2=IBRC(2,K)
      I1=KSEQ(N1,1)
      I2=KSEQ(N2,1)
      ENG1=EN(KSEQ(N1,3))-EN(KSEQ(N1,4))
      ENG2=EN(KSEQ(N2,3))-EN(KSEQ(N2,4))
      IF(I1.EQ.0)GO TO 717
      IF(I1.LE.MULTI(1))LAB1=1
      IF(I1.GT.MULTI(1).AND.I1.LE.MUL2)LAB1=2
      IF(I1.GT.MUL2)LAB1=3
  717 CONTINUE
      IF(I2.EQ.0)GO TO 718
      IF(I2.LE.MULTI(1))LAB2=1
      IF(I2.GT.MULTI(1).AND.I2.LE.MUL2)LAB2=2
      IF(I2.GT.MUL2)LAB2=3
  718 CONTINUE
      IF(I1.EQ.0)GO TO 712
      CH1=ELM(I1)*ELM(I1)*DELTA(N1,1)/(1.+CONV(ENG1,LAB1))
  712 IF(I2.EQ.0)GO TO 713
      CH2=ELM(I2)*ELM(I2)*DELTA(N2,1)/(1.+CONV(ENG2,LAB2))
  713 J1=KSEQ(N1,2)
      IF(J1.EQ.0)GO TO 2
      IFLG=IFLG+1
      LAB1=LAB1+2
      CH1=CH1+ELM(J1)*ELM(J1)*DELTA(N1,2)/(1.+CONV(ENG1,LAB1))
 2    J2=KSEQ(N2,2)
      IF(J2.EQ.0)GO TO 4
      IFLG=IFLG+1
      LAB2=LAB2+2
      CH2=CH2+ELM(J2)*ELM(J2)*DELTA(N2,2)/(1.+CONV(ENG2,LAB2))
 4    U=(CH1/CH2-BRAT(K,1))/BRAT(K,2)
      CHISQ=CHISQ+U*U
      CHILO=CHILO+(BRAT(K,1)*ALOG(CH1/CH2/BRAT(K,1))
     */BRAT(K,2))**2
      IF(IPRM(3).EQ.-1)WRITE(22,10)KSEQ(N1,3),KSEQ(N1,4),
     *KSEQ(N2,3),KSEQ(N2,4),BRAT(K,1),BRAT(K,2),CH1/CH2
     *,-U
 1    CONTINUE
      IF(IPRM(3).EQ.-1)IPRM(3)=0
      RETURN
 11   FORMAT(1X,///10X,36HEXP. AND CALCULATED BRANCHING RATIOS
     *,//5X,3HNS1,5X,3HNF1,5X,3HNS2,5X,3HNF2,5X,10HRATIO(1:2),
     *9X,5HERROR,7X,10HCALC.RATIO,5X,15H(EXP-CAL)/ERROR,//)
 10   FORMAT(5X,3(1I2,6X),1I2,5X,3(1F10.5,5X),5X,1F4.1)
      END
      SUBROUTINE LIMITS
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      DO 1 J=1,MEMAX
      IF(IVAR(J).EQ.0)GO TO 1
      IF(ELM(J).LE.ELMU(J).AND.ELM(J).GE.ELML(J))GO TO 1
      IF(ELM(J).GT.ELMU(J))GO TO 2
      ELM(J)=ELML(J)
      WRITE(22,10)J,ELM(J)
      GO TO 1
 2    ELM(J)=ELMU(J)
      WRITE(22,10)J,ELM(J)
 1    CONTINUE
 10   FORMAT(2X,'Warning - matrix element ',1I3,' reset to ',1F10.6)
      RETURN
      END
      SUBROUTINE SZEREG(LST,LS,L)
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      IF(LST.EQ.LS)RETURN
      LST1=LST
      LSP=LS-1
 100  IA=IY(LST1,L)
      IF(IA.GT.1000)IA=IA/1000
      INX=LST1
      DO 1 K=LST1,LSP
      IB=IY(K+1,L)
      IF(IB.GT.1000)IB=IB/1000
      IA=MIN0(IA,IB)
      IF(IA.EQ.IB)INX=K+1
 1    CONTINUE
      IF(INX.EQ.LST1)GO TO 200
      IH=IY(LST1,L)
      IY(LST1,L)=IY(INX,L)
      IY(INX,L)=IH
      YH=YEXP(L,LST1)
      DYH=DYEX(L,LST1)
      YEXP(L,LST1)=YEXP(L,INX)
      DYEX(L,LST1)=DYEX(L,INX)
      YEXP(L,INX)=YH
      DYEX(L,INX)=DYH
 200  LST1=LST1+1
      IF(LST1.GT.LSP)RETURN
      GO TO 100
      END
      SUBROUTINE SIXEL(RIK,RV,EM,JK,KK,INDX,LU)
      COMPLEX ARM
      COMMON/AZ/ARM(600,7)
      COMMON/ODCH/DEV(500)
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/TRB/ITS
      COMMON/SEL/KVAR(500)
      KK6=KK+5
      RN=DEV(LU)
      AL=(RV-RN)*20./RIK
      IF(ITS.EQ.1.AND.KVAR(INDX).NE.0)WRITE(18,*)LU,INDX,IEXP,AL/EM
      AL1=ABS(AL)
      IF(ITS.EQ.2)WRITE(18,*)LU,INDX,IEXP,AL1
      IF(AL1.LE.ABS(AIMAG(ARM(KK6,JK))))RETURN
      DO 1 J=KK,KK6
      A1=ABS(AIMAG(ARM(J,JK)))
      IF(AL1.LE.A1)GO TO 1
      J1=J+1
      DO 2 L=J1,KK6
      L1=KK6+J1-L
      C1=REAL(ARM(L1-1,JK))
      C2=AIMAG(ARM(L1-1,JK))
 2    ARM(L1,JK)=CMPLX(C1,C2)
      RX=FLOAT(INDX)
      ARM(J,JK)=CMPLX(RX,AL)
      GO TO 11
 1    CONTINUE
 11   RETURN
      END
      SUBROUTINE PRELM(IOP)
      character*3 wrn
      COMMON/HHH/HLM(500)
	common/switch/jzb,ibps
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/COEX2/NMAX,NDIM,NMAX1
	common/resc/elm25(500),elm26(500)
      INX=0
      WRITE(22,100)
      DO 1 J=1,8
      M=MULTI(J)
      IF(M.EQ.0)GO TO 1
      WRITE(22,2)J
      IF(IOP.EQ.1)WRITE(22,101)
      IF(IOP.EQ.2.or.iop.eq.4)WRITE(22,102)
      IF(IOP.EQ.3)WRITE(22,103)
      DO 3 K=1,NMAX
      L=LDNUM(J,K)
      IF(L.EQ.0)GO TO 3
      DO 4 KK=1,L
      INX=INX+1
      GO TO (5,6,7,77)IOP
 5    WRITE(22,10)INX,LEAD(1,INX),LEAD(2,INX),ELM(INX)
      GO TO 4
 6    CONTINUE
 77   continue
      IF(IVAR(INX).EQ.0)GO TO 20
      IF(IVAR(INX).GT.1000)GO TO 21
      if(iop.eq.2)WRITE(22,11)INX,LEAD(1,INX),LEAD(2,INX)
     *,ELM(INX),ELML(INX)
     *,ELMU(INX)
      if(iop.eq.4.and.jzb.eq.25)
     *WRITE(22,11)INX,LEAD(1,INX),LEAD(2,INX)
     *,ELM25(INX),ELML(INX)
     *,ELMU(INX)
       if(iop.eq.4.and.jzb.eq.26)
     * WRITE(22,11)INX,LEAD(1,INX),LEAD(2,INX)
     *,ELM26(INX),ELML(INX)
     *,ELMU(INX)
      GO TO 4
 21   WRITE(22,12)INX,LEAD(1,INX),LEAD(2,INX),ELM(INX),(IVAR(INX)-1000)
      GO TO 4
 20   WRITE(22,13)INX,LEAD(1,INX),LEAD(2,INX),ELM(INX)
      GO TO 4
 7    ISP=LEAD(2,INX)
       PV=(ELMU(INX)-ELML(INX))/100.
       WRN='   '
       IF((ELM(INX)-ELML(INX)).LT.PV)WRN= '*?*'
       IF((ELMU(INX)-ELM(INX)).LT.PV)WRN= '*?*'
      if(jzb.eq.25)STE=elm25(INX)
	if(jzb.eq.26)ste=elm26(inx)
      B=ELM(INX)*ELM(INX)/(2.*SPIN(ISP)+1.)
      IF(LEAD(1,INX).EQ.LEAD(2,INX))B=9999999.
      WRITE(22,11)INX,LEAD(1,INX),LEAD(2,INX),ELM(INX),
     *100.*(ELM(INX)-STE)/STE,B,WRN
 4    CONTINUE
 3    CONTINUE
 1    CONTINUE
 100  FORMAT(2X/40X,15HMATRIX ELEMENTS,//)
 2    FORMAT(5X,14HMULTIPOLARITY=,1I1)
 10   FORMAT(5X,1I3,5X,1I2,5X,1I2,5X,1F10.5)
 11   FORMAT(5X,1I3,5X,1I2,5X,1I2,3(5X,1F10.5),1A3)
 12   FORMAT(5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,
     *10HCOUPLED TO,1X,1I3)
 13   FORMAT(5X,1I3,5X,1I2,5X,1I2,5X,1F10.5,5X,5HFIXED)
 101  FORMAT(4X,5HINDEX,3X,2HNF,5X,2HNS,10X,2HME)
 102  FORMAT(4X,5HINDEX,3X,2HNF,5X,2HNS,10X,2HME,15X,6HLIMITS)
 103  FORMAT(4X,5HINDEX,3X,2HNF,5X,2HNS,10X,2HME,10X,
     *9HPC CHANGE,5X,17HRED. TRANS. PROB.)
      RETURN
      END
        SUBROUTINE RECOIL(ALAB,ATTL,BETA,THETA)
        DIMENSION ALAB(9,9),ATTL(9,9),ATEMP(16)
        HOLD=ALAB(1,1)
      IF(abs(HOLD).lt.1.e-9)RETURN
        CALL ROTATE(ALAB,ATTL,-THETA,7,2)
        ATTL(2,1)=(2./SQRT(15.))*(SQRT(5.)*ATTL(1,1)
     *    -ATTL(3,1))
        ATTL(2,2)=-ATTL(3,2)/SQRT(5.)
        ATTL(4,1)=(4./SQRT(35.))*(3.*ATTL(3,1)
     *    -SQRT(5.)*ATTL(5,1))
        ATTL(4,2)=(8.*SQRT(2.)*ATTL(3,2)
     *    -5.*SQRT(3.)*ATTL(5,2))/SQRT(35.)
        ATTL(4,3)=(2./SQRT(7.))*(2.*ATTL(3,3)
     *    -SQRT(3.)*ATTL(5,3))
        ATTL(4,4)=-ATTL(5,4)
        ATTL(6,1)=(10./SQRT(11.))*(ATTL(5,1)-(3./SQRT(13.))*ATTL(7,1))
        ATTL(6,2)=(1./SQRT(11.))*(4.*SQRT(6.)*ATTL(5,2)-5.*
     *  SQRT(35./13.)*ATTL(7,2))
        ATTL(6,3)=SQRT(4./11.)*(SQRT(21.)*ATTL(5,3)-10.*SQRT(2./13.)*
     *  ATTL(7,3))
        ATTL(6,4)=SQRT(1./11.)*(8.*ATTL(5,4)-15.*SQRT(3./13.)*
     *  ATTL(7,4))
        ATTL(6,5)=SQRT(4./11.)*(3.*ATTL(5,5)-5.*SQRT(5./13.)*
     *  ATTL(7,5))
        ATTL(6,6)=-ATTL(7,6)*SQRT(25./13.)
        ATTL(8,1)=(56./SQRT(195.))*ATTL(7,1)
        ATTL(8,2)=(32./SQRT(65.))*ATTL(7,2)
        ATTL(8,3)=(8.*SQRT(3./13.))*ATTL(7,3)
        ATTL(8,4)=(16.*SQRT(2./39.))*ATTL(7,4)
        ATTL(8,5)=(8.*SQRT(11./65.))*ATTL(7,5)
        ATTL(8,6)=(16.*SQRT(2./65.))*ATTL(7,6)
        ATTL(8,7)=(8./SQRT(15.))*ATTL(7,7)
        DO 150  L=2,8,2
        DO 149  M=1,L
        ATTL(L,M)=BETA * ATTL(L,M)
  149   CONTINUE
  150   CONTINUE
        BETASQ=BETA*BETA
        IF (BETASQ .LT. 1.0E-10)  GO TO 199
        I1=0
        DO 400 I=1,7,2
        DO 400 J=1,I
        I1=I1+1
        ATEMP(I1)=ATTL(I,J)
 400    CONTINUE
        DUM=(2./5.)*SQRT(5.)*ATEMP(1)-(10./7.)*ATEMP(2)+
     *  (12./35.)*SQRT(5.)*ATEMP(5)
        ATTL(3,1)=ATEMP(2)+BETASQ*DUM
        DUM=-(17./14.)*ATEMP(3)+(2./7.)*SQRT(6.)*ATEMP(6)
        ATTL(3,2)=ATEMP(3)+BETASQ*DUM
        DUM=-(4./7.)*ATEMP(4)+(2./7.)*SQRT(3.)*ATEMP(7)
        ATTL(3,3)=ATEMP(4)+BETASQ*DUM
        DUM=(8./7.)*SQRT(5.)*ATEMP(2)-(380./77.)*ATEMP(5)
     *  +(100./11.)*SQRT(1./13.)*ATEMP(10)
        ATTL(5,1)=ATEMP(5)+BETASQ*DUM
        DUM=(20./21.)*SQRT(6.)*ATEMP(3)-(723./154.)*ATEMP(6)+
     *  (20./11.)*SQRT(70./39.)*ATEMP(11)
        ATTL(5,2)=ATEMP(6)+BETASQ*DUM
        DUM=(20./21.)*SQRT(3.)*ATEMP(4)-(306./77.)*ATEMP(7)+
     *  (40./11.)*SQRT(14./39.)*ATEMP(12)
        ATTL(5,3)=ATEMP(7)+BETASQ*DUM
        DUM=-(61./22.)*ATEMP(8)+(40./11.)*SQRT(3./13.)*ATEMP(13)
        ATTL(5,4)=ATEMP(8)+BETASQ*DUM
        DUM=-(12./11.)*ATEMP(9)+(20./11.)*SQRT(5./13.)*ATEMP(14)
        ATTL(5,5)=ATEMP(9)+BETASQ*DUM
        DUM=(210./11.)*SQRT(1./13.)*ATEMP(5)-(574./55.)*ATEMP(10)
        ATTL(7,1)=ATEMP(10)+BETASQ*DUM
        DUM=(14./11.)*SQRT(210./13.)*ATEMP(6)-(1121./110.)*ATEMP(11)
        ATTL(7,2)=ATEMP(11)+BETASQ*DUM
        DUM=(28./11.)*SQRT(42./13.)*ATEMP(7)-(104./11.)*ATEMP(12)
        ATTL(7,3)=ATEMP(12)+BETASQ*DUM
        DUM=(84./11.)*SQRT(3./13.)*ATEMP(8)-(181./22.)*ATEMP(13)
        ATTL(7,4)=ATEMP(13)+BETASQ*DUM
        DUM=(42./11.)*SQRT(5./13.)*ATEMP(9)-(358./55.)*ATEMP(14)
        ATTL(7,5)=ATEMP(14)+BETASQ*DUM
        ATTL(7,6)=ATEMP(15)*(1.-(43./10.)*BETASQ)
        ATTL(7,7)=ATEMP(16)*(1.-(8./5.)*BETASQ)
        ATTL(9,1)=(672./5.)*SQRT(1./221.)*ATEMP(10)*BETASQ
        ATTL(9,2)=(144./5.)*SQRT(21./221.)*ATEMP(11)*BETASQ
        ATTL(9,3)=36.*SQRT(12./221.)*ATEMP(12)*BETASQ
        ATTL(9,4)=24.*SQRT(22./221.)*ATEMP(13)*BETASQ
        ATTL(9,5)=(144./5.)*SQRT(11./221.)*ATEMP(14)*BETASQ
        ATTL(9,6)=(72./5.)*SQRT(2./17.)*ATEMP(15)*BETASQ
        ATTL(9,7)=(24./5.)*SQRT(7./17.)*ATEMP(16)*BETASQ
  199   CALL ROTATE(ATTL,ALAB,THETA,9,1)
        TEST=ABS(1.0-ALAB(1,1)/HOLD)
        IF(TEST.GT.1.0E-07)GO TO 333
        GOTO 51
  333   WRITE(22,335)TEST
  51    CONTINUE
        RETURN
  335   FORMAT(18H ERROR IN ROTATION,1X,1E10.3/)
        END
        SUBROUTINE ROTATE(ALAB,ATTL,THETA,K2,KD)
        DIMENSION ALAB(9,9),ATTL(9,9)
        IF(ABS(THETA).GT..01)GO TO 1
        DO 2 J=1,9
        DO 2 K=1,9
  2     ATTL(J,K)=ALAB(J,K)
        RETURN
   1    CONTINUE
        DJARG=THETA
        DO 111 KA=1,K2,KD
        IDJ=KA-1
        DO 110 KAPPA=1,KA
        IDMP=KAPPA-1
        SUM=0.0
        DO 105  KAPRI=1,KA
        IDM=KAPRI-1
        DKKK=DJMM(DJARG,IDJ,IDM,IDMP)
        SUM=SUM+DKKK*ALAB(KA,KAPRI)
  105   CONTINUE
        IF (KA .EQ. 1)  GO TO 100
        DO 106  KAPRI=2,KA
        IDM=-KAPRI+1
        DKKK=DJMM(DJARG,IDJ,IDM,IDMP)
        SUM=SUM+DKKK*ALAB(KA,KAPRI)*(-1.0)**(KAPRI-1)
  106   CONTINUE
  100   ATTL(KA,KAPPA)=SUM
  110   CONTINUE
  111   CONTINUE
        RETURN
        END
        SUBROUTINE YLM1(THETA,YLMR)
        DIMENSION YLMR(9,9),ST(9)
        CT=COS(THETA)
        CTSQ=CT*CT
        ST(1)=SIN(THETA)
        DO 10  I=2,9
        J=I-1
        ST(I)=ST(J)*ST(1)
  10    CONTINUE
        DO 30 L=2,9
        DO 29 M=1,9
        YLMR(L,M)=0.0
  29    CONTINUE
  30    CONTINUE
        YLMR(2,2)=-SQRT(6.)/2.
        YLMR(2,1)=SQRT(3.)*CT
        YLMR(3,3)=SQRT(30.)/4.
        YLMR(3,2)=-(SQRT(30.)/2.)*CT
        YLMR(3,1)=(SQRT(5.)/2.)*(3.*CTSQ-1.)
        YLMR(4,4)=-SQRT(35.)/4.
        YLMR(4,3)=(SQRT(210.)/4.)*CT
        YLMR(4,2)=-(SQRT(21.)/4.)*(5.*CTSQ-1.)
        YLMR(4,1)=(SQRT(7.)/2.)*CT*(5.*CTSQ-3.)
        YLMR(5,5)=3.*SQRT(70.)/16.
        YLMR(5,4)=-(3.*SQRT(35.)/4.)*CT
        YLMR(5,3)=(3.*SQRT(10.)/8.)*(7.*CTSQ-1.)
        YLMR(5,2)=-(3.*SQRT(5.)/4.)*CT*(7.*CTSQ-3.)
        YLMR(5,1)=(3./8.)*((35.*CTSQ-30.)*CTSQ+3.)
        YLMR(6,6)=-3.*SQRT(77.)/16.
        YLMR(6,5)=(3.*SQRT(770.)/16.)*CT
        YLMR(6,4)=-(SQRT(385.)/16.)*(9.*CTSQ-1.)
        YLMR(6,3)=(SQRT(2310.)/8.)*CT*(3.*CTSQ-1.)
        YLMR(6,2)=-(SQRT(330.)/16.)*((21.*CTSQ-14.)*CTSQ+1.)
        YLMR(6,1)=(SQRT(11.)/8.)*CT*((63.*CTSQ-70.)*CTSQ+15.)
        YLMR(7,7)=SQRT(3003.)/32.
        YLMR(7,6)=-(3.*SQRT(1001.)/16.)*CT
        YLMR(7,5)=(3.*SQRT(182.)/32.)*(11.*CTSQ-1.)
        YLMR(7,4)=-(SQRT(1365.)/16.)*CT*(11.*CTSQ-3.)
        YLMR(7,3)=(SQRT(1365.)/32.)*((33.*CTSQ-18.)*CTSQ+1.)
        YLMR(7,2)=-(SQRT(546.)/16.)*CT*((33.*CTSQ-30.)*CTSQ+5.)
        YLMR(7,1)=(SQRT(13.)/16.)*(((231.*CTSQ-315.)*CTSQ
     1  +105.)*CTSQ-5.)
        YLMR(8,8)=-3.*SQRT(1430.)/64.
        YLMR(8,7)=(3.*SQRT(5005.)/32.)*CT
        YLMR(8,6)=-(3.*SQRT(770.)/64.)*(13.*CTSQ-1.)
        YLMR(8,5)=(3.*SQRT(770.)/32.)*(13.*CTSQ-3.)*CT
        YLMR(8,4)=-(3.*SQRT(70.)/64.)*((143.*CTSQ-66.)*CTSQ+3.)
        YLMR(8,3)=(3.*SQRT(35.)/32.)*((143.*CTSQ-110.)*CTSQ+15.)*CT
        YLMR(8,2)=-(SQRT(210.)/64.)*(((429.*CTSQ-495.)*CTSQ+135.)*
     1  CTSQ-5.)
        YLMR(8,1)=(SQRT(15.)/16.)*(((429.*CTSQ-693.)*CTSQ+315.)*
     1  CTSQ-35.)*CT
        YLMR(9,9)=3.*SQRT(24310.)/256.
        YLMR(9,8)=-(3.*SQRT(24310.)/64.)*CT
        YLMR(9,7)=(SQRT(7293.)/64.)*(15.*CTSQ-1.)
        YLMR(9,6)=-(3.*SQRT(34034.)/64.)*(5.*CTSQ-1.)*CT
        YLMR(9,5)=(3.*SQRT(2618.)/128.)*((65.*CTSQ-26.)*CTSQ+1.)
        YLMR(9,4)=-(SQRT(39270.)/64.)*((39.*CTSQ-26.)*CTSQ+3.)*CT
        YLMR(9,3)=(3.*SQRT(595.)/64.)*(((143.*CTSQ-143.)*CTSQ+33.)*
     1  CTSQ-1.)
        YLMR(9,2)=-(3.*SQRT(34.)/64.)*(((715.*CTSQ-1001.)*CTSQ+385.)
     1  *CTSQ-35.)*CT
        YLMR(9,1)=(SQRT(17.)/128.)*((((6435.*CTSQ-12012.)*CTSQ+6930.)
     1  *CTSQ-1260.)*CTSQ+35.)
        DO 100 L=2,9
        YLMR(L,1)=YLMR(L,1)*.0795774715
        DO 99 M=2,L
        YLMR(L,M)=YLMR(L,M)*ST(M-1)*.0795774715
99      CONTINUE
100     CONTINUE
        RETURN
        END
      SUBROUTINE FIINT(FI0,FI1,AT,IXS)
      DIMENSION AT(28)
      IF(IXS.EQ.0)GO TO 1
      DO 2 M=2,7
      JS=M/2
      MM=M-1
      WSP=(SIN(MM*FI1)-SIN(MM*FI0))/MM
      JS=JS*7+M
      JF=M+21
      DO 3 J=JS,JF,7
 3    AT(J)=AT(J)*WSP
 2    CONTINUE
      WSP=FI1-FI0
 1    IF(IXS.EQ.0)WSP=6.283185308
      DO 4 J=1,4
      JS=(J-1)*7+1
 4    AT(JS)=AT(JS)*WSP
      RETURN
      END
      SUBROUTINE FIINT1(FI0,FI1,ALAB,IXS)
      DIMENSION ALAB(9,9)
      IF(IXS.EQ.0)GO TO 1
      DO 2 M=2,9
      MM=M-1
      WSP=(SIN(MM*FI1)-SIN(MM*FI0))/MM
      DO 3 J=1,9
  3   ALAB(J,M)=ALAB(J,M)*WSP
 2    CONTINUE
      WSP=FI1-FI0
 1    IF(IXS.EQ.0)WSP=6.283185308
      DO 4 J=1,9
 4    ALAB(J,1)=ALAB(J,1)*WSP
      RETURN
      END
      SUBROUTINE TAPMA(LX,ISKE,ISKO,ISKF,NFLR,IDR,NCO,NFT,ENB)
      COMMON/VLIN/XV(51),YV(51),ZV(20),DSG(20),DSE(20),DS
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/YTEOR/YGN(500),YGP(500),IFMO
      NFT=0
      NFILT=0
      REWIND 14
      IF(ISKE.EQ.0)GO TO 2
 1    READ(14,*)NE,NTT,EMN,EMX,TMN,TMX,NA,TMX,TMX,TMX
      NFIL=NE*NTT*NA
      NFILT=NFILT+NFIL
      DO 3 J=1,NFIL
 3    READ(14,*)LX1,ENB,TTA,NG,DS,(YGN(K),K=1,IDR)
      IF(NFILT.EQ.ISKE)GO TO 2
      GO TO 1
 2    IF(NCO.EQ.0)RETURN
      READ(14,*)NE,NTT,EMN,EMX,TMN,TMX,NA,TMX,TMX,TMX
      IF(ISKO.EQ.0)GO TO 5
      DO 4 J=1,ISKO
 4    READ(14,*)LX1,ENB,TTA,NG,DS,(YGN(K),K=1,IDR)
 5    CONTINUE
      DO 6 J=1,NFLR
      JS=(J-1)*IDR+1
      JF=JS+IDR-1
      READ(14,*)LX1,ENB,TTA,NG1,DS,(ZETA(K),K=JS,JF)
      IF(LX1.NE.LX)NFT=1
      IF(NFT.EQ.1)GO TO 99
      XV(J)=TTA/57.2957795
      IF(ISKF.EQ.0.OR.J.EQ.NFLR)GO TO 6
      DO 7 JJ=1,ISKF
 7    READ(14,*)LX1,EN0,TTA,NG,DS,(YGN(K),K=1,IDR)
 6    CONTINUE
      RETURN
 99   WRITE(22,98)
 98   FORMAT(10X///10X,15HTAPE READ ERROR/
     *10X,11HJOB ABORTED)
      RETURN
      END
      FUNCTION SIMIN(NP,H,Y)
      DIMENSION Y(101)
      IF(NP.LT.3)GO TO 2
      IK=NP-2
      SM=Y(1)+Y(NP)
      DO 1 IN=1,IK
      EE=IN/2.
 1    SM=SM+2.*Y(IN+1)/(1.+INT(EE)-EE)
      SIMIN=SM*H/3.
      RETURN
 2    IF(NP.EQ.1)GO TO 3
      SIMIN=(Y(1)+Y(2))*H/2.
      RETURN
 3    SIMIN=Y(1)
      RETURN
      END
      SUBROUTINE MIXUP
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/XRA/SE
      DO 2 K=1,MEMAX
      IF(IVAR(K).EQ.0.OR.IVAR(K).GT.999)GO TO 2
      ELM(K)=ELML(K)+RNDM(SE)*(ELMU(K)-ELML(K))
 2    CONTINUE
      DO 3 K=1,MEMAX
      IF(IVAR(K).LT.999)GO TO 3
      K1=IVAR(K)-1000
      IF(abs(ELMU(K1)).lt.1.e-9)GO TO 4
      ELM(K)=ELM(K1)*SA(K)
      GO TO 3
 4    ELM(K)=0.
 3    CONTINUE
      RETURN
      END
      FUNCTION FXIS1(I,N)
      COMMON/CXI/XI(500)
      GO TO (1,2,2,1,2,2,1)N
 1    FXIS1=-SIGN(1.,XI(I))
      RETURN
 2    FXIS1=1.
      RETURN
      END
      FUNCTION FXIS2(I,N)
      COMMON/CXI/XI(500)
      GO TO(1,2,2,1,2,2,1)N
 1    FXIS2=1.
      RETURN
 2    FXIS2=-SIGN(1.,XI(I))
      RETURN
      END
      SUBROUTINE PODZIEL(I,J)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/APRCAT/QAPR(500,2,7),IAPR(500,2),ISEX(75)
      COMMON/APRX/LERF,IDIVE(50,2)
      IF(I.EQ.3)GO TO 7
      IF(I.EQ.1)GO TO 2
      L1=IDIVE(J,2)
      IDIVE(J,2)=L1+1
 7    L2=IDIVE(J,2)
      IF(I.EQ.3)L1=1
      DO 5 K=1,LP2
      DO 5 L=1,7
 5    QAPR(K,2,L)=QAPR(K,2,L)*L1/L2
      IF(I.EQ.2)WRITE(22,11)J,IDIVE(J,1),L2
      IF(I.EQ.3)GO TO 8
      RETURN
 2    L1=IDIVE(J,1)
      IDIVE(J,1)=L1+1
 8    L2=IDIVE(J,1)
      IF(I.EQ.3)L1=1
      DO 6 K=1,LP2
      DO 6 L=1,7
 6    QAPR(K,1,L)=QAPR(K,1,L)*L1/L2
      IF(I.EQ.1)WRITE(22,11)J,L2,IDIVE(J,2)
      RETURN
 11   FORMAT(5X,5H*****,1X,25HEXP(A) EXPANSION FAILURE!,1X,5H*****
     */5X,10HEXPERIMENT,1X,1I2,3X,15HNEW SUBDIVISION,1X,
     *1H(,1I1,1H,,1I1,1H))
      END
      SUBROUTINE KLOPOT(K,RLR)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,LP8,
     *LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/CCOUP/ZETA(50000),LZETA(8)
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON/SEL/KVAR(500)
      REWIND 14
      REWIND 18
      IPF=1
      LNGT=0
      INDX=1
      REWIND 15
      REWIND 17
      DO 1 I=1,MEMAX
      ELM(I)=0.
      ELMU(I)=0.
  1   ELML(I)=0.
      IEXH=1
 700  G=0.
      D=0.
 200  READ(15,*)IEX,A,B,C,E
      IF(IEX.NE.IEXH)GO TO 400
      G=G+E*A/C/C
      D=D+A*A/C/C
      LNGT=LNGT+1
      GO TO 200
 400  EP(IEXH)=G/D
      TLBDG(IEXH)=G
      VINF(IEXH)=D
      IEXH=IEX
      BACKSPACE 15
      IF(IEX.NE.0)GO TO 700
      REWIND 15
      IEXP=1
 3    G1=0.
      G2=0.
      INH=INDX
      IEXH=IEXP
 5    READ(18,*)LU,INDX,IEXP,AL
      IF(INDX.EQ.0)GO TO 601
      READ(17,*)NI,NF,SGM,AL1
      READ(15,*)IEX,A,B,C,E
      IF(IEXP.EQ.1.OR.IPF.EQ.1)GO TO 800
 801  READ(15,*)IEX,A,B,C,E
      IF(IEXP.NE.IEX)GO TO 801
      IPF=1
 800  CONTINUE
      IF(INDX.NE.INH)GO TO 601
      DY=AL*AL1/B/C
      G1=E*DY+G1
      G2=-2.*DY*A+G2
      WRITE(14,*)INDX,IEXP,NI,NF,DY,A,E,C
      GO TO 5
 601  LOC=(IEXH-1)*LP2+INH
      IPF=0
      ZETA(LOC)=(VINF(IEXH)*G1+TLBDG(IEXH)*G2)/VINF(IEXH)/VINF(IEXH)
      INH=INDX
      REWIND 15
      BACKSPACE 17
      BACKSPACE 18
      IF(INDX.NE.0)GO TO 3
      WRITE(14,*)INDX,IEXP,NI,NF,DY,A,E,B
      REWIND 14
      REWIND 17
 606  READ(14,*)INDX,IEXP,NI,NF,DY,A,E,B
      IF(INDX.EQ.0)GO TO 603
      LOC=(IEXP-1)*LP2+INDX
      SGM=(E-A*EP(IEXP))/B
      U=DY*EP(IEXP)*B+A*ZETA(LOC)/B
      WRITE(17,*)INDX,IEXP,SGM,NI,NF,U
      GO TO 606
 603  WRITE(17,*)INDX,IEXP,SGM,NI,NF,U
      REWIND 17
      LL=0
      CH=0.
  607 READ(17,*)INDX,IEXP,SGM,NI,NF,U
      IF(INDX.EQ.0)GO TO 608
      LL=LL+1
      CH=CH+SGM*SGM
      UX=2.*SGM*U
      ELM(INDX)=ELM(INDX)+UX
      ELMU(INDX)=ELMU(INDX)+ABS(UX)
      GO TO 607
 608  CONTINUE
      WRITE(22,111)
      WRITE(22,211)CH/LL
      DO 212 I=1,NEXPT
 212  WRITE(22,213)I,EP(I)
      WRITE(22,217)
      DO 11 I=1,MEMAX
      IF(KVAR(I).EQ.0)GO TO 11
      RL=ALOG10(ELMU(I)/ABS(ELM(I)))
      IF(RL.GE.RLR)ELML(I)=1.
      WRITE(22,112)I,RL,ELMU(I)/LNGT
 11   CONTINUE
      WRITE(22,115)
      DO 12 I=1,MEMAX
      IF(KVAR(I).EQ.0)GO TO 12
      LC=0
      IF(ELML(I).LT..5)GO TO 12
      REWIND 17
 13   READ(17,*)INDX,IEXP,SGM,NI,NF,AL
      IF(INDX.EQ.0)GO TO 14
      IF(INDX.NE.I)GO TO 13
      LC=LC+1
      IY(LC,1)=IEXP
      IY(LC,2)=NI
      IY(LC,3)=NF
      CORF(LC,1)=SGM
      CORF(LC,2)=AL
      GO TO 13
 14   CONTINUE
      NP=0
      NM=0
      DO 15 J=1,LC
      U=CORF(J,1)*CORF(J,2)*2.
      IF(ABS(U)/ELMU(I).LT..05)GO TO 15
      IF(U.LT.0.)NM=NM+1
      IF(U.GT.0.)NP=NP+1
 15   CONTINUE
      WRITE(22,113)I,NP,NM
      DO 16 L=1,K
      UMP=0.
      UMM=0.
      DO 17 J=1,LC
      U=2.*CORF(J,1)*CORF(J,2)
      IF(U.LT.0.)GO TO 18
      IF(U.LT.UMP)GO TO 19
      UMP=U
      JP=J
 19   GO TO 17
 18   IF(U.GT.UMM)GO TO 17
      UMM=U
      JM=J
 17   CONTINUE
      WRITE(22,114)IY(JP,1),IY(JP,2),IY(JP,3),CORF(JP,1),
     *CORF(JP,2),UMP,IY(JM,1),IY(JM,2),IY(JM,3),
     *CORF(JM,1),CORF(JM,2),UMM
      CORF(JP,1)=0.
      CORF(JM,1)=0.
 16   CONTINUE
 12   CONTINUE
      RETURN
 111  FORMAT(2X////40X,20HTROUBLESHOOTING ROUT,
     *25HINE HAS BEEN ACTIVATED...//5X,
     *31HLOCAL MINIMUM ANALYSIS FOLLOWS://)
 112   FORMAT(6X,1I3,18X,1F4.1,20X,1E7.2)
 115   FORMAT(2X////40X,35HANALYSIS OF SIGNIFICANT DEPENDENCES
     *//)
 113  FORMAT(1X/5X,10(1H*),5X,4HM.E.,1X,1I3,5X,1I3,1X,
     *19HPOSITIVE COMPONENTS,20X,1I3,1X,
     *19HNEGATIVE COMPONENTS///
     *30X,8HPOSITIVE,52X,8HNEGATIVE//5X,3HEXP,2X,
     *10HTRANSITION,2X,5HSIGMA,3X,10HDERIVATIVE,3X,
     *17HD(SIGMA**2)/D(ME),4X,1HI,1X,3HEXP,2X,
     *10HTRANSITION,2X,5HSIGMA,3X,10HDERIVATIVE,
     *3X,17HD(SIGMA**2)/D(ME))
 114  FORMAT(6X,1I2,3X,1I2,2H--,1I2,5X,1F4.1,4X,
     *1E9.2,7X,1E9.2,9X,1HI,2X,1I2,3X,1I2,2H--,
     *1I2,5X,1F4.1,4X,1E9.2,7X,1E9.2)
 211  FORMAT(2X//5X,29HCHISQ FOR FIRST GE(LI)S ONLY ,
     *31HWITH INDEPENDENT NORMALIZATION=,1E12.4//
     *5X,15HNORM.CONSTANTS://)
 213  FORMAT(5X,4HEXP.,1X,1I2,5X,2HC=,1E14.6)
 217  FORMAT(1X//5X,4HM.E.,20X,2HRL,20X,8HSTRENGTH,//)
      END
      SUBROUTINE MIXR(NW,IPSW,CHI,CHILO)
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/MIXD/DMIXE(20,2),DMIX(20),IMIX(20),NDL
      COMMON/LOGY/LNY,INTR,IPS1
      IF(NDL.EQ.0)RETURN
      NW=NW+NDL
      DO 1 I=1,NDL
      IT=IMIX(I)
      INX=KSEQ(IT,1)
      INX1=KSEQ(IT,2)
      IF(ABS(ELM(INX1)).LT.1.E-5)ELM(INX1)=1.E-5
      DL=DMIX(I)*ELM(INX)/ELM(INX1)
      IF(IPSW.EQ.1)DMIX(I)=DL
      CHI=CHI+(DL-DMIXE(I,1))**2/DMIXE(I,2)/DMIXE(I,2)
      IF(LNY.EQ.1)CHILO=CHILO+(DMIXE(I,1)*ALOG(
     *ABS(DL/DMIXE(I,1)))/DMIXE(I,2))**2
   1  CONTINUE
      IF(IPSW.EQ.0)RETURN
      WRITE(22,2)
 2    FORMAT(1X//10X,19HE2/M1 MIXING RATIOS/10X,
     *10HTRANSITION,10X,9HEXP.DELTA,10X,10HCALC.DELTA,10X,5HSIGMA/)
      DO 3 I=1,NDL
      DL=(DMIX(I)-DMIXE(I,1))/DMIXE(I,2)
      IT=IMIX(I)
 3    WRITE(22,4)KSEQ(IT,3),KSEQ(IT,4),DMIXE(I,1),DMIX(I),DL
      RETURN
 4    FORMAT(10X,1I2,3H---,1I2,14X,1F7.2,12X,1F7.2,13X,1F5.2)
      END
      SUBROUTINE COORD(WTH,WPH,WTHH,NAA,IFW,PFI,WPI,WTLB,LZ,TYY,TZZ)
      DIMENSION PFI(101),WPI(11,2)
      COMMON/VLIN/XV(51),YV(51),ZV(20),DSG(20),DSE(20),DS
      COMMON/KIN/EPS(50),EROOT(50),IEXP,IAXS(50),FIEX(50,2)
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
       COMMON/SECK/ISKIN(50)
      DATA RADE/57.2957795/
      IF(IFW.NE.0) GO TO 6
      TYY=WTH-WTHH
      TZZ=WTH+WTHH
 6    XTH=WTH/RADE
      XPH=WPH/RADE
      XTHH=WTHH/RADE
      ZL=TAN(XTHH)
      ZA=COS(XTH)
      ZA1=SIN(XTH)
      ZB=COS(XTHH)
      RMASS=XA1(LZ)/XA
      IF(IZ1(LZ).LT.0) RMASS=1./RMASS
      IF(IFW.EQ.2) GO TO 7
      WS=(TZZ-TYY)/(NAA+1)
      IF(IFW.EQ.1) WS=(TZZ-TYY)/(NAA-1)
 7    DO 5 I=1,NAA
      IF(IFW.NE.2) GO TO 27
      XAA=ABS(WTLB)/RADE
      IF(WTLB.GT.0.) GO TO 20
      IF(IZ1(LZ).LT.0) GO TO 32
      IF(XA1(LZ).LE.XA) GO TO 26
      GO TO 34
 32   IF(XA.LE.XA1(LZ)) GO TO 26
 34   IF(ISKIN(LZ).EQ.0) GO TO 28
 26   TTCM=XAA+TASIN(RMASS*SIN(XAA))
      XAA=(3.14159265-TTCM)/2.
      GO TO 20
 28   continue
      TTCM=XAA-TASIN(RMASS*SIN(XAA))
      XAA=ABS(TTCM)/2.
      GO TO 20
 27   IF(IFW.EQ.0) YV(I)=TYY+I*WS
      XAA=(TYY+WS*(I-1))/RADE
      IF(IFW.EQ.1.AND.(I.EQ.1.OR.I.EQ.NAA)) GO TO 15
      IF(IFW.EQ.0) XAA=YV(I)/RADE
 20   continue
      GI=(ZA-COS(XAA)/ZB)/(ZL*ZA1)
      GA=TACOS(GI)
      WPA=ATAN(ZL*SIN(GA)/(ZA1+ZL*COS(GA)*ZA))
      WPA=ABS(WPA)
      IF(IFW.EQ.2) GO TO 25
      IF(IFW.EQ.1) GO TO 10
      WPI(I,1)=(XPH-WPA)*RADE
      WPI(I,2)=(XPH+WPA)*RADE
      GO TO 5
 10   PFI(I)=2.*WPA*RADE
      GO TO 5
 15   PFI(I)=0.
      GO TO 5
 25   FIEX(LZ,1)=(XPH-WPA)
      FIEX(LZ,2)=(XPH+WPA)
 5    CONTINUE
      IF(WTLB.GE.0..OR.IFW.NE.0) GO TO 30
      DO 50 I=1,NAA
      XAA=YV(I)/RADE
      THETB=ATAN(SIN(2.*XAA)/(RMASS-COS(2.*XAA)))*RADE
      IF(THETB.LT.0.) THETB=180.+THETB
      YV(I)=-1.*THETB
      WPI(I,1)=WPI(I,1)+180.
      WPI(I,2)=WPI(I,2)+180.
 50   CONTINUE
 30   CONTINUE
      RETURN
      END
      SUBROUTINE CHMEM(NW,CHI,CHILO)
      COMMON/ME2D/NAMX,IAMX(100),IAMY(100,2),EAMX(100,2)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      IF(NAMX.EQ.0) RETURN
      NW=NW+NAMX
      DO 10 IA=1,NAMX
      IB=IAMX(IA)
      IF(IAMY(IA,1).NE.IAMY(IA,2)) GO TO 12
      DI=(ELM(IB)-EAMX(IA,1))/EAMX(IA,2)
      CHILO=CHILO+(ALOG(ABS(ELM(IB)/EAMX(IA,1)))*
     *ABS(EAMX(IA,1))/EAMX(IA,2))**2
      CHI=CHI+DI*DI
      GO TO 10
 12   DI=(ELM(IB)-EAMX(IA,1))/EAMX(IA,2)
      CHILO=CHILO+(ALOG(ABS(ELM(IB)/EAMX(IA,1)))*
     *ABS(EAMX(IA,1))/EAMX(IA,2))**2
      CHI=CHI+DI*DI
 10   CONTINUE
      RETURN
      END
      SUBROUTINE PTICC(IDR)
      COMMON/CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON/COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/LEV/TAU(75),KSEQ(500,4)
      WRITE(22,100)
      DO 10 L=1,IDR
      IINX=KSEQ(L,1)
      NI=KSEQ(L,3)
      NF=KSEQ(L,4)
      ENET=EN(NI)-EN(NF)
      CONE2=CONV(ENET,2)
      if(abs(spin(ni)-spin(nf)).gt.2.)cone2=0.
      CONM1=0.
      CONE1=0.
      IF(IINX.LE.MULTI(1))CONE1=CONV(ENET,1)
      IF(ABS(SPIN(NI)-SPIN(NF)).LT.2.) CONM1=CONV(ENET,4)
      WRITE(22,101) NI,NF,SPIN(NI),SPIN(NF),ENET,CONE1,CONE2
     *,CONM1
 10   CONTINUE
 100  FORMAT(1X//20X,31HCALCULATED INTERNAL CONVERSION ,
     >29HCOEFFICIENTS FOR E1,E2 AND M1//5X,2HNI,5X,2HNF,
     >7X,2HII,8X,2HIF,9X,11HENERGY(MEV),6X,7HICC(E1),8X,
     >7HICC(E2),8X,
     >7HICC(M1))
 101  FORMAT(5X,I2,5X,I2,7X,F4.1,6X,F4.1,9X,F6.4,8X,E9.4,6X,
     >E9.4,6X,E9.4)
      RETURN
      END
      FUNCTION RNDM(SE)
      IF(SE.GT.32000.)GO TO 1
 2    SE=SE*SE
      U=ALOG10(SE)
      I=IFIX(U)+1
      T=SE/(10.**I)
      R=SQRT(SQRT(SQRT(T)))
      P=SQRT(SQRT(SQRT(.1)))
      RXDM=(R-P)/(1.-P)
      RXDM=10.*RXDM
      AI=FLOAT(IFIX(RXDM))
      RNDM=RXDM-AI
      RETURN
 1    SE=100.*T+.511
      GO TO 2
      END
      SUBROUTINE KONTUR(IDR,CHIS0,CHIL,IFBF,INPO,JJ,SH,BTEN,REM)
      DIMENSION F(3),BTEN(1200)
      COMMON/VLIN/XV(51),YV(51),ZV(20),DSG(20),DSE(20),DS
      COMMON/DFTB/DEVD(500),DEVU(500)
      COMMON/COMME/ELM(500),ELMU(500),ELML(500),SA(500)
      COMMON/CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/HHH/HLM(500)
      COMMON/ILEWY/NWR
      COMMON/LOGY/LNY,INTR,IPS1
      LNY=0
      H=.05*ABS(HLM(JJ))
      IF(INPO.NE.-1)H=ABS(SH)
 555  INTR=0
      SAJJ=ABS(SA(JJ))
      DO 222 L=1,MEMAX
      ELM(L)=HLM(L)
 222  SA(L)=SA(L)/SAJJ
      YV(1)=0.
      XV(1)=HLM(JJ)
      F(3)=1.
      I=1
 111  CONTINUE
      ITL=0
      V=ELMU(JJ)-ELM(JJ)
      IF(SA(JJ).LT.0.)V=ELM(JJ)-ELML(JJ)
      IF(H.GT.V)ITL=1
      IF(H.GT.V)H=V
      I=I+1
      F(1)=F(3)
      DO 2 J=1,MEMAX
  2   ELM(J)=.5*H*SA(J)+ELM(J)
      CALL LIMITS
      CALL FTBM(3,CHIS1,IDR,1,CHILO,BTEN)
      IF(CHIS1.GT.CHIS0)GO TO 66
      IF(INPO.EQ.-1)WRITE(22,52)JJ,ELM(JJ),CHIS1
      IF(CHIS1.GT.CHIL.OR.INPO.EQ.-1)GO TO 66
      IFBF=1
      IX=1
      CHIL=CHIS1
      WRITE(22,51)CHIL
      GO TO 112
  66  CONTINUE
      WW=.5*(CHIS0-CHIS1)*NWR
      IF(WW.GE.REM)GO TO 888
      F(2)=EXP(WW)
      IF(I.EQ.2.AND.F(2).LT..1.AND.abs(XV(1)-HLM(JJ)).lt.1e-9)GO TO 444
      DO 3 J=1,MEMAX
  3   ELM(J)=ELM(J)+.5*SA(J)*H
      V=ELM(JJ)
      CALL LIMITS
      IF(abs(V-ELM(JJ)).gt.1.e-6)ITL=1
      CALL FTBM(3,CHIS2,IDR,1,CHILO,BTEN)
      IF(CHIS2.GT.CHIS0)GO TO 44
      IF(INPO.EQ.-1)WRITE(22,52)JJ,ELM(JJ),CHIS2
      IF(CHIS2.GT.CHIL.OR.INPO.EQ.-1)GO TO 44
      IFBF=1
      IX=2
      CHIL=CHIS2
      WRITE(22,51)CHIL
      GO TO 112
 44   CONTINUE
      WW=.5*(CHIS0-CHIS2)*NWR
      IF(WW.GT.REM)GO TO 888
      F(3)=EXP(WW)
      IF(ITL.EQ.1)WRITE(22,500)JJ
      IF(I.NE.2)GO TO 200
      IF(ITL.EQ.1)GO TO 200
      IF(F(3).LT..1.AND.abs(XV(1)-HLM(JJ)).lt.1.e-9)GO TO 444
      IF(F(1).GT.F(2).AND.F(2).GT.F(3))GO TO 200
      IF(F(1).LT.F(2).AND.F(2).GT.F(3))GO TO 300
      I=1
      XV(1)=ELM(JJ)
      IF(INPO.EQ.-1)H=2.*H
      GO TO 111
 300  D1=F(2)-F(1)
      D2=F(3)-F(1)
      AC=(D2-4.*D1)*H/(D2-2.*D1)/4.
      DO 301 L=1,MEMAX
 301  ELM(L)=(ELM(L)-H*SA(L))+AC*SA(L)
      CALL LIMITS
      XV(1)=ELM(JJ)
      I=1
      CALL FTBM(3,CHIS1,IDR,1,CHILO,BTEN)
      WW=.5*(CHIS0-CHIS1)*NWR
      IF(WW.GE.REM)GOTO 888
      F(3)=EXP(WW)
      GO TO 111
 200  Y=YV(I-1)
      YV(I)=RK4(Y,H,F)
      XV(I)=ELM(JJ)
      IF(NWR*(CHIS2-CHIS0).LT.2..AND.INPO.EQ.-1)H=2.*H
      IF(ITL.EQ.1)GO TO 33
      IF(F(3).LT.1.E-3)GO TO 33
      GO TO 1
  52  FORMAT(5X,4HELM(,1I3,2H)=,1F10.6,5X,6HCHISQ=,1E12.4)
  51  FORMAT(10X,50HBETTER POINT FOUND...MATRIX ELEMENTS WRITTEN ON 17
     *,3X,6HCHISQ=,1E12.4)
 112  REWIND 17
      DO 55 L=1,MEMAX
  55  WRITE(17,*)ELM(L)
      IF(IX.EQ.1)GO TO 66
      IF(IX.EQ.2)GO TO 44
  1   GO TO 111
 444  H=H/2.
      GO TO 555
 33   C=YV(I)
      M=0
      DO 6 L=1,I
      YV(L)=1.00001-YV(L)/C
      IF(M.EQ.0.AND.YV(L).LT..317)M=L
 6    CONTINUE
      X=(XV(M)-XV(M-1))*(.317-YV(M))/(YV(M-1)-YV(M))
      T=XV(M)-X-HLM(JJ)
      IF(T.GE.0.)DEVU(JJ)=T
      IF(T.LT.0.)DEVD(JJ)=T
 500  FORMAT(5X,11HWARNING-ME(,1I3,1H),5X
     *,32HINTEGRATION STOPPED AT THE LIMIT)
      RETURN
 888  WRITE(22,889)JJ
 889  FORMAT(5X,13H** WARNING **,/,2X,3HME=,1I3,2X,
     *59HTOO FAR FROM THE MINIMUM TO CARRY OUT THE ERROR ESTIMATION!,/)
      RETURN
      END
      FUNCTION RK4(Y,H,F)
      DIMENSION F(3)
      RK4=Y+H*(F(1)+4.*F(2)+F(3))/6.
      RETURN
      END
      SUBROUTINE QFIT(QUI,TAU1,TAU2,ENG,XL1,CF,NL,IND)
      COMMON/DIMX/DIX(4),ODL(200)
      DIMENSION TAU1(10),ENG(10),TAU2(10,7),XL1(7),QUI(8,10),CF(8,2)
      CALL GAMATT(QUI,TAU1,TAU2,XL1,NL)
      IND1=5
      IF(IND.EQ.4)IND1=6
      IF(IND.EQ.5)IND1=7
      DO 1 K=1,8
      CO=QUI(K,IND)
      CN=QUI(K,10)
      CM=QUI(K,IND1)
      CA=(ENG(IND1)-ENG(IND))**2
      CB=(ENG(10)-ENG(IND))**2
      D=CA*(CO-CN)-CB*(CO-CM)
      D1=CA*CM*(CO-CN)-CB*CN*(CO-CM)
      D2=CA*CB*(CN-CM)
      CF(K,1)=D1/D
      CF(K,2)=D2/D
  1   CONTINUE
      RETURN
      END
      SUBROUTINE GAMATT(QUI,TAU1,TAU2,XL1,NL)
      COMMON/DIMX/DIX(4),ODL(200)
      DIMENSION TAU1(10),TAU2(10,7),XL1(7),THING3(10)
     *,Q(9),QUI(8,10)
      DO 3 I=1,10
      I1=1
      THING3(I)=0.
 2    THING1=-TAU2(I,I1)*XL1(I1)+THING3(I)
      I1=I1+1
      THING3(I)=THING1
      IF(I1.LE.NL) GO TO 2
  3   CONTINUE
      DO 4 I=1,10
      TAU=TAU1(I)
      THING=THING3(I)
      CALL GCF(TAU,THING,Q)
      DO 5 K=2,9
 5    QUI(K-1,I)=Q(K)
 4    CONTINUE
      END
      SUBROUTINE GCF(TAU,THING,Q)
      COMMON/DIMX/A,R,XL,D,ODL(200)
      DIMENSION F(101),B(4),Q(9)
      B(1)=ATAN2(A,D+XL)
      B(2)=ATAN2(A,D)
      B(3)=ATAN2(R,D+XL)
      B(4)=ATAN2(R,D)
      DO 100 K=1,9
      Q(K)=0.0
      DO 100 J=1,3
      YL=B(J)
      YU=B(J+1)
      DL=(YU-YL)/100.
      DO 90 M=1,101
      XM=YL+DL*(M-1)
      GO TO (10,20,30),J
 10   EX=TAU*(A-(D+XL)*TAN(XM))/SIN(XM)
      GO TO 40
 20   EX=-TAU*XL/COS(XM)
      GO TO 40
 30   EX=TAU*(D*TAN(XM)-R)/SIN(XM)
 40   F(M)=SIN(XM)*(1-EXP(EX))*EXP(THING/COS(XM))
      IF(J.EQ.1) F(M)=F(M)*EXP(-TAU*(A/SIN(XM)-D/COS(XM)))
      GO TO (90,50,60,70,80,81,82,83,84),K
 50   F(M)=F(M)*COS(XM)
      GO TO 90
 60   F(M)=F(M)*(1.5*COS(XM)**2-0.5)
      GO TO 90
 70   F(M)=F(M)*(2.5*COS(XM)**3-1.5*COS(XM))
      GO TO 90
 80   F(M)=F(M)*(4.375*COS(XM)**4-3.75*COS(XM)**2+.375)
      GO TO 90
 81   F(M)=F(M)*((63.*COS(XM)**5-70.*COS(XM)**3+15.)/8.)
      GO TO 90
 82   F(M)=F(M)*((21.*COS(XM)**2*(11.*COS(XM)**4-15.*COS(XM)**
     *2+5.)-5.)/16.)
      GO TO 90
 83   F(M)=F(M)*(429.*COS(XM)**7-693.*COS(XM)**5+315.*COS(XM)**3
     *-35.*COS(XM))/16.
      GO TO 90
 84   F(M)=F(M)*(6435.*COS(XM)**8-12012.*COS(XM)**6+6930.*COS(XM)
     ***4-1260.*COS(XM)**2+35.)/128.
 90   CONTINUE
      EV=0.0
      OD=0.0
      DO 95 M=2,98,2
      EV=EV+F(M)
 95   OD=OD+F(M+1)
      FINT=DL/3.*(F(1)+4.*(EV+F(100))+2.*OD+F(101))
      Q(K)=Q(K)+FINT
 100  CONTINUE
      DO 101 I=1,8
 101  Q(I+1)=Q(I+1)/Q(1)
      Q(1)=Q(1)/2.
      RETURN
      END
      COMPLEX FUNCTION TCEXP(Z)
      COMPLEX Z
      A=REAL(Z)
      B=AIMAG(Z)
      A=EXP(A)
      C=A*COS(B)
      D=A*SIN(B)
      TCEXP=CMPLX(C,D)
      RETURN
      END
      FUNCTION TCABS(Z)
      COMPLEX Z
      A=REAL(Z)
      B=AIMAG(Z)
      IF(ABS(A).LT.1.E-16)A=0.
      IF(ABS(B).LT.1.E-16)B=0.
      TCABS=SQRT(A*A+B*B)
      RETURN
      END
      FUNCTION TASIN(X)
      TEST=ABS(X)-1.
      IF(ABS(TEST).LT.1.E-9)GO TO 1
      DOL=SQRT(1.-X*X)
      WAR=X/DOL
      TASIN=ATAN(WAR)
      RETURN
   1  TASIN=1.570796327
      IF(X.LT.0.)TASIN=-1.570796327
      RETURN
      END
      FUNCTION TACOS(X)
      TACOS=1.570796327-TASIN(X)
      RETURN
      END
      SUBROUTINE OPENF
      CHARACTER NAME*60,OPT1*20,OPT2*20
	common/switch/jzb,ibps
 1    READ(jzb,*)I,J,K
      IF (I.EQ.0) RETURN
      IF (J.EQ.1) OPT1='OLD'
      IF (J.EQ.2) OPT1='NEW'
      IF (J.EQ.3) OPT1='UNKNOWN'
      IF (K.EQ.1) OPT2='FORMATTED'
      IF (K.EQ.2) OPT2='UNFORMATTED'
      READ(jzb,1000)NAME
 1000 FORMAT(A)
      if(i.ne.25.and.i.ne.26)
     *OPEN(I,IOSTAT=K,FILE=NAME,STATUS=OPT1,FORM=OPT2)
      if(i.eq.25.or.i.eq.26)k=0
c      IF (K.EQ.0) WRITE(6,1030) 'OPENED ',NAME
c 1030 FORMAT (1X,2A)
c      WRITE(6,1010) ' IO-num = ',I,OPT1,OPT2
c 1010 FORMAT (1X,A,I4,2(1x,A))
      IF (K.EQ.0) GO TO 1
c      WRITE (6,1020) 'PROBLEMS OPENING ',NAME,K
c 1020 FORMAT(A,A,I6)
      RETURN
      END
      subroutine effix(ipd,en,effi)
      dimension xx(51),yy(51)
      common/efcal/abc(8,10),akavka(8,200),thick(200,7)
	effi=1.e-6
      en=en+1.e-24
      enl=alog(en)
      do 1 i=1,10
      ll=11-i
      j=ll
      if(enl.ge.abc(8,ll))go to 2
      j=-1
  1   continue
  2   continue
      if(j.eq.-1)effi=1.e-10
      if(j.eq.-1)return
      if(j.eq.1.or.j.eq.10)go to 3
      if(j.eq.9)go to 5
      xx(1)=abc(8,j)
      xx(2)=abc(8,j+1)
      xx(3)=abc(8,j+2)
      go to 6
  3   s=0.
      do 10 l=1,7
      if(abs(thick(ipd,l)).lt.1.e-9)go to 10
      t=exp(abc(l,j))
      d=thick(ipd,l)
      s=s+t*d
 10   continue
      go to 20
  5   xx(1)=abc(8,8)
      xx(2)=abc(8,9)
      xx(3)=abc(8,10)
  6   continue
      s=0.
      do 11 l=1,7
      if(abs(thick(ipd,l)).lt.1.e-9)go to 11
      if(j.eq.9)go to 12
      yy(1)=abc(l,j)
      yy(2)=abc(l,j+1)
      yy(3)=abc(l,j+2)
      go to 13
  12  yy(1)=abc(l,8)
      yy(2)=abc(l,9)
      yy(3)=abc(l,10)
  13  continue
      call lagran(xx,yy,3,0,enl,t,1,1)
      s=s+exp(t)*thick(ipd,l)
  11  continue
  20  effi=exp(-s)
c FITEFF or GREMLIN check
      if(akavka(8,ipd).le.-999.)go to 1501
      if(akavka(5,ipd).gt.0..and.akavka(5,ipd).lt.10.)goto1301
      if(akavka(5,ipd).ge.10.)go to 1500
c GREMLIN
      w=alog(20.*en)
      pw=akavka(1,ipd)+akavka(2,ipd)*w+akavka(3,ipd)*w*w+
     *akavka(4,ipd)*w*w*w
      effi=effi*exp(pw)
      if(abs(akavka(5,ipd)).lt.1.e-9)go to 21
      n=ifix(akavka(6,ipd)+.1)
      pw=w**n
      w=akavka(5,ipd)/pw
      effi=effi*exp(w)
  21  continue
      if(abs(akavka(8,ipd)).lt.1.e-9)return
      w=(akavka(7,ipd)-1000.*en)/akavka(8,ipd)
      pw=exp(w)
      if(abs(pw-1.).lt.1.e-6)write(22,23)
      effi=effi/(1.-pw)
  23  format(5x,'***** CRASH - EFFIX *****')
      return
c FITEFF eff. calib. by P.Olbratowski use
c PJN@2000    
 1301 w=alog(en/akavka(5,ipd))
      pw=akavka(2,ipd)*w
      if(en.lt.akavka(5,ipd))pw=pw+w*w*(akavka(3,ipd)+w*akavka(4,ipd))
      effi=effi*exp(pw)*akavka(1,ipd)
      return
c     JAERI calibration - TC, Nov.2000 
 1500 w=alog(en/.511)
      effi=exp(akavka(1,ipd)+akavka(2,ipd)*w-exp(akavka(3,ipd)+
     *akavka(4,ipd)*w))
      return
C     LEUVEN CALIBRATION
 1501 effi=akavka(1,ipd)
      enek=1000.*en
	w=alog(enek) 
      do 1502 ixi=1,6
	effi=effi+akavka(ixi+1,ipd)*w**ixi	
 1502 continue
      effi=exp(effi)
	return	
      end
      SUBROUTINE ADHOC(OPH,IDR,NFD,NTAP,IYR)
      CHARACTER*4 OPH
      common/clust/iclust(50,200),lastcl(50,20),sumcl(20,500)
     *,irawex(50)
      common/switch/jzb,ibps
      common/cccds/ndst(50)
      COMMON/IDENT/BEQ
      common/efcal/abc(8,10),akavka(8,200),thick(200,7)
      COMMON/DIMX/DIX(4),ODL(200)
      COMMON/TRA/DELTA(500,3),ENDEC(500),ITMA(50,200),ENZ(200)
      COMMON/LIFE/NLIFT
      COMMON/MIXD/DMIXE(20,2),DMIX(20),IMIX(20),NDL
      COMMON/ME2D/NAMX,IAMX(100),IAMY(100,2),EAMX(100,2)
      COMMON/LIFE1/LIFCT(50),TIMEL(2,50)
      COMMON/ERRAN/KFERR
      COMMON/MGN/LP1,LP2,LP3,LP4,LP6,LP7,
     *LP8,LP9,LP10,LP11,LP12,LP13,LP14
      COMMON/SECK/ISKIN(50)
      COMMON/BRNCH/BRAT(50,2),IBRC(2,50),NBRA
      COMMON/YEXPT/YEXP(32,1500),IY(1500,32),CORF(1500,32),DYEX(32,1500)
     *,NYLDE(50,32),UPL(32,50),YNRM(32,50),IDRN,ILE(32)
      COMMON/YTEOR/YGN(500),YGP(500),IFMO
      COMMON/LEV/TAU(75),KSEQ(500,4)
      COMMON/CCC/eg(50),cc(50,5),NANG(200),Q(3,200,8),NICC,
     *AGELI(50,200,2)
      COMMON/GGG/G(7)
      COMMON /CLCOM/LAMDA(8),LEAD(2,500),LDNUM(8,75),LAMMAX,MULTI(8)
      COMMON /COEX/EN(75),SPIN(75),ACCUR,DIPOL,ZPOL,ISO,ACCA
      COMMON/MINNI/IMIN,LNORM(50)
      COMMON/CX/NEXPT,IZ,XA,IZ1(50),XA1(50),EP(50),TLBDG(50),VINF(50)
      COMMON /CEXC/MAGEXC,MEMAX,LMAXE,MEMX6,IVAR(500)
      COMMON/PRT/IPRM(20)
      COMMON/CLCOM0/IFAC(75)
      COMMON/CLCOM9/ERR
      COMMON/COEX2/NMAX,NDIM,NMAX1
      COMMON/THTAR/ITTE(50)
      COMMON/SKP/JSKIP(50)
      COMMON/TRB/ITS
      IOSR=0
      read(jzb,*)IFMO
      read(jzb,*)NICC,NISTR
      read(jzb,*)(EG(JICC),JICC=1,NICC)
      IYR=1
      DO 951 JIC=1,NISTR
      read(jzb,*)ISRT1
      IF(ISRT1.GT.6)ISRT1=ISRT1-3
 951  read(jzb,*)(CC(JICC,ISRT1),JICC=1,NICC)
      read(jzb,*)(NANG(JICC),JICC=1,NEXPT)
      REWIND 9
      READ(9,*)NFD
      DO 954 JICC=1,NFD
      READ(9,*)ODL(JICC)
      READ(9,*)ENZ(JICC)
      DO 954 ISRT1=1,8
  954 READ(9,*)(Q(LICC,JICC,ISRT1),LICC=1,3)
      DO 952 JIC=1,NEXPT
      JUF=NANG(JIC)
      IF(JUF.LT.0)GO TO 955
      read(jzb,*)(ITMA(JIC,JICC),JICC=1,JUF)
      read(jzb,*)(AGELI(JIC,JICC,1),JICC=1,JUF)
      read(jzb,*)(AGELI(JIC,JICC,2),JICC=1,JUF)
      GO TO 952
  955 JUF=IABS(JUF)
      DO 956 JICC=1,JUF
      AGELI(JIC,JICC,1)=AGELI(JIC-1,JICC,1)
      AGELI(JIC,JICC,2)=AGELI(JIC-1,JICC,2)
      ITMA(JIC,JICC)=ITMA(JIC-1,JICC)
  956 CONTINUE
      IF(OPH.NE.'GOSI')NANG(JIC)=IABS(NANG(JIC))
  952 CONTINUE
      CALL SEQ(IDR)
      DO 953 JIC=1,NEXPT
      JUF=NANG(JIC)
      JUF=IABS(JUF)
      DO 953 JICC=1,JUF
      DO 953 LXT=1,2
 953  AGELI(JIC,JICC,LXT)=AGELI(JIC,JICC,LXT)*.0174532925
      TAU(1)=1.E+25
      read(jzb,*)NS1,NS2
      DO 735 LI=1,IDR
      IF(KSEQ(LI,3).EQ.NS1.AND.KSEQ(LI,4).EQ.NS2)GO TO 736
 735  CONTINUE
 736  IDRN=LI
      IF(OPH.NE.'GOSI')RETURN
      DO 737 LI=1,NEXPT
      JUF=NANG(LI)
      IF(JUF.LT.0)GO TO 957
      read(jzb,*)ndst(li)
      ndas=ndst(li)
      read(jzb,*)(UPL(JICC,LI),JICC=1,ndas)
      read(jzb,*)(YNRM(JICC,LI),JICC=1,ndas)
      GO TO 737
  957 JUF=IABS(JUF)
      NANG(LI)=JUF
      NDST(LI)=NDST(LI-1)
      DO 958 JICC=1,JUF
      UPL(JICC,LI)=UPL(JICC,LI-1)
  958 YNRM(JICC,LI)=YNRM(JICC,LI-1)
 737  CONTINUE
      read(jzb,*)NTAP
      IF(NTAP.EQ.0)GO TO 2266
      IPRI=IPRM(2)
      CALL READY(IDR,NTAP,IPRI)
      NDTP=0
      DO 87 IEXP1=1,NEXPT
      JUF=ndst(IEXP1)
      DO 87 IUF=1,JUF
 87   NDTP=NDTP+NYLDE(IEXP1,IUF)
      NVARE=0
      DO 88 IEXP1=1,MEMAX
      IF(IVAR(IEXP1).EQ.1.OR.IVAR(IEXP1).EQ.2)NVARE=NVARE+1
 88   CONTINUE
      WRITE(22,89)NDTP,NVARE
 89   FORMAT(1X//5X,1I4,1X,19HEXPERIMENTAL YIELDS,10X,
     *1I3,1X,28HMATRIX ELEMENTS TO BE VARIED///)
 2266 read(jzb,*)NBRA,WBRA
      IF(ITS.NE.2)GO TO 2623
      REWIND 18
      WRITE(18,*)MEMAX
 2623 CONTINUE
      IF(NBRA.EQ.0)GO TO 7321
      WRITE(22,431)
 431  FORMAT(40X,16HBRANCHING RATIOS,//5X,3HNS1,5X,3HNF1,
     *5X,3HNS2,5X,3HNF2,5X,10HRATIO(1:2),9X,5HERROR)
      DO 730 LB=1,NBRA
      read(jzb,*)NS1,NS2,NS3,NS4,BRAT(LB,1),BRAT(LB,2)
      BRAT(LB,2)=BRAT(LB,2)/(SQRT(WBRA)+1.E-10)
      WRITE(22,432)NS1,NS2,NS3,NS4,BRAT(LB,1),BRAT(LB,2)
 432  FORMAT(5X,1I2,6X,1I2,6X,1I2,6X,1I2,5X,1F10.5,5X,1F10.5)
      DO 731 LI=1,IDR
      IF(KSEQ(LI,3).EQ.NS3.AND.KSEQ(LI,4).EQ.NS4)GO TO 733
      IF(KSEQ(LI,3).EQ.NS1.AND.KSEQ(LI,4).EQ.NS2)GO TO 732
      GO TO 731
 732  IBRC(1,LB)=LI
      GO TO 731
 733  IBRC(2,LB)=LI
 731  CONTINUE
      IF(ITS.NE.2)GO TO 730
      N1=IBRC(1,LB)
      N2=IBRC(2,LB)
      WRITE(18,*)KSEQ(N1,1),KSEQ(N2,1)
      WRITE(18,*)KSEQ(N1,1),KSEQ(N2,2)
      WRITE(18,*)KSEQ(N1,1),KSEQ(N1,2)
      WRITE(18,*)KSEQ(N2,1),KSEQ(N1,2)
      WRITE(18,*)KSEQ(N2,1),KSEQ(N2,2)
      IF(KSEQ(N1,2).NE.0.AND.KSEQ(N2,2).NE.0)WRITE(18,*)KSEQ(N1,2),
     *KSEQ(N2,2)
 730  CONTINUE
      WRITE(22,433)WBRA
 433  FORMAT(5X,38HBRANCHING RATIOS ARE TAKEN WITH WEIGHT,2X,1E14.6)
 7321 read(jzb,*)NLIFT,WLF
      IF(NLIFT.EQ.0)GO TO 7388
      WRITE(22,7325)
      DO 7323 ILFT=1,NLIFT
      read(jzb,*)LIFCT(ILFT),TIMEL(1,ILFT),TIMEL(2,ILFT)
      TIMEL(2,ILFT)=TIMEL(2,ILFT)/(SQRT(WLF)+1.E-10)
 7323 WRITE(22,7324)LIFCT(ILFT),TIMEL(1,ILFT),TIMEL(2,ILFT)
      WRITE(22,7326)WLF
 7325 FORMAT(1X///30X,15HLIFETIMES(PSEC)///5X,
     *5HLEVEL,9X,8HLIFETIME,5X,5HERROR/)
 7324 FORMAT(7X,1I2,6X,1F10.2,3X,1F10.2)
 7326 FORMAT(1X/10X,31HLIFETIMES ARE TAKEN WITH WEIGHT,
     *2X,1E14.6)
 7388 read(jzb,*)NDL,WDL
      IF(NDL.EQ.0)GO TO 7451
      WRITE(22,7404)
      DO 7400 LI=1,NDL
      read(jzb,*)NS1,NS2,DMIXE(LI,1),DMIXE(LI,2)
      DMIXE(LI,2)=DMIXE(LI,2)/(SQRT(WDL)+1.E-10)
      WRITE(22,7405)NS1,NS2,DMIXE(LI,1),DMIXE(LI,2)
      DO 7401 LB=1,IDR
      IF(KSEQ(LB,3).EQ.NS1.AND.KSEQ(LB,4).EQ.NS2)GO TO 7402
      GO TO 7401
 7402 IMIX(LI)=LB
      DMIX(LI)=.8326*(EN(NS1)-EN(NS2))
      IF(ITS.EQ.2)WRITE(18,*)KSEQ(LB,1),KSEQ(LB,2)
 7401 CONTINUE
 7400 CONTINUE
      WRITE(22,7406) WDL
 7406 FORMAT(/10X,41HE2/M1 MIXING RATIOS ARE TAKEN WITH WEIGHT,2X,
     11E14.6)
 7451 CONTINUE
      IF(ITS.EQ.2)WRITE(18,*)IOSR,IOSR
      read(jzb,*)NAMX,WAMX
      IF(NAMX.EQ.0)RETURN
      WRITE(22,7453)
 7453 FORMAT(1X//30X,30HEXPERIMENTAL MATRIX ELEMENT(S)///
     110X,10HTRANSITION,10X,7HMAT.EL.,10X,5HERROR/)
      DO 7460 IAX=1,NAMX
      read(jzb,*)llia,NS1,NS2,EAMX(IAX,1),EAMX(IAX,2)
      IAMY(IAX,1)=NS1
      IAMY(IAX,2)=NS2
      EAMX(IAX,2)=EAMX(IAX,2)/(SQRT(WAMX)+1.E-10)
      WRITE(22,7405) NS1,NS2,EAMX(IAX,1),EAMX(IAX,2)
      iamx(iax)=mem(ns1,ns2,llia)
 7460 CONTINUE
      WRITE(22,7465) WAMX
 7465 FORMAT(/10X,40H MATRIX ELEMENT(S) ARE TAKEN WITH WEIGHT,
     12X,1E14.6)
 7404 FORMAT(1X//20X,32HEXPERIMENTAL E2/M1 MIXING RATIOS///
     *10X,10HTRANSITION,12X,5HDELTA,10X,5HERROR/)
 7405 FORMAT(10X,1I2,3H---,1I2,14X,1F9.4,8X,1F9.4)
      RETURN
      END
      function elmt(xi1,xi2,lam,nb1,nb2,xk1,xk2,xm1,xm2,xm3)
      la=lam
      if(la.gt.6)la=la-6
      xlam=float(la)
      i1=ifix(2.*xi1)
      i2=ifix(2.*xi2)
      llam=2*la
      k1=ifix(2.*xk1)
      k2=ifix(2.*xk2)
      fac=sqrt(2.*xi1+1.)*sqrt(2.*xi2+1.)
C-----In-band matrix element
      if(nb1.ne.nb2)go to 1
C-----K=0
      if(k1.ne.0)go to 2
      elmt=fac*wthrej(i1,llam,i2,0,0,0)*xm1
      return
C-----In band, K.ne.0
  2   pha1=(-1.)**((i1-llam+k1)/2)
      pha2=(-1.)**((k1+i1)/2+1)*pha1
      elmt=fac*(pha1*wthrej(i1,llam,i2,k1,0,-k1)*xm1+
     *pha2*wthrej(i1,llam,i2,-k1,2*k1,-k1)*xm2)
      return
C-----Interband, K-allowed             
C-----One K=0
  1   continue
      if(abs(k1-k2).ge.llam)go to 4
      if(k1.ne.0.and.k2.ne.0)go to 3 
      ipha=(i1-llam+k2)/2
      if(k2.eq.0)ipha=((i2-llam+k1)/2)
      pha1=(-1.)**ipha
      elmt=fac*pha1*wthrej(i1,llam,i2,0,k2,-k2)*xm1
      if(k2.eq.0)elmt=fac*pha1*wthrej(i2,llam,i1,0,k1,-k1)*xm1       
      if(k1.ne.0.or.k2.ne.0)elmt=elmt*sqrt(2.)
      return
C-----Both K's non-zero
  3   continue
      pha1=(-1.)**((i1-llam+k2)/2)
      pha2=(-1.)**((i1+k1)/2)*pha1
      elmt=fac*(pha1*wthrej(i1,llam,i2,k1,k2-k1,-k2)*xm1+
     *pha2*wthrej(i1,llam,i2,-k1,k1+k2,-k2)*xm2)  
      return
C-----Forbidden and K1-K2=lambda, Mikhailov formula
  4   continue
      addt=0.
      if(k1.ne.1)go to 11
      addt=(-1.)**((i1+1)/2)*(i1+1)/2.*xm3
 11   continue     
      xn=abs(xk1-xk2)-xlam
      n=ifix(xn+.1)
      if(n.eq.0)go to 5
      if(n.eq.1)go to 6
      s1=xi1-xk1      
      s2=xi1+xk1+1.
      do 7 l=1,n
      s1=s1*(xi1-xk1-float(l))
  7   s2=s2*(xi1+xk2+1.+float(l))
      fct=sqrt(s1*s2)
      go to 9
  6   fct=sqrt((xi1-xk1)*(xi1+xk1+1.))
      go to 9
  5   fct=1.
  9   pha1=(-1.)**ifix((xi1-xlam+xk2)+.1)
      elmt=fac*pha1*fct*wthrej(i1,llam,i2,k2-llam,llam,-k2)*
     *(xm1+xm2*(xi2*(xi2+1.)-xi1*(xi1+1.))+addt)
      return
      end
      SUBROUTINE OPENF1
      CHARACTER NAME*60,OPT1*20,OPT2*20
 1    READ *,I,J,K
      IF (I.EQ.0) RETURN
      IF (J.EQ.1) OPT1='OLD'
      IF (J.EQ.2) OPT1='NEW'
      IF (J.EQ.3) OPT1='UNKNOWN'
      IF (K.EQ.1) OPT2='FORMATTED'
      IF (K.EQ.2) OPT2='UNFORMATTED'
      READ 1000,NAME
 1000 FORMAT(A)
      OPEN(I,IOSTAT=K,FILE=NAME,STATUS=OPT1,FORM=OPT2)
c      IF (K.EQ.0) WRITE(6,1030) 'OPENED ',NAME
c 1030 FORMAT (1X,2A)
c      WRITE(6,1010) ' IO-num = ',I,OPT1,OPT2
c 1010 FORMAT (1X,A,I4,2(1x,A))
      IF (K.EQ.0) GO TO 1
c      WRITE (6,1020) 'PROBLEMS OPENING ',NAME,K
c 1020 FORMAT(A,A,I6)
      RETURN
      END
C-----------------------------------------------------------------------
