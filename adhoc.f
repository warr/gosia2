
C----------------------------------------------------------------------
C SUBROUTINE ADHOC
C
C Called by: GOSIA
C Calls: READY, SEQ
C
C Purpose: to handle the OP,YIEL option.
C
C Uses global variables:
C      AGELI  - angles of the Ge detectors
C      BRAT   - branching ratio and its error
C      CC     - conversion coefficients for different energies and multipolarities
C      DMIX   - 0.8326 * gamma energy
C      DMIXE  - mixing ratio and its error
C      EAMX   - known matrix elements and their error
C      EG     - energies for conversion coefficients
C      EN     - energy of level
C      ENZ    - depends on absorber
C      IAMX   - index of matrix element for known matrix element
C      IAMY   - level indices of pair of levels for which matrix element is known
C      IBRC   - index branching ratios
C      IDRN   - index of normalising transition for yields
C      IFMO   - include correction to angular distance for finite recoil distance.
C      IMIX   - decay associated with known mixing ratio
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      ITMA   - identify detectors according to OP,GDET
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      IVAR   - indicates a limit or correlation is set
C      JZB    - unit to read from
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LIFCT  - index for lifetimes
C      MEMAX  - number of matrix elements
C      NAMX   - number of known matrix elements
C      NANG   - number of gamma-ray detectors for each experiment
C      NBRA   - number of branching ratios
C      NDL    - number of mixing ratios
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NICC   - number of conversion coefficients
C      NLIFT  - number of lifetimes
C      NYLDE  - number of yields
C      ODL    - results of OP,GDET calculation
C      Q      - solid angle attenuation coefficients
C      TAU    - lifetime in picoseconds
C      TIMEL  - lifetimes and their errors
C      UPL    - upper limits for all gamma detectors
C      YNRM   - relative normalization factors for gamma detectors
C
C Formal parameters:
C      Oph    - this indicates the option (GOSI, STAR etc.)
C      Idr    - number of decays
C      Nfd    - number of physical detectors
C      Ntap   - unit of yield file
C      Iyr    - flag set here
C
C Here we parse the input of the OP,YIEL command and store the values.

      SUBROUTINE ADHOC(Oph,Idr,Nfd,Ntap,Iyr)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AGELI , BRAT , CC , CORF , DELTA , DIPOL ,
     &       DIX , DMIX , DMIXE , DYEX , EAMX , EG , EN , ENDEC , ENZ ,
     &       EP , ODL , Q
      REAL*8 SPIN , TAU , TIMEL , TLBDG , UPL , VINF , wamx , wbra ,
     &       wdl , wlf , XA , XA1 , YEXP , YGN , YGP , YNRM , ZPOL
      INTEGER*4 IAMX , IAMY , iax , IBRC , Idr , IDRN , iexp1 , IFMO ,
     &          ILE , ilft , IMIX , iosr , ipri , IPRM , ISO , isrt1 ,
     &          ITMA , ITS , iuf , IVAR
      INTEGER*4 IY , Iyr , IZ , IZ1 , jic , jicc , juf , KSEQ , lb ,
     &          li , licc , LIFCT , llia , LMAXE , lxt , MAGEXC , MEM ,
     &          MEMAX , MEMX6 , n1
      INTEGER*4 n2 , NAMX , NANG , NBRA , ndas , NDL , NDST , ndtp ,
     &          NEXPT , Nfd , NICC , nistr , NLIFT , ns1 , ns2 , ns3 ,
     &          ns4 , Ntap , nvare , NYLDE
      INTEGER*4 IBPS , IUNIT3 , JZB ! For gosia2
      INTEGER*4 ISPL ! Added for spline
      CHARACTER*4 Oph
      COMMON /CCCDS / NDST(50)
      COMMON /DIMX  / DIX(4) , ODL(200)
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) ,
     &                ENZ(200)
      COMMON /LIFE  / NLIFT
      COMMON /MIXD  / DMIXE(20,2) , DMIX(20) , IMIX(20) , NDL
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      COMMON /LIFE1 / LIFCT(50) , TIMEL(2,50)
      COMMON /BRNCH / BRAT(50,2) , IBRC(2,50) , NBRA
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) ,
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) ,
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /YTEOR / YGN(1500) , YGP(1500) , IFMO
      COMMON /LEV   / TAU(75) , KSEQ(1500,4)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) ,
     &                NICC , NANG(200) , ISPL
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA ,
     &                ISO
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) ,
     &                TLBDG(50) , VINF(50)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      COMMON /PRT   / IPRM(20)
      COMMON /TRB   / ITS
      COMMON /SWITCH/ JZB , IBPS , IUNIT3 ! For gosia2

C     Read OP,YIEL parameters
      iosr = 0
      READ (JZB,*) IFMO ! IFLAG
      READ (JZB,*) NICC , nistr ! N1, N2
      READ (JZB,*) (EG(jicc),jicc=1,NICC) ! E1,E2...
      Iyr = 1
      DO jic = 1 , nistr
        READ (JZB,*) isrt1 ! I1
         IF ( isrt1.GT.6 ) isrt1 = isrt1 - 3
         READ (JZB,*) (CC(jicc,isrt1),jicc=1,NICC) ! CC(I1,1)...CC(I1,N1)
      ENDDO
      READ (JZB,*) (NANG(jicc),jicc=1,NEXPT) ! NANG(I)...NANG(NEXPT)

C     Read file for gamma-ray energy dependence of Ge solid-angle attenuation
C     coefficients Q
      REWIND 9
      READ (9,*) Nfd
      DO jicc = 1 , Nfd
         READ (9,*) ODL(jicc) ! DIX(4) - distance from target to front of detector
         READ (9,*) ENZ(jicc) ! Depends on absorber
         DO isrt1 = 1 , 8
            READ (9,*) (Q(licc,jicc,isrt1),licc=1,3)
         ENDDO
      ENDDO

C     Read detector identities, theta and phi
      DO jic = 1 , NEXPT ! For each experiment
         juf = NANG(jic)
         IF ( juf.LT.0 ) THEN ! If NANG < 0 use previous values
            juf = ABS(juf) ! Number of detector angles
            DO jicc = 1 , juf ! For each detector angle
               AGELI(jic,jicc,1) = AGELI(jic-1,jicc,1) ! theta same as previous detector
               AGELI(jic,jicc,2) = AGELI(jic-1,jicc,2) ! phi same as previous detector
               ITMA(jic,jicc) = ITMA(jic-1,jicc)
            ENDDO
            IF ( Oph.NE.'GOSI' ) NANG(jic) = ABS(NANG(jic))
         ELSE
            READ (JZB,*) (ITMA(jic,jicc),jicc=1,juf) ! IP(1)...IP(NANG(I))
            READ (JZB,*) (AGELI(jic,jicc,1),jicc=1,juf) ! Theta Ge det
            READ (JZB,*) (AGELI(jic,jicc,2),jicc=1,juf) ! Phi Ge det
         ENDIF
      ENDDO ! Loop jic on experiments

C     Call SEQ to calculate "chronological" order of levels, so we can
C     account for feeding
      CALL SEQ(Idr)

C     Convert angles into radians
      DO jic = 1 , NEXPT ! For each experiment
         juf = NANG(jic)
         juf = ABS(juf) ! Number of detector angles
         DO jicc = 1 , juf ! For each detector angle
            DO lxt = 1 , 2 ! 1 is theta, 2 is phi
               AGELI(jic,jicc,lxt) = AGELI(jic,jicc,lxt)*.0174532925 ! 0.017452925 = pi / 180
            ENDDO
         ENDDO
      ENDDO ! Loop jic on experiments

C     Set normalising transition
      TAU(1) = 1.E+25 ! Initialise ground-state lifetime to 1E25 picoseconds
      READ (JZB,*) ns1 , ns2 ! NS1, NS2
      DO li = 1 , Idr ! Search through decays for right pair of levels
         IF ( KSEQ(li,3).EQ.ns1 .AND. KSEQ(li,4).EQ.ns2 ) GOTO 100
      ENDDO
 100  IDRN = li ! Index of normalising transition
      IF ( Oph.NE.'GOSI' ) RETURN

C     Read upper limits and relative normalisation factors
      DO li = 1 , NEXPT ! Loop on experiments
         juf = NANG(li)
         IF ( juf.LT.0 ) THEN ! If NANG < 0, use same as previous
            juf = ABS(juf)
            NANG(li) = juf ! Number of detector angles
            NDST(li) = NDST(li-1) ! Number of datasets same as previous
            DO jicc = 1 , juf ! For each detector angle
               UPL(jicc,li) = UPL(jicc,li-1) ! Upper limits same as previous
               YNRM(jicc,li) = YNRM(jicc,li-1) ! Relative normalisation same as previous
            ENDDO
         ELSE
            READ (JZB,*) NDST(li) ! NDST
            ndas = NDST(li)
            READ (JZB,*) (UPL(jicc,li),jicc=1,ndas) ! UPL1...N
            READ (JZB,*) (YNRM(jicc,li),jicc=1,ndas) ! YNRM1...N
         ENDIF
      ENDDO ! Loop li on experiments

C     Read file for experimental yields
      READ (JZB,*) Ntap ! NTAP
      IF ( Ntap.NE.0 ) THEN
         ipri = IPRM(2)
         CALL READY(Idr,Ntap,ipri) ! Read yields from unit Ntap
         ndtp = 0
         DO iexp1 = 1 , NEXPT ! Loop on experiments
            juf = NDST(iexp1) ! Number of datasets
            DO iuf = 1 , juf ! Loop on datasets
               ndtp = ndtp + NYLDE(iexp1,iuf)
            ENDDO
         ENDDO ! Loop iexp1 on experiments

C        Count free variables
         nvare = 0
         DO iexp1 = 1 , MEMAX ! For each matrix element
            IF ( IVAR(iexp1).EQ.1 .OR. IVAR(iexp1).EQ.2 )
     &           nvare = nvare + 1
         ENDDO
         WRITE (22,99001) ndtp , nvare
99001    FORMAT (1X//5X,1I4,1X,'EXPERIMENTAL YIELDS',10X,1I3,1X,
     &           'MATRIX ELEMENTS TO BE VARIED'///)
      ENDIF ! IF ( Ntap.NE.0 )

C     Read branching ratios
      READ (JZB,*) NBRA , wbra ! NBRA, WBRA
      IF ( ITS.EQ.2 ) THEN
         REWIND 18
         WRITE (18,*) MEMAX
      ENDIF
      IF ( NBRA.NE.0 ) THEN
         WRITE (22,99002)
99002    FORMAT (40X,'BRANCHING RATIOS',//5X,'NS1',5X,'NF1',5X,'NS2',5X,
     &           'NF2',5X,'RATIO(1:2)',9X,'ERROR')
         DO lb = 1 , NBRA ! I1,I2,I3,I4,B,DB repeated NBRA times
            READ (JZB,*) ns1 , ns2 , ns3 , ns4 , BRAT(lb,1) , BRAT(lb,2)
            BRAT(lb,2) = BRAT(lb,2)/(SQRT(wbra)+1.E-10) ! Relative error
            WRITE (22,99003) ns1 , ns2 , ns3 , ns4 , BRAT(lb,1) ,
     &                       BRAT(lb,2)
99003       FORMAT (5X,1I2,6X,1I2,6X,1I2,6X,1I2,5X,1F10.5,5X,1F10.5)
            DO li = 1 , Idr ! Search decays for these pairs of levels
               IF ( KSEQ(li,3).EQ.ns3 .AND. KSEQ(li,4).EQ.ns4 ) THEN
                  IBRC(2,lb) = li ! Decay index for first pair
               ELSEIF ( KSEQ(li,3).EQ.ns1 .AND. KSEQ(li,4).EQ.ns2 ) THEN
                  IBRC(1,lb) = li ! Decay index for second pair
               ENDIF
            ENDDO
            IF ( ITS.EQ.2 ) THEN
               n1 = IBRC(1,lb) ! Decay of first pair
               n2 = IBRC(2,lb) ! Decay of second pair
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

C     Read lifetimes
      READ (JZB,*) NLIFT , wlf ! NL, WL
      IF ( NLIFT.NE.0 ) THEN
         WRITE (22,99005)
99005    FORMAT (1X///30X,'LIFETIMES(PSEC)'///5X,'LEVEL',9X,'LIFETIME',
     &           5X,'ERROR'/)
         DO ilft = 1 , NLIFT ! INDEX, T, DT repeated NL times
            READ (JZB,*) LIFCT(ilft) , TIMEL(1,ilft) , TIMEL(2,ilft)
            TIMEL(2,ilft) = TIMEL(2,ilft)/(SQRT(wlf)+1.E-10) ! Relative error
            WRITE (22,99006) LIFCT(ilft) , TIMEL(1,ilft) , TIMEL(2,ilft)
99006       FORMAT (7X,1I2,6X,1G10.2,3X,1G10.2)
         ENDDO
         WRITE (22,99007) wlf
99007    FORMAT (1X/10X,'LIFETIMES ARE TAKEN WITH WEIGHT',2X,1E14.6)
      ENDIF

C     Read known mixing ratios
      READ (JZB,*) NDL , wdl ! NDL, WDL
      IF ( NDL.NE.0 ) THEN
         WRITE (22,99008)
99008    FORMAT (1X//20X,'EXPERIMENTAL E2/M1 MIXING RATIOS'///10X,
     &           'TRANSITION',12X,'DELTA',10X,'ERROR'/)
         DO li = 1 , NDL ! IS, IF, DELTA, ERROR repeated NDL times
            READ (JZB,*) ns1 , ns2 , DMIXE(li,1) , DMIXE(li,2)
            DMIXE(li,2) = DMIXE(li,2)/(SQRT(wdl)+1.E-10)
            WRITE (22,99012) ns1 , ns2 , DMIXE(li,1) , DMIXE(li,2)
            DO lb = 1 , Idr ! Search through decays for right pair of levels
               IF ( KSEQ(lb,3).EQ.ns1 .AND. KSEQ(lb,4).EQ.ns2 ) THEN
                  IMIX(li) = lb ! Decay index
                  DMIX(li) = .8326*(EN(ns1)-EN(ns2)) ! 0.8326 * energy of gamma
                  IF ( ITS.EQ.2 ) WRITE (18,*) KSEQ(lb,1) , KSEQ(lb,2)
               ENDIF
            ENDDO
         ENDDO
         WRITE (22,99009) wdl
99009    FORMAT (/10X,'E2/M1 MIXING RATIOS ARE TAKEN WITH WEIGHT',2X,
     &           1E14.6)
      ENDIF
      IF ( ITS.EQ.2 ) WRITE (18,*) iosr , iosr

C     Read known matrix elements
      READ (JZB,*) NAMX , wamx ! NAMX, WAMX
      IF ( NAMX.EQ.0 ) RETURN
      WRITE (22,99010)
99010 FORMAT (1X//30X,'EXPERIMENTAL MATRIX ELEMENT(S)'///10X,
     &        'TRANSITION',10X,'MAT.EL.',10X,'ERROR'/)

      DO iax = 1 , NAMX ! LAMBDA, INDEX1, INDEX2, ME, DME repeated NAMX times
         READ (JZB,*) llia , ns1 , ns2 , EAMX(iax,1) , EAMX(iax,2)
         IAMY(iax,1) = ns1 ! Level index
         IAMY(iax,2) = ns2 ! Level index
         EAMX(iax,2) = EAMX(iax,2)/(SQRT(wamx)+1.E-10) ! Relative error of ME
         WRITE (22,99012) ns1 , ns2 , EAMX(iax,1) , EAMX(iax,2)
         IAMX(iax) = MEM(ns1,ns2,llia) ! Index to matrix element
      ENDDO
      WRITE (22,99011) wamx
99011 FORMAT (/10X,' MATRIX ELEMENT(S) ARE TAKEN WITH WEIGHT',2X,1E14.6)

99012 FORMAT (10X,1I2,'---',1I2,14X,1F9.4,8X,1F9.4)
      END
