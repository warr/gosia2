 
C----------------------------------------------------------------------
C SUBROUTINE SEQ
C
C Called by: ADHOC
C Calls:     CONV, DECAY, F, GF, LEADF, MEM
C
C Purpose: in order to calculate the yields, we need to start with the highest
C level and calculate its yield, so as to work out the feeding for the lower
C levels and take this into account, gradually working our way down to the
C ground state.
C
C Uses global variables:
C      DELTA  - \delta_\lambda: index 1 = electric^2, 2 = magnetic^2, 3 = cross term
C      EN     - energy of level
C      ENDEC  - energy difference for each matrix element
C      FP     - F coefficient * DELTA^2
C      GKP    - Gk * DELTA^2
C      KLEC   - number of decays for each level
C      KSEQ   - indices for each decay (level1, level2, matrix element, multipolarity + 10)
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LP2    - maximum number of matrix elements (500)
C      LP3    - maximum number of levels (75)
C      MULTI  - number of matrix elements having a given multipolarity
C      NMAX   - number of levels
C      NMAX1  - number of levels with decays
C      SPIN   - spin of level
C      TAU    - normally lifetime in picoseconds (here it is used for energies, however)
C
C Formal parameters:
C      Idr    - returns number of items in KSEQ array.
C
C We store the order in the KSEQ array of common block LEV.
C
C Note that in the code, a multipolarity 1 = E1, 2 = E2 ... 6 = E6, 7 = M1,
C 8 = M2.
 
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
      COMMON /TRA   / DELTA(1500,3) , ENDEC(1500) , ITMA(50,200) , 
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

      DO l = 1 , LP3 ! LP3 = 75 (number of levels)
         KLEC(l) = 0 ! Initialise KLEC to zero
      ENDDO

      DO k = 1 , LP2 ! LP2 = 500 (number of matrix elements)
         DO j = 1 , 3
            DO l = 1 , 4
               FP(l,k,j) = 0.
               IF ( j.NE.3 ) GKP(l,k,j) = 0.
            ENDDO
            DELTA(k,j) = 0.
         ENDDO
      ENDDO

C     Store the energies in TAU array
      DO n = 1 , NMAX
         TAU(n) = EN(n)
      ENDDO

      DO n = 1 , NMAX ! Loop on levels
C        Find level with highest energy
         emax = 0.
         DO j = 1 , NMAX ! Loop on levels
            IF ( TAU(j).GE.emax ) THEN
               emax = TAU(j)
               jsave = j
            ENDIF
         ENDDO
         DO is = 1 , NMAX ! Loop on levels
            DO la = 1 , 8 ! Loop on multipolarities
               IF ( la.LE.3 .OR. la.EQ.7 .OR. la.EQ.8 ) THEN ! E3, M1, M2
                  ld = LDNUM(la,is) ! Number of levels connected to this one with this multipolarity
                  IF ( ld.NE.0 ) THEN
                     DO ir = 1 , ld ! For each level ir connected to level is with multipolarity la
                        m = LEADF(is,ir,la)
                        IF ( m.EQ.jsave .OR. is.EQ.jsave ) THEN
                           IF ( is.NE.jsave .OR. EN(m).LT.EN(is) ) THEN
                              IF ( m.NE.jsave .OR. EN(is).LT.EN(m) )
     &                             THEN
                                 indx = MEM(is,m,la) ! Matrix element from level is to level m with multipolarity la
                                 idecay = idecay + 1
                                 KSEQ(idecay,1) = m       ! Level
                                 KSEQ(idecay,2) = is      ! Level
                                 KSEQ(idecay,3) = indx    ! Matrix element
                                 KSEQ(idecay,4) = la + 10 ! Multipolarity + 10
                                 IF ( EN(m).LE.EN(is) ) THEN ! If the levels are degenerate, swap order
                                    KSEQ(idecay,1) = is
                                    KSEQ(idecay,2) = m
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO ! Loop on levels ir
                  ENDIF
               ENDIF
            ENDDO ! Loop on multipolarity la
         ENDDO ! Loop on levels is
         TAU(jsave) = -1.
      ENDDO ! Loop on levels n

C     Now for each decay, calculate transition amplitudes for each
C     multipolarity
      DO l = 1 , idecay ! For each decay
         istr1 = 0
         IF ( KSEQ(l,4).LT.10 ) GOTO 200 ! KSEQ(l,4) is 10 + multipolarity
         istr2 = 0
         n = KSEQ(l,1) ! Initial level
         m = KSEQ(l,2) ! Final level
         inx = KSEQ(l,3) ! Matrix element
         la = KSEQ(l,4) - 10 ! Multipolarity
         ega = EN(n) - EN(m)    ! ega = E_\gamma
         twoi = 1./SQRT(2.*SPIN(n)+1.)
         spini = SPIN(n) + .001
         spinf = SPIN(m) + .001
         egs = SQRT(ega)*twoi   ! egs = \sqrt{E_\gamma \over 2 I_1 + 1}
         js = l + 1
         la1 = 0
         inx1 = 0
         DO j = js , idecay ! For each decay
            IF ( KSEQ(j,4).GE.10 ) THEN ! KSEQ(1,4) is 10 + multipolarity
               n1 = KSEQ(j,1) ! Initial level
               m1 = KSEQ(j,2) ! Final level
               IF ( n1.EQ.n .AND. m1.EQ.m ) THEN ! Decays involving the same pair of levels
                  inx1 = KSEQ(j,3) ! Matrix element
                  la1 = KSEQ(j,4) - 10 ! Multipolarity
                  KSEQ(j,4) = KSEQ(j,4) - 10 ! Subtract ten to indicate we have handled this one
               ENDIF
            ENDIF
         ENDDO ! Loop on decays j
         KSEQ(l,4) = KSEQ(l,4) - 10 ! Subtract ten to indicate we have handled this one
         Idr = Idr + 1
         mule = 0
         mulm = 0
         nob = 1
 50      IF ( la.LE.3 ) THEN
            IF ( la.EQ.1 ) THEN
               DELTA(Idr,1) = 398.77393*ega*egs ! E1
               mule = 1
               istr1 = 1
            ELSEIF ( la.EQ.2 ) THEN
               DELTA(Idr,1) = 3.5002636*egs*ega*ega ! E2
               mule = 2
               istr1 = 2
            ELSEIF ( la.EQ.3 ) THEN
               DELTA(Idr,1) = 0.023891302*ega*ega*ega*egs ! E3
               mule = 3
               istr1 = 3
            ELSE
               GOTO 100
            ENDIF
            GOTO 150
         ENDIF
 100     la = la - 6
         IF ( la.EQ.2 ) THEN
            DELTA(Idr,2) = 0.036806836*ega*ega*egs ! M2
            mulm = 2
            istr2 = 5
         ELSE
            DELTA(Idr,2) = 4.1932861*ega*egs ! M1
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
         ENDEC(Idr) = EN(n) - EN(m) ! Energy difference between levels
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
         ENDDO ! Loop on mk
         DELTA(Idr,1) = DELTA(Idr,1)*(1.+CONV(ega,istr1))
         DELTA(Idr,2) = DELTA(Idr,2)*(CONV(ega,istr2)+1.)
         KLEC(n) = KLEC(n) + 1 ! Increment KLEC for initial level
 200  ENDDO ! Loop on decays l

      NMAX1 = 0
      DO n = 1 , NMAX ! For each level count those which have decays
         IF ( KLEC(n).NE.0 ) NMAX1 = NMAX1 + 1
      ENDDO
      END
