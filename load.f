 
C----------------------------------------------------------------------
C SUBROUTINE LOAD
C
C Called by: FTBM, GOSIA
C Calls:     LSLOOP
C
C Purpose: calculates various parameters, xi, psi which are stored in
C variables XI (common CXI) and PSI (common PCOM).
C
C Uses global variables:
C      CAT    - substates of levels (n_level, J, m)
C      DIPOL  - E1 polarization parameter
C      EMMA   - Controls number of magnetic substates in full coulex calc.
C      EN     - energy of level
C      EP     - bombarding energy
C      ERR    - error flag
C      IPATH  - index of substate in level with same m as substate Irld
C      ISHA   - is half-integer spin
C      ISMAX  - number of substates used
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of not-investated nucleus
C      LAMDA  - list of multipolarities to calculate
C      LAMMAX - number of multipolarities to calculate
C      LDNUM  - number of matrix elements with each multipolarity populating level
C      LEAD   - pair of levels involved in each matrix element
C      LMAX   - ground-state spin + 1
C      LMAXE  - maximum multipolarity needed for calculation
C      LP7    - maximum number of zeta coefficients (45100)
C      LP10   - maximum number of magnetic substates (1200)
C      LZETA  - index into ZETA array for zeta for a given multipolarity
C      MAGA   - number of magnetic substates in approximate calculation
C      MAGEXC - flag: 0 means no magnetic excitations, 1 means with mag. exc.
C      MEMAX  - number of matrix elements
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C      PSI    - psi coefficients
C      QAPR   - approximate Coulomb amplitudes
C      SPIN   - spin of level
C      VINF   - speed of projectile at infinity
C      XA     - A of investigated nucleus
C      XA1    - A of not-investated nucleus
C      XI     - xi coupling constants
C      ZPOL   - dipole term
C
C Formal parameters:
C      Iexp   - Number of experiment
C      Ient   - Flag : 1, 2, 3 (read only)
C      Icg    - Flag : 1, 2 (read only)
C      Polm   - (read only)
C      Joj    - index of substate (write only)
 
      SUBROUTINE LOAD(Iexp,Ient,Icg,Polm,Joj)
      IMPLICIT NONE
      REAL*8 a1 , a2 , aaz2 , aaz3 , aazz , ah , cpsi , dep , eta , 
     &       etan , Polm , pp1 , pp2
      REAL*8 ppp , rlam , ssqrt , szet , wrt , wrtm , z1 , z2 , zet , 
     &       zsqa
      INTEGER*4 i , i1 , i2 , i3 , Icg , Ient , Iexp , ir , is , ispi , 
     &          ispo
      INTEGER*4 jj , jjj , Joj , la , lam , lam1 , ld , m , m1
      INTEGER*4 m2 , mstop , n , n2 , n3 , nn , nz
      INCLUDE 'clcom.inc'
      INCLUDE 'mgn.inc'
      INCLUDE 'coex.inc'
      INCLUDE 'pspin.inc'
      INCLUDE 'cexc.inc'
      INCLUDE 'pcom.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'clm.inc'
      INCLUDE 'clcom8.inc'
      INCLUDE 'clcom9.inc'
      INCLUDE 'coex2.inc'
      INCLUDE 'cx.inc'
      INCLUDE 'cexc0.inc'
      INCLUDE 'cxi.inc'
      INCLUDE 'caux0.inc'
      INCLUDE 'aprcat.inc'
      INCLUDE 'pth.inc'
      DIMENSION etan(100) , cpsi(8)
      
      LMAX = INT(SPIN(1)+1.1) ! ground-state spin + 1
      IF ( Ient.EQ.1 ) THEN
         ISHA = 0
         ispi = INT(SPIN(1)+.51)
         ispo = INT(SPIN(1)+.49)
         IF ( ispi.NE.ispo ) ISHA = 1 ! Half-integer spin
         z1 = DBLE(ABS(IZ1(Iexp)))
         z2 = DBLE(IZ)
         a1 = XA1(Iexp)
         a2 = XA
         ZPOL = DIPOL*EP(Iexp)*a2/(z2*z2*(1.+a1/a2))
         IF ( IZ1(Iexp).LE.0 ) THEN
            ah = a1
            a1 = a2
            a2 = ah
         ENDIF

C        Calculate xi and store it in XI in common CXI
C        The value 6.349770 is 197.33/1.44*sqrt(2/931.49).
C        i.e. hbar c / e^2 * sqrt(2 / amu).
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

C        Calculate C_\lambda \over (s Z_1 Z_2)^\lambda 
         aazz = 1./(1.+a1/a2)/z1/z2
         cpsi(1) = 5.169286*aazz

C        Electric excitations up to LMAXE
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

C        Magnetic excitations
         IF ( MAGEXC.NE.0 ) THEN ! If magnetic excitations are required
            aazz = VINF(Iexp)/95.0981942
            cpsi(7) = aazz*cpsi(1)
            IF ( LAMMAX.NE.8 ) cpsi(8) = aazz*cpsi(2)
         ENDIF

C        Calculate psi and store in PSI in common PCOM
         zsqa = z1*SQRT(a1)
         i3 = 1
         ppp = 1. + a1/a2
         DO i1 = 1 , LAMMAX ! For each calculated multipolarity
            lam = LAMDA(i1)
            lam1 = lam
            IF ( lam.GT.6 ) lam1 = lam - 6
            DO n2 = 1 , NMAX ! For each level
               nn = LDNUM(lam,n2) ! Number of matrix elements connected to this level by this multipolarity
               IF ( nn.NE.0 ) THEN
                  n3 = LEAD(1,i3)
                  pp1 = EP(Iexp) - ppp*EN(n3)
                  DO m1 = 1 , nn ! For each matrix element connected to level with this multipolarity
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

C     Initialise NSTART and NSTOP arrays
      DO n = 1 , NMAX ! For each level
         NSTART(n) = 0
         NSTOP(n) = 0
      ENDDO

      is = 1
      NSTART(1) = 1
      DO n = 1 , NMAX ! For each level
         wrt = Polm - EMMA(Iexp)
         wrtm = Polm + EMMA(Iexp)
         IF ( Icg.EQ.2 ) wrt = Polm - DBLE(MAGA(Iexp))
         IF ( Icg.EQ.2 ) wrtm = Polm + DBLE(MAGA(Iexp))
         IF ( wrtm.LT.-SPIN(n) ) THEN
            NSTART(n) = 0 ! Level has no interesting substates
         ELSE
            IF ( ABS(wrt).GT.SPIN(n) ) wrt = -SPIN(n)
            IF ( wrtm.GT.SPIN(n) ) wrtm = SPIN(n)
            mstop = INT(wrtm-wrt+1.01)
            DO i = 1 , mstop ! For each substate
               CAT(is,1) = n               ! Number of level
               CAT(is,2) = SPIN(n)         ! Spin of level
               CAT(is,3) = wrt + DBLE(i-1) ! m quantum number of substate
               IF ( n.EQ.1 .AND. ABS(CAT(is,3)-Polm).LT.1.E-6 ) Joj = is
               is = is + 1
            ENDDO ! Loop on substates i
         ENDIF
         NSTART(n+1) = is ! First substate of level n+1
         NSTOP(n) = is - 1 ! Last substate of level n
      ENDDO ! Loop on levels n

      ISMAX = is - 1 ! ISMAX is the number of substates used
      IF ( ISMAX.LE.LP10 ) THEN ! LP10 is max. number of substates (1200)
         IF ( Ient.EQ.3 ) RETURN
         nz = 0
         DO jj = 1 , 7
            DO jjj = 1 , MEMAX ! For each matrix element
               QAPR(jjj,1,jj) = 0.
               QAPR(jjj,2,jj) = 0.
            ENDDO
         ENDDO

C        Initialise pointers to ZETA array
         DO i = 1 , 8
            LZETA(i) = 0
         ENDDO

         DO i1 = 1 , LAMMAX ! For each multipolarity
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
                  n = CAT(ir,1) ! number of level for substate ir
                  IF ( Icg.NE.1 ) THEN
                     IF ( MAGA(Iexp).EQ.0 .AND. ir.NE.IPATH(n) ) GOTO 10
                     IF ( ABS(ir-IPATH(n)).GT.1 ) GOTO 10
                  ENDIF
                  ld = LDNUM(la,n) ! Number of matrix elements connected to this level by this multipolarity
                  IF ( ld.EQ.0 ) THEN
                     ir = ir + NSTOP(n) - NSTART(n)
                  ELSE
                     CALL LSLOOP(ir,n,nz,ld,lam,la,ssqrt,Icg,Iexp)
                  ENDIF
                  GOTO 10
               ENDIF
            ENDIF
         ENDDO ! Loop over multipolarity

C        Make sure we haven't reached the limit of zeta coefficients          
         IF ( nz.GT.LP7 ) THEN ! LP7 is 45100
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

      ERR = .TRUE. ! Set error flag
      END
