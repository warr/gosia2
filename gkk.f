 
C----------------------------------------------------------------------
C SUBROUTINE GKK
C
C Called by: GKVAC
C Calls:     ATS, WSIXJ, XSTATIC
C
C Purpose: calculate time-dependent deorientation coefficients
C
C Uses global variables:
C      AKS    - <\alpha_k> values
C      AVJI   - average J  (this is G(1) in GOSIA)
C      DQ     - width of gaussian distribution
C      FILE   - K          (this is G(6) in GOSIA)
C      GAMMA  - Gamma      (this is G(2) in GOSIA)
C      GFAC   - g          (this is G(5) in GOSIA)
C      GKI    - G_k for a single level
C      IBYP   - flag to indicate whether we calculate <\alpha_k>
C      POWER  - x          (this is G(7) in GOSIA)
C      QCEN   - center of gaussian distribution
C      SUM    - sum over 6-j symbol squared
C      TIMEC  - Tau_C      (this is G(4) in GOSIA)
C      XLAMB  - Lambda*    (this is G(3) in GOSIA)
C      XNOR   - normalisation factor
C
C Formal parameters:
C      Iz     - Z of nucleus
C      Beta   - v/c
C      Spin   - spin of state
C      Time   - lifetime of state
C      Il     - index into AKS array
C
C We start by calling XSTATIC to calculate the static part. This calculates
C QCEN (the centre of the gaussian charge state distribution), DQ (the
C gaussian width of this distribution) and XNOR (the normalization parameter
C such that the sum over probabilities is one).
C
C We calculate:
C <a_k> = \sum_l p(J_1) \sum_F {(2 F + 1)^2 \over 2 J_1 + 1} *
C                              {\sixj{F F k I I J_1}}^2.
C
C We include a correction to take into account the effect of nuclear lifetimes
C which are comparable to the mean time between random reorientations \tau_c.
C
C Note that certain values have defaults:
C AVJI = 3, GAMMA = 0.02, XLAMB = 0.0345, TIMEC = 3.5, GFAC = Z/A,
C FIEL = 6E-6 and POWER = 0.7, which are set in GOSIA, where they are treated
C as an array called G in the order of the values in the GGG common block.
C However, the user may change them using the VAC suboption of the CONT option
C of OP,COUL or OP,GOSI.
C
C Note that WSIXJ requires all its parameters to be doubled, so it can handle
C half-integers properly.
C
C The function ATS is used to determine the ground-state spin for a given
C element.

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
         CALL XSTATIC(Iz,inq,ifq,Beta) ! inq and ifq are range of integral
         l = 0
         DO i = 1 , 6
            AKS(i,Il) = 0.
         ENDDO
 50      IF ( imean.EQ.1 ) inq = 1
         IF ( imean.EQ.1 ) ifq = 1

         DO j = inq , ifq
            l = l + 1
            nz = Iz - j
            xji = ATS(nz) ! Ground-state spin of atom
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
 100     ENDDO ! Loop on j
         imean = imean + 1
         IF ( imean.EQ.1 ) GOTO 50
      ENDIF

      hmean = FIEL*Iz*(Beta**POWER) ! Mean magnetic field in fluctuating state
      wsp = 4789.*GFAC*hmean/AVJI ! 4789 is the nuclear magneton
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
         up = up*XLAMB*Time + 1.       ! numerator
         down = Time*(xlam+XLAMB) + 1. ! denominator = r
         GKI(k) = up/down
         alp = 9.*xlam*xlam + 8.*xlam*TIMEC*(w2-xlam*xlam)
         alp = SQRT(alp) - 3.*xlam
         alp = alp/4./xlam/TIMEC                      ! alp is p
         upc = xlam*Time*(down-2.*alp*alp*Time*TIMEC) ! numerator
         dwc = (down+alp*Time)*(down+2.*alp*Time)     ! denominator
         ccf = 1. + upc/dwc                           ! ccf is correction factor
         GKI(k) = GKI(k)*ccf
      ENDDO
      END
