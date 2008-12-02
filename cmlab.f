 
C----------------------------------------------------------------------
C SUBROUTINE CMLAB
C
C Called by: GOSIA
C Calls:     TASIN
C
C Purpose: calculate for center of mass frame
C
C Uses global variables:
C      BETAR  - recoil beta
C      DSIGS  - dsigma for each experiment
C      EN     - level energies
C      EP     - bombarding energy
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 - 1)
C      ERR    - error flag
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      ISKIN  - kinematic flag (0,1)
C      IZ     - Z of investigated nucleus
C      IZ1    - Z of not-investated nucleus
C      NCM    - calculate kinematics assuming this state for final state (default = 2)
C      NEXPT  - number of experiments
C      NMAX   - number of level energies
C      TETACM - theta of particle detector in center of mass frame
C      TLBDG  - theta of particle detector in lab frame (in degrees)
C      VINF   - speed of projectile at infinity
C      XA     - A of investigated nucleus
C      XA1    - A of not-investated nucleus
C      TREP   - theta of recoiling nucleus (in radians)
C
C Formal parameters:
C      Ii     - experiment number (or zero for all experiments)
C      Dsig   - dsigma
C      Tetrn  - theta of recoiling nucleus in lab frame (in radians)

      SUBROUTINE CMLAB(Ii,Dsig,Tetrn)
      IMPLICIT NONE
      REAL*8 a1 , a2 , ared , d2a , dista , dists , Dsig , emax , epmin
      REAL*8 r3 , TASIN , tau , taup , tcmdg , tcmrad , Tetrn , tlbrad , 
     &       tmxdg , z1 , z2 , zcmdg , zcmrad , zlbrad
      INTEGER*4 iflaa , Ii , lexp , lexp0 , lexp1 , n
      INCLUDE 'clcom9.inc'
      INCLUDE 'seck.inc'
      INCLUDE 'prt.inc'
      INCLUDE 'tcm.inc'
      INCLUDE 'brec.inc'
      INCLUDE 'caux0.inc'
      INCLUDE 'coex2.inc'
      INCLUDE 'cx.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'coex.inc'
      DATA r3/0./

      lexp0 = 1
      lexp1 = NEXPT
      IF ( Ii.NE.0 ) lexp0 = Ii
      IF ( Ii.NE.0 ) lexp1 = Ii
      DO lexp = lexp0 , lexp1 ! For each experiment
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
C
C        dists is Cline's estimate of the maximum safe bombarding energy
         dists = 1.44*(a1+a2)*z1*z2/((a1**.33333+a2**.33333)*1.25+5.)/a2
C        dista is 0.05 * distance of closest approach for head-on collisions
         dista = 0.0719949*(1.0+a1/a2)*z1*z2/EP(lexp)
C        d2a is the distance of closest approach for head-on collisions in fm
C        q^2/4/pi/epsilon_0 * (1+a1/a2) * Z1 * Z2 / Ep. For Ep in MeV and d2a
C        in fm, q^2/4/pi/epsilon_0 = 1.44
         d2a = 20.0*dista ! = 1.44 * (1.0+a1/a2)*z1*z2/EP(lexp)
C        VINF is the initial velocity of the incoming projectile (at infinity)
C        VINF = sqrt(2 * EP / 931.494028 * A1) : 931.494028 = 1 AMU
         VINF(lexp) = 0.0463365*SQRT(EP(lexp)/a1)

C        If IPRM(1) we want extra printout
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

C        Final kinetic energy \v{E} = E_P - \Delta E (1 + m_P / m_T)
C        Here we set ared = (1 + m_P / m_T)
C        The maximum excitation energy corresponds to \v{E} = 0, so
C        \DeltaE = E_P \over {1 + m_P / m_T) = E_P/ared
C        We check that there are no states defined which are higher than this.

         tlbrad = TLBDG(lexp)/57.2957795 ! Theta of detector to radians
         ared = 1.0 + a1/a2 ! reduced mass
         emax = EP(lexp)/ared ! Maximum excitation energy
         DO n = 1 , NMAX ! For each level
            IF ( EN(n).GT.emax ) GOTO 50 ! Give error if energy of state too high
         ENDDO

C        Gosia calculates assuming the kinematics for all states are approximately
C        those corresponding to the state NCM (by default NCM = 2 : the first excited
C        state). So for this energy we calculate the \v{E} and store it in epmin.
C        We also calculate tau defined as (a1/a2)*sqrt(E_P / \v{E}) for this value
C        of \v{E}.
C        A value of tau less than 1 corresponds to normal kinematics, so the full
C        range of theta in the centre of mass system corresponds to the full range
C        in the lab system. However, for tau greater than 1 (i.e. inverse kinematics)
C        there are two possible values for the lab angle for a given centre of mass
C        angle and there is a maximum lab angle, which can be attained: tmxdg given
C        by SIN(tmxdg) = 1 / tau.
         epmin = EP(lexp) - EN(NCM)*ared
         taup = SQRT(EP(lexp)/epmin)
         tau = taup*a1/a2
         IF ( tau.LE.1.0 ) GOTO 100 ! No limit on scattering angle
         tmxdg = TASIN(1.0/tau)*57.2957795 ! Maximum lab angle in degrees
         IF ( tmxdg.GE.TLBDG(lexp) ) GOTO 100 ! Within limit of scattering angle

         WRITE (22,99007) tmxdg , lexp
99007    FORMAT (1X,'ERROR- MAXIMUM SCATTERING ANGLE IS ',F7.2,
     &           ' DEGREES',' FOR EXPERIMENT ',1I2)
         GOTO 200 ! Error

 50      WRITE (22,99008) emax , lexp
99008    FORMAT (1X,'ERROR- MAXIMUM EXCITATION ENERGY IS ',F8.4,' MEV',
     &           ' FOR EXPERIMENT ',1I2)
         GOTO 200 ! Error

C        Calculate centre of mass angle
 100     tcmrad = tlbrad + TASIN(tau*SIN(tlbrad)) ! In radians
         tcmdg = tcmrad*57.2957795 ! and in degrees

C        In inverse kinematics, for a given lab angle, there are two solutions
C        for the centre of mass angle.
         IF ( tau.GT.1.0 ) THEN ! Inverse kinematics
            IF ( IPRM(1).EQ.1 ) THEN
               IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99009)
     &              tcmdg , lexp
99009          FORMAT (5X,'SECOND POSSIBLE CM SCATTERING ANGLE IS',F7.2,
     &                 ' DEGREES FOR EXPERIMENT ',1I2)
            ENDIF
            IF ( ISKIN(lexp).NE.1 ) THEN ! If ISKIN is set, take the second solution
               tcmdg = 180. + 2.*TLBDG(lexp) - tcmdg
               tcmrad = tcmdg/57.2957795
            ENDIF
         ENDIF

C        EPS is "epsilon" the eccentricity parameter.
         EPS(lexp) = 1./SIN(tcmrad/2.)
         TETACM(lexp) = tcmrad
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99010) tcmdg , 
     &           EPS(lexp)
99010       FORMAT (5X,'CM SCATTERING ANGLE',1X,1F10.3,1X,'DEG',5X,
     &              'EPSILON',1X,1F10.4)
         ENDIF

C        If Z1 is negative, we are interested in target excitations, but if it
C        is positive, we want the projectile excitation, so calculate the lab
C        recoil energy of appropriate particle and store it in BETAR (we will
C        convert this to beta of the recoil later)
         IF ( IZ1(lexp).GT.0 ) BETAR(lexp) = a1*a2/(a1+a2)
     &        **2*(1.+taup*taup-2.*taup*COS(tcmrad))*epmin
         IF ( IZ1(lexp).LT.0 ) BETAR(lexp) = (a2/(a1+a2))
     &        **2*(1.+tau*tau+2.*tau*COS(tcmrad))*epmin

C        More additional printout
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99011)
     &           BETAR(lexp)
99011       FORMAT (5X,'RECOIL ENERGY(MEV)',2X,1F10.4)
         ENDIF

C        This is the beta of the recoiling particle of interest (target or projectile
C        depending on sign of Z1, which is used as a flag)
         BETAR(lexp) = .0463365*SQRT(BETAR(lexp)/XA) ! 0.0463365 = sqrt(2/931.494028)
         IF ( IPRM(1).EQ.1 ) THEN
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99012)
     &           BETAR(lexp)
99012       FORMAT (5X,'RECOIL BETA',2X,1E14.6)
            IF ( Ii.EQ.0 .AND. IPRM(10).EQ.1 ) WRITE (22,99013) EP(lexp)
     &           /(dists*.5*(1.+EPS(lexp)))
99013       FORMAT (5X,'BOMBARDING ENERGY=',1F10.3,1X,
     &              'OF SAFE BOMBARDING ENERGY AT THIS ANGLE')
         ENDIF

C        iflaa = 0 when projectile detected, = 1 when target detected
C        r3 is the Jacobian dOmega/domega
         IF ( iflaa.NE.1 ) THEN ! Projectile detected
            IF ( ABS(tcmdg-180.).LT.1.E-5 ) THEN
               r3 = (1.-tau)**2
            ELSE
               r3 = SIN(tlbrad)/SIN(tcmrad)
               r3 = r3*r3*ABS(COS(tcmrad-tlbrad))
               r3 = 1./r3
            ENDIF
         ENDIF

C        Calculate the values for the target. In the centre of mass system, the
C        target and projectile angles differ by 180 degrees
         zcmdg = 180. - tcmdg ! Target angle in degrees in cm system
         zcmrad = zcmdg/57.2957795 ! and in radians
         zlbrad = ATAN(SIN(zcmrad)/(COS(zcmrad)+taup)) ! target theta in lab (radians)

C        iflaa = 0 when projectile detected, = 1 when target detected
C        r3 is the Jacobian dOmega/domega
         IF ( iflaa.NE.0 ) THEN ! Target detected, but theta is for projectile!
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

C        Now calculate dsigma
         Dsig = 250.*r3*SQRT(EP(lexp)/(EP(lexp)-ared*EN(NCM)))
     &          *dista*dista*(EPS(lexp))**4
         EROOT(lexp) = SQRT(EPS(lexp)*EPS(lexp)-1.)
         DSIGS(lexp) = Dsig
         Tetrn = zlbrad
         IF ( IZ1(lexp).LT.0. ) Tetrn = tlbrad
         TREP(lexp) = Tetrn
      ENDDO ! Loop over experiments lexp

      IPRM(10) = 0 ! Turn off printing so we don't write things twice
      RETURN

C     An error has occured, so set error flag and return
 200  ERR = .TRUE. ! Set error flag
      END
