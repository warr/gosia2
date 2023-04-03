
C----------------------------------------------------------------------
C SUBROUTINE MINI
C
C Called by: GOSIA
C Calls:     FTBM, LIMITS
C
C Purpose: perform the minimization
C
C Uses global variables:
C      CHIS11 - chi squared
C      CORF   - internal correction factors
C      DEVD   -
C      DEVU   -
C      DLOCK  - limit of derivative below which matrix element fixed if LOCKS=1
C      ELM    - matrix elements
C      ELMH   - matrix element
C      GRAD   - partial derivative of chi squared wrt. matrix element
C      HLMLM  - old value of matrix element or chi squared
C      ICS    - read internal correction factors from file rather than recalculating
C      IFBFL  - calculate derivatives with forward-backward method
C      INTR   - flag to swap chisqr and log(chisqr)
C      IPRM   - printing flags (see suboption PRT of OP,CONT)
C      IPS1   - terminate after calculating and writing correction factors
C      ITAK2  -
C      IUNIT3 - unit for TAPE3 (gosia2)
C      IVAR   - fixed, correlated or free flag
C      JENTR  - flag set to 0 normally, 1 in OP,ERRO
C      KFERR  - error flag for minimization
C      KVAR   -
C      LFL1   -
C      LNY    - use logs to calculate chi squared
C      LOCKF  - fix those with most significat derivative
C      LOCKS  - lock flag. if LOCKS=1, fix at first stage of minimisation
C      LP4    - 1500
C      LP6    - 32
C      MEMAX  - number of matrix elements
C      NLOCK  - number of matrix elements to lock
C      NWR    - number of datapoints used in fit
C      SA     - ratio of elements for correlated elements
C
C Formal parameters:
C      Chisq  - chi squared of minimization
C      Chiok  - desired chi squared
C      Nptl   - number of iterations allowed
C      Conv   - parameter for convergence test
C      Imode  - mode of minimization
C      Idr    -
C      Xtest  -
C      Ips    -
C      Is     - generate input for sigma program flag
C      Jjh    -
C      Bten   -
C
C FTBM does the main calculation and LIMITS makes sure the matrix elements
C don't go outside the limits specified by the user.

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
      INTEGER*4 IBPS, IUNIT3, JZB ! For gosia2
      DIMENSION ipm(10) , Bten(1200) , gradp(1500)
      COMMON /DUMM  / GRAD(1500) , HLMLM(1500) , ELMH(1500)
      COMMON /ILEWY / NWR
      COMMON /CH1T  / CHIS11
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 ,
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /UWAGA / ITAK2
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) ,
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) ,
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /DFTB  / DEVD(1500) , DEVU(1500)
      COMMON /PRT   / IPRM(20)
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      COMMON /SEL   / KVAR(1500)
      COMMON /FIT   / LOCKF , NLOCK , IFBFL , LOCKS , DLOCK
      COMMON /ERRAN / KFERR
      COMMON /LOGY  / LNY , INTR , IPS1
      COMMON /ERCAL / JENTR , ICS
      COMMON /SWITCH/ JZB , IBPS , IUNIT3 ! For gosia2
      DATA chirf/0./,dm/0./,sumg2/0./

C     Initialise gradp to zero for each matrix element
      DO i = 1 , MEMAX
         gradp(i) = 0.
      ENDDO

C     Initialise some parameters
      icount = 0
      lnm = 0
      LNY = 0
      INTR = 0
      metf = 0
      LFL1 = 0
      ncall = 0
      ITAK2 = 0
      chil = 1.D38

C     Handle the different modes
C     Imode = IJKL, where
C        I=1 => fast approximation to calculate chi squared and its partial derivatives
C        I=2 => full Coulomb excitation formalism, but derivatives with fast approximation
C
C        J=0 => steepest descent minimization
C        J=1 => gradient minimization
C
C        K=0 => absolute changes of matrix elements
C        K=1 => relative changes
C
C        L=0 => yields, branching ratios used to calculate chi squared
C        L=1 => logs used to claculate chi squared
      IF ( Imode.LT.2000 ) THEN ! Fast approximation for chi squared and derivatives
         icl1 = 0
         icl2 = 3
         IF ( Imode.GE.1100 ) metf = 1
         IF ( (Imode-1000-100*metf).GE.10 ) lnm = 1
         IF ( (Imode-1000-100*metf-10*lnm).EQ.1 ) LNY = 1 ! Use logs
         IF ( JENTR.EQ.1 ) GOTO 200 ! If we are in OP,ERRO, jump
         IF ( ICS.NE.0 ) THEN ! Read correction factors from file, rather than recalculating
            REWIND 11
            DO jnm = 1 , LP4 ! LP4 is 1500
               READ (11) (CORF(jnm,kh2),kh2=1,LP6) ! LP6 is 32
            ENDDO
            ICS = 0
            GOTO 200
         ENDIF
      ELSE ! Full Coulomb excitation formalism for chi squared, fast approx for derivatives
         icl1 = 1
         IF ( Imode.GE.2100 ) metf = 1
         IF ( (Imode-2000-100*metf).GE.10 ) lnm = 1
         IF ( (Imode-2000-100*metf-10*lnm).EQ.1 ) LNY = 1 ! Use logs
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

C     Call FTBM to perform a single calculation
 100  CALL FTBM(0,chiss,Idr,0,chl,Bten)

C     Write correction factors
      REWIND 11
      DO jnm = 1 , LP4
         WRITE (11) (CORF(jnm,kh2),kh2=1,LP6)
      ENDDO

      IF ( IPS1.EQ.0 ) RETURN ! If IPS1 = 0, terminate after writing correction factors

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
         IF ( uxa.LT.Chiok ) GOTO 600
 350     ino = 1
         IF ( metf.EQ.1 ) ipas = ipas + 1
         IF ( IFBFL.EQ.1 ) ino = 2 ! IFBFL = 1 means use forward-backward method
         DO jjj = 1 , ino
            DO jnm = 1 , MEMAX
               GRAD(jnm) = 0.
               IF ( IVAR(jnm).EQ.1 .OR. IVAR(jnm).EQ.2 ) THEN
                  DO jcoup = 1 , MEMAX
                     ELMH(jcoup) = ELM(jcoup)
                  ENDDO
                  DO jcoup = 1 , MEMAX
                     IF ( jnm.NE.jcoup ) THEN
                        IF ( IVAR(jcoup).LT.1000 ) GOTO 355
                        jcp = IVAR(jcoup) - 1000
                        IF ( jcp.NE.jnm ) GOTO 355
                        IF ( IVAR(jnm).EQ.0 ) GOTO 355
                     ENDIF
                     flt = 1.01
                     IF ( jjj.EQ.2 ) flt = .99
                     ELM(jcoup) = ELMH(jcoup)*flt
 355              ENDDO
                  CALL FTBM(3,chis13,Idr,ncall,chx,Bten)
                  IF ( jjj.EQ.1 ) HLMLM(jnm) = chis13
                  IF ( IFBFL.NE.1 .OR. jjj.NE.1 ) THEN
                     IF ( jjj.EQ.2 ) chis12 = chis13
                     GRAD(jnm) = 100.*(HLMLM(jnm)-chis12)/ELMH(jnm)
                     IF ( IFBFL.EQ.1 ) GRAD(jnm) = GRAD(jnm)/2. ! Forward-backward
                     IF ( lnm.EQ.1 ) GRAD(jnm) = GRAD(jnm)
     &                    *ABS(ELMH(jnm))
                  ENDIF
                  DO jcoup = 1 , MEMAX
                     ELM(jcoup) = ELMH(jcoup)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         IF ( KFERR.EQ.1 ) THEN
            GRAD(Jjh) = 0.
            IF ( Is.EQ.1 .AND. icount.EQ.1 ) WRITE (IUNIT3,*) ! For sigma program (unit changed for gosia2)
     &           (NWR*GRAD(jnm),jnm=1,MEMAX)
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
     &              GRAD(jnm)*GRAD(jnm)
            ENDDO
            IF ( sumg2.LT.1.E-10 ) GOTO 700
            sumg2 = SQRT(sumg2)
            DO jnm = 1 , MEMAX
               GRAD(jnm) = GRAD(jnm)/sumg2
            ENDDO
            IF ( metf.NE.0 ) THEN
               dm = 0.
               DO jnm = 1 , MEMAX
                  IF ( IVAR(jnm).EQ.2 .OR. IVAR(jnm).EQ.1 ) dm = dm +
     &                 ELM(jnm)*ELM(jnm)*GRAD(jnm)*GRAD(jnm)
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
               GOTO 350
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
     &           WRITE (22,99010) Chisq
            WRITE (*,99010) Chisq
            IF ( MOD(icount,IPRM(6)).EQ.0 ) THEN
               WRITE (22,99002)
99002          FORMAT (20X,'GRADIENT'//)
               nlinn = MEMAX/10 + 1
               DO jlin = 1 , nlinn
                  jsa = (jlin-1)*10 + 1
                  DO jin = 1 , 10
                     ipm(jin) = jsa + jin - 1
                  ENDDO
                  jst = MIN(jsa+9,MEMAX)
                  jpr = MIN(10,MEMAX-jsa+1)
                  WRITE (22,99003) (ipm(jin),jin=1,jpr)
99003             FORMAT (5X,10(5X,1I3,4X))
                  WRITE (22,99004) (GRAD(jin),jin=jsa,jst)
99004             FORMAT (5X,10(1X,1E10.4,1X)/)
               ENDDO
            ENDIF
         ENDIF
         IF ( chil.LT.Chiok ) GOTO 600 ! We've achieved desired chi square
         DO l = 1 , MEMAX
            HLMLM(l) = ELM(l)
         ENDDO
         DO l = 1 , MEMAX
            IF ( ABS(GRAD(l)).LE.DLOCK .AND. LOCKS.EQ.1 .AND.
     &           icount.EQ.1 .AND. IVAR(l).LE.999 .AND. IVAR(l).NE.0 )
     &           THEN
               IF ( KFERR.NE.1 ) KVAR(l) = 0
               IF ( KFERR.NE.1 ) WRITE (22,99005) l , GRAD(l)
99005          FORMAT (1X,'MATRIX ELEMENT',1X,1I3,1X,'LOCKED',3X,
     &                 'DERIVATIVE=',1E14.6)
               IVAR(l) = 0
            ENDIF
         ENDDO
         istec = 0
      ENDIF

 400  DO j = 1 , MEMAX
         ELMH(j) = ELM(j)
      ENDDO

C     Find steepest gradient
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

 500  DO j = 1 , MEMAX
         ELM(j) = ELMH(j) - ht*GRAD(j)
      ENDDO

      DO j = 1 , MEMAX
         IF ( IVAR(j).GE.1000 ) THEN  ! For correlated elements
            indx1 = IVAR(j) - 1000    ! Index of element to which it is correlated
            ELM(j) = ELM(indx1)*SA(j) ! SA is the ratio we require
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
         GOTO 500
      ELSE
         CALL LIMITS
         CALL FTBM(icl2,Chisq,Idr,ncall,chilo,Bten)
         IF ( Chisq.GE.chil ) THEN
            ht = ht/2.
            IF ( ABS(ht).GE.Conv ) GOTO 500
         ELSE
            chil = Chisq
            sumht = sumht + ht
            IF ( ABS(ht/sumht).GE..01 ) GOTO 400
         ENDIF
         crit = 0.
         DO jjj = 1 , MEMAX
            crit = crit + (ELM(jjj)-HLMLM(jjj))**2
         ENDDO
         crit = SQRT(crit)
         IF ( crit.LT.Conv ) GOTO 700
         IF ( Chisq.GE.Chiok ) THEN
            rfk = chirf/Chisq
            IF ( rfk.LE.Xtest .OR. icount.GE.Nptl ) GOTO 300
            GOTO 100
         ENDIF
      ENDIF

C     Required chi square achieved
 600  chil = Chisq
      IF ( Ips.EQ.0 ) WRITE (22,99006) icount
99006 FORMAT (5X,'AT STEP',1X,1I5,1X,'CHISQ CRITERION FULFILLED')
      IF ( Ips.EQ.0 ) WRITE (22,99010) chil
      RETURN

 700  IF ( LOCKF.EQ.0 ) THEN ! Terminate if convergence satisfied
         IF ( Chisq.GE.chil ) THEN
            DO jjj = 1 , MEMAX
               ELM(jjj) = ELMH(jjj)
            ENDDO
         ENDIF
         IF ( KFERR.EQ.1 ) RETURN
         IF ( Ips.EQ.0 ) WRITE (22,99007) icount , crit
99007    FORMAT (5X,'AT STEP',1X,1I5,'CONVERGENCE ACHIEVED(',1E14.6,')')
         IF ( Ips.EQ.0 ) WRITE (22,99010) MIN(chil,Chisq)
      ELSE ! Fix most significant chi squared derivatives
         DO kkk = 1 , NLOCK ! NLOCK is number of derivatives to fix
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
      ENDIF
      INTR = 0
      RETURN

99010 FORMAT (5X,'*** CHISQ=',1E14.6,1X,'***')
99011 FORMAT (1X/5X,'MATRIX ELEMENT',1X,1I3,1X,'LOCKED!')
      END
