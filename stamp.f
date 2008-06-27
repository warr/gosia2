 
C----------------------------------------------------------------------
C FUNCTION STAMP
C
C Called by: LAIAMP
C Calls:     TRINT
C
C Purpose: Estimate amplitude
C
C Formal parameters:
C      Epsi   - epsilon for this experiment
C      Errt   - sqrt(epsilon^2 - 1) for this experiment
C      Xiv    - value of xi
C      Dw     - step in omega (0.03)
C      W0     - value of omega
C      Lmda   - lambda (1...6 for E1...6 and 1,2 for M1,2)
C      Mua    - mu
C
C Return value:
C      Estimated amplitude
 
      COMPLEX*16 FUNCTION STAMP(Epsi,Errt,Xiv,Dw,W0,Lmda,Mua)
      IMPLICIT NONE
      REAL*8 a , axi , b , bic , bic2 , bis , bis2 , ca , cb , cia , 
     &       cib , cic , cis , Dw , dwi , Epsi , Errt , ex , exa , fct
      INTEGER*4 la , Lmda , mi , Mua
      REAL*8 sa , sb , sia , sib , W0 , Xiv
      DATA fct/0./

      mi = Mua - 1
      axi = ABS(Xiv) ! Absolute value of xi
      la = Lmda
      IF ( Lmda.EQ.7 ) la = 3

      IF ( axi.LT.1.E-5 ) THEN
         a = -2.*W0
         IF ( la.EQ.3 ) a = -W0
         exa = EXP(a)
         dwi = 3*Dw
         cic = exa*(EXP(dwi)-1.)
         STAMP = DCMPLX(cic,0.D0)
         IF ( la.EQ.2 ) THEN
            IF ( mi.EQ.0 ) fct = 3.*(3.-Epsi*Epsi)/Epsi/Epsi/Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.837117307*Errt/Epsi/Epsi/Epsi/Epsi ! 1.837117307 = sqrt(27/8)
            IF ( mi.EQ.2 ) fct = -3.674234613*Errt*Errt/Epsi/Epsi/Epsi/ ! 3.674234613 = sqrt(27/2)
     &                           Epsi
         ELSEIF ( la.EQ.3 ) THEN
            fct = -1.414213562*Errt/Epsi/Epsi ! 1.414213562 = sqrt(2)
         ELSE
            IF ( mi.EQ.0 ) fct = 1./Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.414213562*Errt/Epsi/Epsi ! 1.414213562 = sqrt(2)
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
            STAMP = DCMPLX(cic,cis)
         ELSE
            dwi = (bic2-cis+ca/a)/2.
            exa = (bis2+cic+sa/a)/2.
            STAMP = DCMPLX(dwi,exa)
         ENDIF
         IF ( la.EQ.2 ) THEN
            IF ( mi.EQ.0 ) fct = .75*(3.-Epsi*Epsi)*axi*axi/Epsi/Epsi
            IF ( mi.EQ.1 ) fct = 1.837117307*Errt*axi*axi/Epsi/Epsi ! 1.837117307 = sqrt(27/8)
            IF ( mi.EQ.2 ) fct = -.9185586535*Errt*Errt*axi*axi/Epsi/ ! 0.9185586535 = sqrt(27/32)
     &                           Epsi
         ELSEIF ( la.EQ.3 ) THEN
            fct = -.3535533905*Errt*axi*axi ! 0.3535533907 = sqrt(1/8)
         ELSE
            IF ( mi.EQ.0 ) fct = .5*axi/Epsi
            IF ( mi.EQ.1 ) fct = .3535533907*Errt*axi/Epsi ! 0.3535533907 = sqrt(1/8)
         ENDIF
      ENDIF

      STAMP = STAMP*fct
      STAMP = CONJG(STAMP)
      END
