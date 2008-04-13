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
