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
