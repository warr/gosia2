 
C----------------------------------------------------------------------
C SUBROUTINE COORD
C
C Called by: GOSIA
C Calls:     TACOS, TASIN
C
C Purpose: calculate geometry for circular detector
C
C Uses global variables:
C      FIEX   - phi range of particle detector
C      ISKIN  - kinematic flag
C      IZ1    - Z of non-investigated nucleus
C      XA     - A of investigated nucleus
C      XA1    - A of non-investigated nucleus
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C
C Formal parameters:
C      Wth    - theta of centre of detector
C      Wph    - phi of centre of detector
C      Wthh   - half angle subtended
C      Naa    - number of mesh points
C      Ifw    - flag
C      Pfi    -
C      Wpi    -
C      Wtlb   -
C      Lz     - experiment number
C      Tyy    -
C      Tzz    -
 
      SUBROUTINE COORD(Wth,Wph,Wthh,Naa,Ifw,Pfi,Wpi,Wtlb,Lz,Tyy,Tzz)
      IMPLICIT NONE
      REAL*8 DS , DSE , DSG , EP , EPS , EROOT , FIEX , ga , gi , Pfi , 
     &       rade , rmass , TACOS , TASIN , thetb , TLBDG , ttcm , Tyy , 
     &       Tzz , VINF
      REAL*8 wpa , Wph , Wpi , ws , Wth , Wthh , Wtlb , XA , XA1 , xaa , 
     &       xph , xth , xthh , XV , YV , za , za1 , zb , zl , ZV
      INTEGER*4 i , IAXS , IEXP , Ifw , ISKIN , IZ , IZ1 , Lz , Naa , 
     &          NEXPT
      DIMENSION Pfi(101) , Wpi(11,2)
      COMMON /VLIN  / XV(51) , YV(51) , ZV(20) , DSG(20) , DSE(20) , DS
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /SECK  / ISKIN(50)
      DATA rade/57.2957795/ ! pi / 2

      IF ( Ifw.EQ.0 ) THEN
         Tyy = Wth - Wthh
         Tzz = Wth + Wthh
      ENDIF
      xth = Wth/rade
      xph = Wph/rade
      xthh = Wthh/rade
      zl = TAN(xthh)
      za = COS(xth)
      za1 = SIN(xth)
      zb = COS(xthh)
      rmass = XA1(Lz)/XA
      IF ( IZ1(Lz).LT.0 ) rmass = 1./rmass
      IF ( Ifw.NE.2 ) THEN
         ws = (Tzz-Tyy)/(Naa+1)
         IF ( Ifw.EQ.1 ) ws = (Tzz-Tyy)/(Naa-1)
      ENDIF
      DO i = 1 , Naa
         IF ( Ifw.NE.2 ) THEN
            IF ( Ifw.EQ.0 ) YV(i) = Tyy + i*ws
            xaa = (Tyy+ws*(i-1))/rade
            IF ( Ifw.EQ.1 .AND. (i.EQ.1 .OR. i.EQ.Naa) ) THEN
               Pfi(i) = 0.
               GOTO 100
            ELSE
               IF ( Ifw.EQ.0 ) xaa = YV(i)/rade
            ENDIF
         ELSE
            xaa = ABS(Wtlb)/rade
            IF ( Wtlb.GT.0. ) GOTO 50
            IF ( IZ1(Lz).LT.0 ) THEN
               IF ( XA.LE.XA1(Lz) ) GOTO 20
            ELSEIF ( XA1(Lz).LE.XA ) THEN
               GOTO 20
            ENDIF
            IF ( ISKIN(Lz).EQ.0 ) THEN
               ttcm = xaa - TASIN(rmass*SIN(xaa))
               xaa = ABS(ttcm)/2.
               GOTO 50
            ENDIF
 20         ttcm = xaa + TASIN(rmass*SIN(xaa))
            xaa = (3.14159265-ttcm)/2.
         ENDIF
 50      gi = (za-COS(xaa)/zb)/(zl*za1)
         ga = TACOS(gi)
         wpa = ATAN(zl*SIN(ga)/(za1+zl*COS(ga)*za))
         wpa = ABS(wpa)
         IF ( Ifw.EQ.2 ) THEN
            FIEX(Lz,1) = (xph-wpa)
            FIEX(Lz,2) = (xph+wpa)
         ELSEIF ( Ifw.EQ.1 ) THEN
            Pfi(i) = 2.*wpa*rade
         ELSE
            Wpi(i,1) = (xph-wpa)*rade
            Wpi(i,2) = (xph+wpa)*rade
         ENDIF
 100  ENDDO
      IF ( Wtlb.LT.0. .AND. Ifw.EQ.0 ) THEN
         DO i = 1 , Naa
            xaa = YV(i)/rade
            thetb = ATAN(SIN(2.*xaa)/(rmass-COS(2.*xaa)))*rade
            IF ( thetb.LT.0. ) thetb = 180. + thetb
            YV(i) = -1.*thetb
            Wpi(i,1) = Wpi(i,1) + 180.
            Wpi(i,2) = Wpi(i,2) + 180.
         ENDDO
      ENDIF
      END
