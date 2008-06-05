 
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
C      Wth    - theta of centre of detector (degrees)
C      Wph    - phi of centre of detector (degrees)
C      Wthh   - half angle subtended (degrees)
C      Naa    - number of theta divisions
C      Ifw    - flag: 0 for meshpoints, 1 for subdivisions, 2 for pin diodes
C      Pfi    - phi range for each theta value
C      Wpi    - phi range of detector
C      Wtlb   - angle of particle detector in theta (degrees) in lab frame
C      Lz     - experiment number
C      Tyy    - lower limit of theta (degrees)
C      Tzz    - upper limit of theta (degrees)
 
      SUBROUTINE COORD(Wth,Wph,Wthh,Naa,Ifw,Pfi,Wpi,Wtlb,Lz,Tyy,Tzz)
      IMPLICIT NONE
      REAL*8 ga , gi , Pfi , rade , rmass , TACOS , TASIN , thetb , 
     &       ttcm , Tyy , Tzz
      REAL*8 wpa , Wph , Wpi , ws , Wth , Wthh , Wtlb , xaa , xph , 
     &       xth , xthh , za , za1 , zb , zl
      INTEGER*4 i , Ifw , Lz , Naa
      DIMENSION Pfi(101) , Wpi(100,2)
      INCLUDE 'vlin.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'cx.inc'
      INCLUDE 'seck.inc'
      DATA rade/57.2957795/ ! 180 / pi
      DATA ws/0./

      IF ( Ifw.EQ.0 ) THEN
         Tyy = Wth - Wthh ! Lower limit of theta
         Tzz = Wth + Wthh ! Upper limit of theta
      ENDIF

C     Convert to radians
      xth = Wth/rade ! theta of centre of detector in radians
      xph = Wph/rade ! phi of centre of detector in radians
      xthh = Wthh/rade ! half angle subtended in radians

C     pre-calculate trigonometric functions
      zl = TAN(xthh)
      za = COS(xth)
      za1 = SIN(xth)
      zb = COS(xthh)

      rmass = XA1(Lz)/XA ! Mass ratio for this experiment
      IF ( IZ1(Lz).LT.0 ) rmass = 1./rmass

C     Calculate size of each division (ws)
      IF ( Ifw.NE.2 ) THEN ! Unless we are using the pin diode option
         ws = (Tzz-Tyy)/(Naa+1)
         IF ( Ifw.EQ.1 ) ws = (Tzz-Tyy)/(Naa-1)
      ENDIF

      DO i = 1 , Naa ! Loop over theta divisions
         IF ( Ifw.NE.2 ) THEN ! Not pin diode option
            IF ( Ifw.EQ.0 ) YV(i) = Tyy + i*ws ! theta value for this step in degrees
            xaa = (Tyy+ws*(i-1))/rade ! and in radians
            IF ( Ifw.EQ.1 .AND. (i.EQ.1 .OR. i.EQ.Naa) ) THEN
               Pfi(i) = 0.
               GOTO 100
            ELSE
               IF ( Ifw.EQ.0 ) xaa = YV(i)/rade
            ENDIF
         ELSE ! Pin diode option
            xaa = ABS(Wtlb)/rade ! Detector angle theta in lab frame in radians
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
            ENDIF ! End of pin diode option
 20         ttcm = xaa + TASIN(rmass*SIN(xaa))
            xaa = (3.14159265-ttcm)/2.
         ENDIF

 50      gi = (za-COS(xaa)/zb)/(zl*za1)
         ga = TACOS(gi)
         wpa = ATAN(zl*SIN(ga)/(za1+zl*COS(ga)*za))
         wpa = ABS(wpa)
         IF ( Ifw.EQ.2 ) THEN ! Pin diode option
            FIEX(Lz,1) = (xph-wpa) ! phi min
            FIEX(Lz,2) = (xph+wpa) ! phi max
         ELSEIF ( Ifw.EQ.1 ) THEN ! Interpolation option
            Pfi(i) = 2.*wpa*rade
         ELSE ! Meshpoint option
            Wpi(i,1) = (xph-wpa)*rade ! Lower phi limit
            Wpi(i,2) = (xph+wpa)*rade ! Upper phi limit
         ENDIF
 100  ENDDO ! Loop on theta divisions i

C     If a negative value of theta was specified for a meshpoint value,
C     we use the target angle
      IF ( Wtlb.LT.0. .AND. Ifw.EQ.0 ) THEN
         DO i = 1 , Naa ! For each theta division
            xaa = YV(i)/rade ! theta in radians
            thetb = ATAN(SIN(2.*xaa)/(rmass-COS(2.*xaa)))*rade
            IF ( thetb.LT.0. ) thetb = 180. + thetb
            YV(i) = -1.*thetb
            Wpi(i,1) = Wpi(i,1) + 180.
            Wpi(i,2) = Wpi(i,2) + 180.
         ENDDO
      ENDIF
      END
