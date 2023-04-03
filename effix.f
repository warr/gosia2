
C----------------------------------------------------------------------
C SUBROUTINE EFFIX
C
C Called by: CEGRY, GOSIA2
C Calls:     LAGRAN, SPLNER
C
C Purpose: calculate the efficiency of the detector at a given energy.
C
C Uses global variables:
C      ABC    - absorption coefficients
C      AKAVKA - efficiency curve parameters
C      THICK  - thickness of each absorber type
C
C Formal parameters:
C      Ipd    - detector number
C      En     - gamma-ray energy
C      Effi   - efficiency
C
C Note that it uses LAGRAN or SPLNER according to the ISPL flag to
C interpolate between the data points given by the user.
C
C The efficiency curve parameters are those of GREMLIN plus an extra control
C flag:
C     AKAVKA(1) = a0
C     AKAVKA(2) = a1
C     AKAVKA(3) = a2
C     AKAVKA(4) = a3
C     AKAVKA(5) = f - for F-factor
C     AKAVKA(6) = N - for F-factor
C     AKAVKA(7) = b - for Woods-saxon factor
C     AKAVKA(8) = c - for Woods-saxon factor
C     AKAVKA(9) = control flag
C
C Efficiency parametrizations (control flag):
C     0  - Gremlin
C     1  - Jaeri
C     2  - Fiteff
C     3  - Leuven
C     4  - Radware

      SUBROUTINE EFFIX(Ipd,En,Effi)
      IMPLICIT NONE
      REAL*8 d , Effi , En , enl , pw , s , t , w , xx , yy
      INTEGER*4 i , Ipd , j , l , ll , n
      DIMENSION xx(101) , yy(101)
      INCLUDE 'efcal.inc'
      INCLUDE 'ccc.inc'

      Effi = 1.E-6
      En = En + 1.E-24
      enl = LOG(En)
      DO i = 1 , 10
         ll = 11 - i
         j = ll
         IF ( enl.GE.ABC(8,ll) ) GOTO 100
         j = -1
      ENDDO
 100  IF ( j.EQ.-1 ) Effi = 1.E-10
      IF ( j.EQ.-1 ) RETURN
      IF ( j.EQ.1 .OR. j.EQ.10 ) THEN
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THICK(Ipd,l)).GE.1.E-9 ) THEN
               t = EXP(ABC(l,j))
               d = THICK(Ipd,l)
               s = s + t*d
            ENDIF
         ENDDO
      ELSE
         IF ( j.EQ.9 ) THEN
            xx(1) = ABC(8,8)
            xx(2) = ABC(8,9)
            xx(3) = ABC(8,10)
         ELSE
            xx(1) = ABC(8,j)
            xx(2) = ABC(8,j+1)
            xx(3) = ABC(8,j+2)
         ENDIF
         s = 0.
         DO l = 1 , 7
            IF ( ABS(THICK(Ipd,l)).GE.1.E-9 ) THEN
               IF ( j.EQ.9 ) THEN
                  yy(1) = ABC(l,8)
                  yy(2) = ABC(l,9)
                  yy(3) = ABC(l,10)
               ELSE
                  yy(1) = ABC(l,j)
                  yy(2) = ABC(l,j+1)
                  yy(3) = ABC(l,j+2)
               ENDIF
               IF ( ISPL.EQ.0 ) CALL LAGRAN(xx,yy,3,0,enl,t,1,1)
               IF ( ISPL.EQ.1 ) CALL SPLNER(xx,yy,3,enl,t,1)
               s = s + EXP(t)*THICK(Ipd,l)
            ENDIF
         ENDDO
      ENDIF
      Effi = EXP(-s)

C     Branch according to type of calibration
      IF ( (AKAVKA(8,Ipd).LE.-999.) .OR. (AKAVKA(9,Ipd).EQ.3.) ) THEN
         GOTO 1003 ! Leuven
      ELSEIF ( AKAVKA(9,Ipd).EQ.4. ) THEN
         GOTO 1004 ! Radware
      ELSEIF ( (AKAVKA(5,Ipd).GT.0. .AND. AKAVKA(5,Ipd).LT.10.) .OR.
     &         (AKAVKA(9,Ipd).EQ.2.) ) THEN
         GOTO 1002 ! Fiteff
      ELSEIF ( (AKAVKA(5,Ipd).LT.10.) .AND. (AKAVKA(9,Ipd).NE.1.) ) THEN
         GOTO 1000 ! Gremlin
      ENDIF
      GOTO 1001 ! Jaeri

C-----------------------------------------------------------------
C     GREMLIN efficiency calibration
 1000 w = LOG(20.*En) ! E0 = 50 keV, so w = LOG(En/E0) with En in MeV
      pw = AKAVKA(1,Ipd) + AKAVKA(2,Ipd)*w + AKAVKA(3,Ipd)
     &     *w*w + AKAVKA(4,Ipd)*w*w*w
      Effi = Effi*EXP(pw)
      IF ( ABS(AKAVKA(5,Ipd)).GE.1.E-9 ) THEN ! F-factor
         n = INT(AKAVKA(6,Ipd)+.1)
         pw = w**n
         w = AKAVKA(5,Ipd)/pw
         Effi = Effi*EXP(w)
      ENDIF
      IF ( ABS(AKAVKA(8,Ipd)).LT.1.E-9 ) RETURN
      w = (AKAVKA(7,Ipd)-1000.*En)/AKAVKA(8,Ipd) ! Woods-saxon factor
      pw = EXP(w)
      IF ( ABS(pw-1.).LT.1.E-6 ) WRITE (22,99001)
99001 FORMAT (5x,'***** CRASH - EFFIX *****')
      Effi = Effi/(1.+pw) ! Older versions of gosia have a minus sign here, which is wrong
                          ! because it is not what is done in gremlin (FITFUN) or the gosia manual
      RETURN

C-----------------------------------------------------------------
C     JAERI efficiency calibration - TC, Nov.2000
 1001 w = LOG(En/.511)
      Effi = EXP(AKAVKA(1,Ipd)+AKAVKA(2,Ipd)
     &       *w-EXP(AKAVKA(3,Ipd)+AKAVKA(4,Ipd)*w))
      RETURN

C-----------------------------------------------------------------
C     FITEFF efficiency calibration by P.Olbratowski use
C     PJN@2000
 1002 w = LOG(En/AKAVKA(5,Ipd))
      pw = AKAVKA(2,Ipd)*w
      IF ( En.LT.AKAVKA(5,Ipd) ) pw = pw +
     &     w*w*(AKAVKA(3,Ipd)+w*AKAVKA(4,Ipd))
      Effi = Effi*EXP(pw)*AKAVKA(1,Ipd)
      RETURN

C-----------------------------------------------------------------
C     Leuven efficiency calibration
 1003 Effi = AKAVKA(1,Ipd)
      IF ( AKAVKA(8,Ipd).LE.0.01 ) THEN
        w = LOG(1000.D0*En)
      ELSE
        w = LOG(1000.D0*En/AKAVKA(8,Ipd))
      ENDIF
      DO i = 1 , 6
         Effi = Effi + AKAVKA(i+1,Ipd)*w**i
      ENDDO
      Effi = EXP(Effi)
      RETURN

C-----------------------------------------------------------------
C     Radware efficiency calibration
C     PJN@2008
 1004 w = LOG(En/.1)
      Effi = (AKAVKA(2,Ipd)+(AKAVKA(3,Ipd)+AKAVKA(4,Ipd)*w)*w)
     &       **(-AKAVKA(8,Ipd))
      w = LOG(En)
      Effi = (AKAVKA(5,Ipd)+(AKAVKA(6,Ipd)+AKAVKA(7,Ipd)*w)*w)
     &       **(-AKAVKA(8,Ipd)) + Effi
      Effi = AKAVKA(1,Ipd)*EXP(Effi**(-1/AKAVKA(8,Ipd)))
      RETURN

      END
