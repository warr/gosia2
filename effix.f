      SUBROUTINE EFFIX(Ipd,En,Effi)
      IMPLICIT NONE
      REAL*8 ABC , AKAVKA , d , Effi , En , enek , enl , pw , s , t , 
     &       THICK , w , xx , yy
      INTEGER*4 i , Ipd , ixi , j , l , ll , n
      DIMENSION xx(51) , yy(51)
      COMMON /EFCAL / ABC(8,10) , AKAVKA(8,200) , THICK(200,7)
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
               CALL LAGRAN(xx,yy,3,0,enl,t,1,1)
               s = s + EXP(t)*THICK(Ipd,l)
            ENDIF
         ENDDO
      ENDIF
      Effi = EXP(-s)
c FITEFF or GREMLIN check
      IF ( AKAVKA(8,Ipd).LE.-999. ) THEN
C     LEUVEN CALIBRATION
         Effi = AKAVKA(1,Ipd)
         enek = 1000.*En
         w = LOG(enek)
         DO ixi = 1 , 6
            Effi = Effi + AKAVKA(ixi+1,Ipd)*w**ixi
         ENDDO
         Effi = EXP(Effi)
         GOTO 99999
      ELSEIF ( AKAVKA(5,Ipd).GT.0. .AND. AKAVKA(5,Ipd).LT.10. ) THEN
c FITEFF eff. calib. by P.Olbratowski use
c PJN@2000
         w = LOG(En/AKAVKA(5,Ipd))
         pw = AKAVKA(2,Ipd)*w
         IF ( En.LT.AKAVKA(5,Ipd) ) pw = pw + 
     &        w*w*(AKAVKA(3,Ipd)+w*AKAVKA(4,Ipd))
         Effi = Effi*EXP(pw)*AKAVKA(1,Ipd)
         RETURN
      ELSEIF ( AKAVKA(5,Ipd).LT.10. ) THEN
c GREMLIN
         w = LOG(20.*En)
         pw = AKAVKA(1,Ipd) + AKAVKA(2,Ipd)*w + AKAVKA(3,Ipd)
     &        *w*w + AKAVKA(4,Ipd)*w*w*w
         Effi = Effi*EXP(pw)
         IF ( ABS(AKAVKA(5,Ipd)).GE.1.E-9 ) THEN
            n = INT(AKAVKA(6,Ipd)+.1)
            pw = w**n
            w = AKAVKA(5,Ipd)/pw
            Effi = Effi*EXP(w)
         ENDIF
         IF ( ABS(AKAVKA(8,Ipd)).LT.1.E-9 ) RETURN
         w = (AKAVKA(7,Ipd)-1000.*En)/AKAVKA(8,Ipd)
         pw = EXP(w)
         IF ( ABS(pw-1.).LT.1.E-6 ) WRITE (22,99001)
99001    FORMAT (5x,'***** CRASH - EFFIX *****')
         Effi = Effi/(1.-pw)
         RETURN
      ENDIF
c     JAERI calibration - TC, Nov.2000
      w = LOG(En/.511)
      Effi = EXP(AKAVKA(1,Ipd)+AKAVKA(2,Ipd)
     &       *w-EXP(AKAVKA(3,Ipd)+AKAVKA(4,Ipd)*w))
99999 END