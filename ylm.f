      SUBROUTINE YLM(Theta,Ylmr)
      IMPLICIT NONE
      REAL*8 ct , ctsq , EPS , EROOT , FIEX , st , Theta , Ylmr
      INTEGER*4 i , IAXS , IEXP , j , l , lf , m
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      DIMENSION Ylmr(9,9) , st(7)
      ct = COS(Theta)
      ctsq = ct*ct
      IF ( IAXS(IEXP).EQ.0 ) THEN
         Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
         Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
         Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
         GOTO 99999
      ENDIF
      st(1) = SIN(Theta)
      DO i = 2 , 7
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      Ylmr(1,3) = .1089659406
      Ylmr(1,2) = -.2179318812*ct
      Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
      Ylmr(2,5) = .1248361677
      Ylmr(2,4) = -.3530900028*ct
      Ylmr(2,3) = .0943672726*(7.*ctsq-1.)
      Ylmr(2,2) = -.1334554768*ct*(7.*ctsq-3.)
      Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
      Ylmr(3,7) = .1362755124
      Ylmr(3,6) = -.4720722226*ct
      Ylmr(3,5) = .100646136*(11.*ctsq-1.)
      Ylmr(3,4) = -.1837538634*ct*(11.*ctsq-3.)
      Ylmr(3,3) = .0918769316*((33.*ctsq-18.)*ctsq+1.)
      Ylmr(3,2) = -.1162161475*ct*((33.*ctsq-30.)*ctsq+5.)
      Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
      DO l = 1 , 3
         lf = 2*l + 1
         DO m = 2 , lf
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)
         ENDDO
      ENDDO
99999 END