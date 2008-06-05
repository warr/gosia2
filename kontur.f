 
C----------------------------------------------------------------------
C SUBROUTINE KONTUR
C
C Called by: GOSIA
C Calls:     FTBM, LIMITS, RK4
C
C Purpose:
C
C Uses global variables:
C      DEVU   -
C      ELM    - matrix elements
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      HLM    - previous values of matrix elements
C      INTR   - flag to swap chisqr and log(chisqr)
C      IPS1   - terminate after calculating and writing correction factors
C      LNY    - use logs to calculate chi squared
C      MEMAX  - number of matrix elements
C      NWR    - number of datapoints used in fit
C      SA     - ratio of elements for correlated elements
C      XV     - energy meshpoints where we calculate exact Coulex
C      YV     - scattering angle meshpoints where we calculate exact Coulex
C
C Formal parameters:
C      Idr    - number of decays
C      Chis0  -
C      Chil   -
C      Ifbf   -
C      Inpo   -
C      Jj     - matrix element
C      Sh     -
C      Bten   -
C      Rem    - natural log of the largest value the computer can represent
 
      SUBROUTINE KONTUR(Idr,Chis0,Chil,Ifbf,Inpo,Jj,Sh,Bten,Rem)
      IMPLICIT NONE
      REAL*8 ac , Bten , c , Chil , chilo , Chis0 , chis1 , chis2 , d1 , 
     &       d2 , ELM , ELML , ELMU , f , h
      REAL*8 Rem , RK4 , SA , sajj , Sh , t , v , ww , x , y
      INTEGER*4 i , Idr , Ifbf , Inpo , INTR , IPS1 , itl , IVAR , ix , 
     &          j , Jj , l , LMAXE , LNY , m , MAGEXC , MEMAX , MEMX6 , 
     &          NWR
      DIMENSION f(3) , Bten(1200)
      INCLUDE 'vlin.inc'
      INCLUDE 'dftb.inc'
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      INCLUDE 'hhh.inc'
      COMMON /ILEWY / NWR
      COMMON /LOGY  / LNY , INTR , IPS1

      LNY = 0
      h = .05*ABS(HLM(Jj))
      IF ( Inpo.NE.-1 ) h = ABS(Sh)
 100  INTR = 0
      sajj = ABS(SA(Jj)) ! ratio of matrix elements for correlation
      DO l = 1 , MEMAX ! For each matrix element
         ELM(l) = HLM(l)
         SA(l) = SA(l)/sajj
      ENDDO
      YV(1) = 0.
      XV(1) = HLM(Jj)
      f(3) = 1.
      i = 1
 200  itl = 0
      v = ELMU(Jj) - ELM(Jj)
      IF ( SA(Jj).LT.0. ) v = ELM(Jj) - ELML(Jj)
      IF ( h.GT.v ) itl = 1
      IF ( h.GT.v ) h = v
      i = i + 1
      f(1) = f(3)
      DO j = 1 , MEMAX
         ELM(j) = .5*h*SA(j) + ELM(j)
      ENDDO
      CALL LIMITS ! Constrain matrix elements within limits
      CALL FTBM(3,chis1,Idr,1,chilo,Bten)
      IF ( chis1.LE.Chis0 ) THEN
         IF ( Inpo.EQ.-1 ) WRITE (22,99003) Jj , ELM(Jj) , chis1
         IF ( chis1.LE.Chil .AND. Inpo.NE.-1 ) THEN
            Ifbf = 1
            ix = 1
            Chil = chis1
            WRITE (22,99004) Chil
            GOTO 500
         ENDIF
      ENDIF
 300  ww = .5*(Chis0-chis1)*NWR
      IF ( ww.GE.Rem ) GOTO 700
      f(2) = EXP(ww)
      IF ( i.EQ.2 .AND. f(2).LT..1 .AND. ABS(XV(1)-HLM(Jj)).LT.1E-9 )
     &     THEN
         h = h/2.
         GOTO 100
      ELSE
         DO j = 1 , MEMAX ! For each matrix element
            ELM(j) = ELM(j) + .5*SA(j)*h
         ENDDO
         v = ELM(Jj)
         CALL LIMITS ! Constrain matrix elements within limits
         IF ( ABS(v-ELM(Jj)).GT.1.E-6 ) itl = 1
         CALL FTBM(3,chis2,Idr,1,chilo,Bten)
         IF ( chis2.LE.Chis0 ) THEN
            IF ( Inpo.EQ.-1 ) WRITE (22,99003) Jj , ELM(Jj) , chis2
            IF ( chis2.LE.Chil .AND. Inpo.NE.-1 ) THEN
               Ifbf = 1
               ix = 2
               Chil = chis2
               WRITE (22,99004) Chil
               GOTO 500
            ENDIF
         ENDIF
      ENDIF
 400  ww = .5*(Chis0-chis2)*NWR
      IF ( ww.GT.Rem ) GOTO 700
      f(3) = EXP(ww)
      IF ( itl.EQ.1 ) WRITE (22,99001) Jj
99001 FORMAT (5X,'WARNING-ME(',1I3,')',5X,
     &        'INTEGRATION STOPPED AT THE LIMIT')
      IF ( i.EQ.2 ) THEN
         IF ( itl.NE.1 ) THEN
            IF ( f(3).LT..1 .AND. ABS(XV(1)-HLM(Jj)).LT.1.E-9 ) THEN
               h = h/2.
               GOTO 100
            ELSEIF ( f(1).LE.f(2) .OR. f(2).LE.f(3) ) THEN
               IF ( f(1).LT.f(2) .AND. f(2).GT.f(3) ) THEN
                  d1 = f(2) - f(1)
                  d2 = f(3) - f(1)
                  ac = (d2-4.*d1)*h/(d2-2.*d1)/4.
                  DO l = 1 , MEMAX
                     ELM(l) = (ELM(l)-h*SA(l)) + ac*SA(l)
                  ENDDO
                  CALL LIMITS
                  XV(1) = ELM(Jj)
                  i = 1
                  CALL FTBM(3,chis1,Idr,1,chilo,Bten)
                  ww = .5*(Chis0-chis1)*NWR
                  IF ( ww.GE.Rem ) GOTO 700
                  f(3) = EXP(ww)
                  GOTO 200
               ELSE
                  i = 1
                  XV(1) = ELM(Jj)
                  IF ( Inpo.EQ.-1 ) h = 2.*h
                  GOTO 200
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      y = YV(i-1)
      YV(i) = RK4(y,h,f)
      XV(i) = ELM(Jj)
      IF ( NWR*(chis2-Chis0).LT.2. .AND. Inpo.EQ.-1 ) h = 2.*h
      IF ( itl.EQ.1 ) GOTO 600
      IF ( f(3).GE.1.E-3 ) GOTO 200
      GOTO 600
 500  REWIND 17
      DO l = 1 , MEMAX ! For each matrix element
         WRITE (17,*) ELM(l)
      ENDDO
      IF ( ix.EQ.1 ) GOTO 300
      IF ( ix.NE.2 ) GOTO 200
      GOTO 400
 600  c = YV(i)
      m = 0
      DO l = 1 , i
         YV(l) = 1.00001 - YV(l)/c
         IF ( m.EQ.0 .AND. YV(l).LT..317 ) m = l
      ENDDO
      x = (XV(m)-XV(m-1))*(.317-YV(m))/(YV(m-1)-YV(m))
      t = XV(m) - x - HLM(Jj)
      IF ( t.GE.0. ) DEVU(Jj) = t
      IF ( t.LT.0. ) DEVD(Jj) = t
      RETURN
 700  WRITE (22,99002) Jj
99002 FORMAT (5X,'** WARNING **',/,2X,'ME=',1I3,2X,
     &     'TOO FAR FROM THE MINIMUM TO CARRY OUT THE ERROR ESTIMATION!'
     &     ,/)
99003 FORMAT (5X,'ELM(',1I3,')=',1F10.6,5X,'CHISQ=',1E12.4)
99004 FORMAT (10X,'BETTER POINT FOUND...MATRIX ELEMENTS WRITTEN ON 17',
     &        3X,'CHISQ=',1E12.4)
      END
