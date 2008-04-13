      COMPLEX*16 FUNCTION EXPON(Inx,Npt,Isg,Isg1,Ndiv,Kdiv)
      IMPLICIT NONE
      REAL*8 ADB , XI
      INTEGER*4 Inx , Isg , Isg1 , Kdiv , Ndiv , Npt
      COMPLEX*16 expo1 , ci , expox , TCEXP
      COMMON /ADX   / ADB(365)
      COMMON /CXI   / XI(500)
      DATA ci/(0.,1.)/
      expox = TCEXP(ci*XI(Inx)*ADB(Npt)*Isg)
      EXPON = expox
      IF ( Ndiv.NE.0 ) THEN
         expo1 = TCEXP(ci*XI(Inx)*ADB(Npt+Isg1)*Isg)
         EXPON = expox + DBLE(Kdiv)*(expo1-expox)/DBLE(Ndiv)
      ENDIF
      END
