      FUNCTION RK4(Y,H,F)
      IMPLICIT NONE
      REAL*8 F , H , RK4 , Y
      DIMENSION F(3)
      RK4 = Y + H*(F(1)+4.*F(2)+F(3))/6.
      END
