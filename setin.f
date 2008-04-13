      SUBROUTINE SETIN
      IMPLICIT NONE
      REAL*8 ADB , CH , EPS , EROOT , FIEX , SH
      INTEGER*4 IAXS , IEXP , k , LP1 , LP10 , LP11 , LP12 , LP13 , 
     &          LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /HIPER / SH(365) , CH(365)
      COMMON /ADX   / ADB(365)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      DO k = 1 , LP12
         ADB(k) = EPS(IEXP)*SH(k) + .03*(k-1)
      ENDDO
      END
