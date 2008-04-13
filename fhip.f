      SUBROUTINE FHIP
      IMPLICIT NONE
      REAL*8 CH , er , ex , SH , w
      INTEGER*4 j , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , 
     &          LP4 , LP6 , LP7 , LP8 , LP9
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /HIPER / SH(365) , CH(365)
      w = -.03
      DO j = 1 , LP12
         w = w + .03
         ex = EXP(w)
         er = 1./ex
         SH(j) = (ex-er)/2.
         CH(j) = (ex+er)/2.
      ENDDO
      END
