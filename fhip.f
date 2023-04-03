
C----------------------------------------------------------------------
C SUBROUTINE FHIP
C
C Called by: GOSIA
C
C Purpose: generates a table of the hyperbolic funcions sinh and cosh for
C later use. Note that these are in steps of \Delta\omega = 0.03. These are
C stored in the common block HIPER.
C
C Uses global variables:
C      CH     - table of cosh values
C      LP12   - number of steps of omega (365)
C      SH     - table of sinh values
C
C LP12 (from common MGN) is the number of values to calculate. This is set to
C 365  in GOSIA, which is the dimension of the arrays.

      SUBROUTINE FHIP
      IMPLICIT NONE
      REAL*8 CH , er , ex , SH , w
      INTEGER*4 j , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 ,
     &          LP4 , LP6 , LP7 , LP8 , LP9
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 ,
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /HIPER / SH(365) , CH(365)

      w = -.03d0
      DO j = 1 , LP12
         w = w + .03d0
         ex = EXP(w)
         er = 1./ex
         SH(j) = (ex-er)/2.
         CH(j) = (ex+er)/2.
      ENDDO
      END
