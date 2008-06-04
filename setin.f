 
C----------------------------------------------------------------------
C SUBROUTINE SETIN
C
C Called by FTBM, GOSIA
C
C Purpose: calculate the adiabatic parameter:
C \epsilon \sinh(\omega) + \omega
C
C Uses global variables:
C      ADB    - adiabatic function
C      EPS    - epsilon
C      IEXP   - experiment number
C      LP12   - number of steps of omega (365)
C      SH     - table of sinh values
C
C Note that it uses the tables of sinh calculated by FHIP (SH in common
C block HIPER) and that both the sinh table and this table of the adiabatic
C parameter are in steps of \Delta\omega = 0.03. The resulting table of
C adiabatic parameters are stored in ADB in common block ADX.
C
C LP12 (from common MGN) is the number of values to calculate. This is set to
C 365  in GOSIA, which is the dimension of the array.
 
      SUBROUTINE SETIN
      IMPLICIT NONE
      REAL*8 ADB , CH , SH
      INTEGER*4 k , LP1 , LP10 , LP11 , LP12 , LP13 , 
     &          LP14 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /HIPER / SH(365) , CH(365)
      COMMON /ADX   / ADB(365)
      INCLUDE 'kin.inc'
      
      DO k = 1 , LP12
         ADB(k) = EPS(IEXP)*SH(k) + .03*(k-1)
      ENDDO
      END
