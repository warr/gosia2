
C-----------------------------------------------------------------------
C FUNCTION EXPON
C
C Called by: NEWLV
C Calls:     TCEXP
C
C Purpose: calculates the exponential:
C       exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C
C Uses global variables:
C      EXPO   - adiabatic exponential
C      ADB    - adiabatic function
C      XI     - xi coupling coefficients
C
C Formal parameters:
C      Inx    - index in XI array
C      Npt    - index in ADB array (this is omega / DOMEGA)
C      Isg    - phase
C      Isg1   - index for sigma
C      Ndiv   - number of divisions
C      Kdiv   - index for division
C
C Return value:
C      the exponential
C
C ci is sqrt(-1)
C XI (from common block CXI) are the XI coupling constants calculated in
C the function LOAD.
C ADB is the adiabatic parameters \epsilon \sinh(\omega) + \omega calculated
C in the function SETIN.

      COMPLEX*16 FUNCTION EXPON(Inx,Npt,Isg,Isg1,Ndiv,Kdiv)
      IMPLICIT NONE
      INTEGER*4 Inx , Isg , Isg1 , Kdiv , Ndiv , Npt
      COMPLEX*16 expo1 , ci , expox , TCEXP
      INCLUDE 'adx.inc'
      INCLUDE 'cxi.inc'
      DATA ci/(0.,1.)/
      
      expox = TCEXP(ci*XI(Inx)*ADB(Npt)*Isg)
      EXPON = expox
      IF ( Ndiv.NE.0 ) THEN
         expo1 = TCEXP(ci*XI(Inx)*ADB(Npt+Isg1)*Isg)
         EXPON = expox + DBLE(Kdiv)*(expo1-expox)/DBLE(Ndiv)
      ENDIF
      END
