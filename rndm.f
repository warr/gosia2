 
C----------------------------------------------------------------------
C FUNCTION RNDM
C
C Called by: MIXUP
C
C Purpose: Generate a pseudo-random number based on the seed Se
C
C Formal parameters:
C      Se     - seed for random number
C
C It is used to generate random matrix elements as a starting position,
C when OP,RAND is called. The parameter to OP,RAND is the seed here.
 
      REAL*8 FUNCTION RNDM(Se)
      IMPLICIT NONE
      REAL*8 ai , p , r , rxdm , Se , t , u
      INTEGER*4 i
      DATA t/0./
      SAVE t

      IF ( Se.GT.32000. ) Se = 100.*t + .511
      Se = Se*Se
      u = LOG10(Se)
      i = INT(u) + 1
      t = Se/(10.**i)
      r = SQRT(SQRT(SQRT(t)))
      p = SQRT(SQRT(SQRT(.1)))
      rxdm = (r-p)/(1.-p)
      rxdm = 10.*rxdm
      ai = DBLE(INT(rxdm))
      RNDM = rxdm - ai
      END
