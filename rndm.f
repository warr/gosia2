      FUNCTION RNDM(Se)
      IMPLICIT NONE
      REAL*8 ai , p , r , RNDM , rxdm , Se , t , u
      INTEGER*4 i
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
