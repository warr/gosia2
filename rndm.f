
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
C
C This function uses the Lehmer method

      REAL*8 FUNCTION RNDM(Se)
      IMPLICIT NONE
      REAL*8 a, b, Se

      If ( Se.LE.0 ) Se = 1.0D0
      a = 16807.D0
      b = 2147483647.D0
      Se = MOD(a * Se, b)
      RNDM = MOD(a * Se, b) / b
      END
