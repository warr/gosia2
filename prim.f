 
C----------------------------------------------------------------------
C SUBROUTINE PRIM
C
C Called by: FAKP
C
C Purpose: Find out how many times each prime divides a number N
C
C Uses global variables:
C      IP     - table of primes
C      IPI    - multipliers for each prime
C
C Formal parameters:
C      N      - number N
 
      SUBROUTINE PRIM(N)
      IMPLICIT NONE
      INTEGER*4 i , N , nni , nnk
      INCLUDE 'fakul.inc'

      nnk = N
      DO i = 1 , 26
         nni = nnk
         IPI(i) = 0
 50      nni = nni/IP(i)
         IF ( IP(i)*nni.EQ.nnk ) THEN
            IPI(i) = IPI(i) + 1
            nnk = nni
            GOTO 50
         ENDIF
      ENDDO
      END
