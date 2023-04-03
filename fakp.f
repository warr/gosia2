
C----------------------------------------------------------------------
C SUBROUTINE FAKP
C
C Called by: GOSIA
C Calls:     PRIM
C
C Purpose: calculate log of primes and the factoring of primes
C
C Uses global variables:
C      IP     - table of primes
C      IPI    - number of time a number is divisible by each prime in IP
C      KF     - sum of factors of primes
C      PILOG  - table of natural logs of primes

      SUBROUTINE FAKP
      IMPLICIT NONE
      INTEGER*4 i , IP , IPI , k , KF , l
      REAL*8 PILOG , x
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)

C     Calculate log of primes
      DO i = 1 , 26
         x = DBLE(IP(i))
         PILOG(i) = LOG(x)
      ENDDO

C     Initialise KF
      DO l = 1 , 26
         KF(1,l) = 0
         KF(2,l) = 0
      ENDDO

C     Calculate factors of numbers
      DO k = 3 , 101
         CALL PRIM(k-1) ! Puts factors in IPI array
         DO i = 1 , 26
            KF(k,i) = KF(k-1,i) + IPI(i) ! IPI are number of times prime is a factor
         ENDDO
      ENDDO
      END
