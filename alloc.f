
C----------------------------------------------------------------------
C SUBROUTINE ALLOC
C
C Called by: FTBM, GOSIA
C Calls:     RANGEL
C
C Purpose: to calculate and store the ranges of the integration over omega
C for each multipolarity.
C
C Uses global variables:
C      IRA    - range of omega for each multipolarity needed for accuracy
C      LOCQ   - location of collision function in ZETA array
C      LP14   - maximum length of space for collision functions (4900)
C
C Formal parameters:
C      Accur  - accuracy required
C
C We set up the LOCQ array, which indexes the start of the block of collision
C function coefficients for each omega. For a given multipolarity, lambda,
C there are lambda+1 collision functions, which have to be evaluated for a
C set of different omega values. We don't integrate over all possible omega
C values, but estimate how many values we need to achieve the accuracy Accur.
C The function RANGEL calculates how many we need for each multipolarity and
C stores it in IRA.
C
C Later (in SNAKE) we will store the values of the collision functions for
C 2 * IRA(1) values for E1, 3 * IRA(2) values for E2, 4 * IRA(3) values for
C E3... 3 * IRA(8) values for M2 in that order.
C
C We are limited to a maximum of LP14 (=4900) values in total.

      SUBROUTINE ALLOC(Accur)
      IMPLICIT NONE
      REAL*8 Accur
      INTEGER*4 iflag , IRA , j , k , k1 , load , LOCQ , LP1 , LP10 ,
     &          LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , LP6 ,
     &          LP7 , LP8 , LP9
      INTEGER*4 MAXLA
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 ,
     &                LP10 , LP11 , LP12 , LP13 , LP14

C     Call RANGEL to determine the range of the integration over omega, which
C     depends on the accuracy Accur.
      CALL RANGEL(Accur)

C     First zero all the elements
      load = 0
      iflag = 0
      DO j = 1 , 8
         DO k = 1 , 7
            LOCQ(j,k) = 0
         ENDDO
      ENDDO

C     Now store values for E1...E6
      DO k = 1 , 6
         k1 = k + 1
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k)
         ENDDO
      ENDDO

C     And for M1, M2
      DO k = 7 , 8
         k1 = k - 6
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k) ! IRA(k) is the number of omega values needed for requested accuracy
         ENDDO
      ENDDO

      IF ( load.LE.LP14 ) RETURN ! The Q-functions must fit in the last LP14 words of ZETA

      WRITE (22,99001)
99001 FORMAT (5X,'NO SPACE FOR Q FUNCTIONS TABULATION'//5X,
     &        'SORRY,JOB WILL BE BRUTALLY TERMINATED!')
      STOP 'JOB TERMINATED BY ALLOC' ! Added N. Warr Jul2007
      END
