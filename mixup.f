
C----------------------------------------------------------------------
C SUBROUTINE MIXUP
C
C Called by: GOSIA
C Calls:     RNDM
C
C Purpose: set the matrix elements to random values, as a starting value.
C
C Uses global variables:
C      ELM    - matrix elements
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      IVAR   - indicates a limit or correlation is set
C      MEMAX  - number of matrix elements
C      SA     - ratio of elements for correlated elements
C      SE     - seed for random number generator
C
C It is called when the user gives the option OP,RAND
C
C Note that if IVAR = 0, then the matrix element is fixed, so we don't do
C anything here. If it is >= 10000, this means it is correlated to another
C matrix element, so use the correlation to determine the new value, which
C may have been changed when we randomized.

      SUBROUTINE MIXUP
      IMPLICIT NONE
      REAL*8 ELM , ELML , ELMU , RNDM , SA , SE
      INTEGER*4 IVAR , k , k1 , LMAXE , MAGEXC , MEMAX , MEMX6
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      COMMON /XRA   / SE

C     Randomize all that are not fixed or correlated
      DO k = 1 , MEMAX ! For each matrix element
         IF ( IVAR(k).NE.0 .AND. IVAR(k).LE.999 ) ! Not fixed or correlated
     &        ELM(k) = ELML(k) + RNDM(SE)*(ELMU(k)-ELML(k))
      ENDDO

C     Now adjust the correlated elements, since we may have changed the
C     element to which it is correlated.
      DO k = 1 , MEMAX
         IF ( IVAR(k).GE.999 ) THEN ! Correlated
            k1 = IVAR(k) - 1000 ! Index to which it is correlated
            IF ( ABS(ELMU(k1)).LT.1.E-9 ) THEN
               ELM(k) = 0.
            ELSE
               ELM(k) = ELM(k1)*SA(k) ! SA is the ratio we require
            ENDIF
         ENDIF
      ENDDO
      END
