
C----------------------------------------------------------------------
C SUBROUTINE LIMITS
C
C Called by: KONTUR, MINI, PRELM
C
C Purpose: to constrain the matrix elements to within the limits specified by
C the user.
C
C Uses global variables:
C      ELM    - matrix elements
C      ELML   - lower limit on matrix elements
C      ELMU   - upper limit on matrix elements
C      IVAR   - indicates a limit or correlation is set
C      MEMAX  - number of matrix elements

      SUBROUTINE LIMITS
      IMPLICIT NONE
      INTEGER*4 j
      INCLUDE 'cexc.inc'
      INCLUDE 'comme.inc'

      DO j = 1 , MEMAX ! Loop over matrix elements
         IF ( IVAR(j).NE.0 ) THEN ! If not fixed
            IF ( ELM(j).GT.ELMU(j) .OR. ELM(j).LT.ELML(j) ) THEN
               IF ( ELM(j).GT.ELMU(j) ) THEN
                  ELM(j) = ELMU(j)
                  WRITE (22,99001) j , ELM(j)
               ELSE
                  ELM(j) = ELML(j)
                  WRITE (22,99001) j , ELM(j)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
99001 FORMAT (2X,'Warning - matrix element ',1I3,' reset to ',1F10.6)
      END
