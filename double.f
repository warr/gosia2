 
C----------------------------------------------------------------------
C SUBROUTINE DOUBLE
C
C Called by: INTG
C
C Purpose: to double the step size for the integeration in INTG.
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      CAT    - substates of levels (n_level, J, m)
C      ISMAX  - number of substates used
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C
C Formal parameters:
C      Iso    - isotropic flag
 
      SUBROUTINE DOUBLE(Iso)
      IMPLICIT NONE
      INTEGER*4 ir , Iso , j
      COMPLEX*16 ARM , fpom
      INCLUDE 'cexc0.inc'
      COMMON /AZ    / ARM(600,7)
      INCLUDE 'clcom8.inc'
      INCLUDE 'coex2.inc'
      
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX ! Loop over levels
            ir = NSTART(j) - 1 ! Index of first substate of level - 1
 20         ir = ir + 1
            fpom = ARM(ir,2)
            ARM(ir,2) = -8.*ARM(ir,3) + 6.*ARM(ir,2) + 3.*ARM(ir,4)
            ARM(ir,1) = -16.*ARM(ir,1) + 9.*ARM(ir,2) + 9.*fpom - 
     &                  ARM(ir,4)
            ARM(ir,3) = fpom
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
       
      DO j = 1 , ISMAX ! Loop over substates
         fpom = ARM(j,2)
         ARM(j,2) = -8.*ARM(j,3) + 6.*ARM(j,2) + 3.*ARM(j,4)
         ARM(j,1) = -16.*ARM(j,1) + 9.*ARM(j,2) + 9.*fpom - ARM(j,4)
         ARM(j,3) = fpom
      ENDDO
99999 END
