 
C----------------------------------------------------------------------
C SUBROUTINE PATH
C
C Called by: GOSIA
C
C Purpose: Calculate path for each level
C
C Uses global variables:
C      CAT    - substates of levels (n_level, J, m)
C      IPATH  - index of substate in level with same m as substate Irld
C      NMAX   - number of levels
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C
C Formal parameters:
C      Irld   - index into ARM array
 
      SUBROUTINE PATH(Irld)
      IMPLICIT NONE
      REAL*8 CAT , spm , vl
      INTEGER*4 i , IPATH , Irld , ISMAX , isp , ist , j , MAGA , 
     &          NSTART , NSTOP
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      INCLUDE 'coex2.inc'
      COMMON /CLCOM8/ CAT(600,3) , ISMAX

      spm = CAT(Irld,3) ! m quantum number for substate Irld
      DO i = 2 , NMAX ! For each level except ground state
         IPATH(i) = 0
         ist = NSTART(i) ! Index of first substate for level
         IF ( ist.NE.0 ) THEN ! If this is non-zero
            isp = NSTOP(i) ! Index of last substate for level
            DO j = ist , isp ! For each substate of level
               vl = CAT(j,3) ! m quantum number for substate j
               IF ( ABS(vl-spm).LT.1.E-6 ) GOTO 50 ! Jump if they have the same m
            ENDDO
         ENDIF
         GOTO 100
 50      IPATH(i) = j ! Store it
 100  ENDDO
      IPATH(1) = Irld ! Special case of ground state
      END
