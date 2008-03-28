 
C----------------------------------------------------------------------
C SUBROUTINE RESET
C
C Called by: INTG
C
C Purpose: to advance by one step. This means f(n-3) is set to the old value
C          of f(n-2), f(n-2) is set to the old value of f(n-1) and f(n-1) is
C          set to the old value of f(n).
C
C Uses global variables:
C      ARM    - reduced matrix elements
C      CAT    -
C      ISMAX  -
C      NMAX   - number of levels
C      NSTART -
C
C Formal parameters:
C      Iso    -
 
      SUBROUTINE RESET(Iso)
      IMPLICIT NONE
      REAL*8 CAT
      INTEGER*4 ir , ISMAX , Iso , j , NDIM , NMAX , NMAX1 , NSTART , 
     &          NSTOP
      COMPLEX*16 ARM
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /AZ    / ARM(600,7)
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      
      IF ( Iso.EQ.0 ) THEN
         DO j = 1 , NMAX
            ir = NSTART(j) - 1
 20         ir = ir + 1
            ARM(ir,1) = ARM(ir,2)
            ARM(ir,2) = ARM(ir,3)
            ARM(ir,3) = ARM(ir,4)
            IF ( CAT(ir,3).LT.-.1 ) GOTO 20
         ENDDO
         GOTO 99999
      ENDIF
       
      DO j = 1 , ISMAX
         ARM(j,1) = ARM(j,2)
         ARM(j,2) = ARM(j,3)
         ARM(j,3) = ARM(j,4)
      ENDDO
99999 END
