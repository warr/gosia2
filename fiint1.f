
C----------------------------------------------------------------------
C SUBROUTINE FIINT1
C
C Called by: ANGULA
C
C Purpose: integrate over phi in lab frame
C
C Formal parameters:
C      Fi0    - phi_0
C      Fi1    - phi_1
C      Alab   - return value
C      Ixs    - axial symmetry flag

      SUBROUTINE FIINT1(Fi0,Fi1,Alab,Ixs)
      IMPLICIT NONE
      REAL*8 Alab , Fi0 , Fi1 , wsp
      INTEGER*4 Ixs , j , m , mm
      DIMENSION Alab(9,9)

      IF ( Ixs.NE.0 ) THEN
         DO m = 2 , 9
            mm = m - 1
            wsp = (SIN(mm*Fi1)-SIN(mm*Fi0))/mm
            DO j = 1 , 9
               Alab(j,m) = Alab(j,m)*wsp
            ENDDO
         ENDDO
         wsp = Fi1 - Fi0
      ENDIF
      IF ( Ixs.EQ.0 ) wsp = 6.283185308 ! 6.283185308 = 2 * pi
      DO j = 1 , 9
         Alab(j,1) = Alab(j,1)*wsp
      ENDDO
      END
