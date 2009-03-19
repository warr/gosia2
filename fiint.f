 
C----------------------------------------------------------------------
C SUBROUTINE FIINT
C
C Called by: ANGULA
C
C Purpose: integrate over phi in frame of recoiling nucleus
C
C Formal parameters:
C      Fi0    - phi_0
C      Fi1    - phi_1
C      At     - return value
C      Ixs    - axial symmetry flag

      SUBROUTINE FIINT(Fi0,Fi1,At,Ixs)
      IMPLICIT NONE
      REAL*8 At , Fi0 , Fi1 , wsp
      INTEGER*4 Ixs , j , jf , js , m , mm
      DIMENSION At(28)
      
      IF ( Ixs.NE.0 ) THEN
         DO m = 2 , 7
            js = m/2
            mm = m - 1
            wsp = (SIN(mm*Fi1)-SIN(mm*Fi0))/mm
            js = js*7 + m
            jf = m + 21
            DO j = js , jf , 7
               At(j) = At(j)*wsp
            ENDDO
         ENDDO
         wsp = Fi1 - Fi0
      ENDIF
      IF ( Ixs.EQ.0 ) wsp = 6.283185308 ! 6.283185308 = 2 * pi
      DO j = 1 , 4
         js = (j-1)*7 + 1
         At(js) = At(js)*wsp
      ENDDO
      END
