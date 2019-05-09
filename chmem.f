 
C----------------------------------------------------------------------
C SUBROUTINE CHMEM
C
C Called by: FTBM
C
C Purpose: compare fitted matrix elements with known ones and calculate the
C          effect on the chi squared
C
C Uses global variables:
C      ELM    - matrix elements
C      EAMX   - known matrix elements and their error
C      NAMX   - number of known matrix elements
C      IAMX   - index of matrix element for known matrix element
C      IAMY   - level indices of pair of levels for which matrix element is known
 
      SUBROUTINE CHMEM(Nw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , di
      INTEGER*4 ia , ib , Nw
      INCLUDE 'me2d.inc'
      INCLUDE 'comme.inc'

      IF ( NAMX.EQ.0 ) RETURN
      Nw = Nw + NAMX
      DO ia = 1 , NAMX
         ib = IAMX(ia)
         IF ( IAMY(ia,1).NE.IAMY(ia,2) ) THEN
            di = (ELM(ib)-EAMX(ia,1))
            IF ( di.LT.0 ) THEN
               di=di/EAMX(ia,2)
               Chilo = Chilo +
     &                 (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &                 /EAMX(ia,2))**2
            ELSE
               di=di/EAMX(ia,3)
               Chilo = Chilo +
     &                 (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &                 /EAMX(ia,3))**2
            ENDIF
            Chi = Chi + di*di
         ELSE
            di = (ELM(ib)-EAMX(ia,1))
            IF ( di.LT.0 ) THEN
              di=di/EAMX(ia,2)
              Chilo = Chilo +
     &                (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &                /EAMX(ia,2))**2
            ELSE
              di=di/EAMX(ia,3)
              Chilo = Chilo +
     &                (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &                /EAMX(ia,3))**2
            ENDIF
            Chi = Chi + di*di
         ENDIF
      ENDDO
      END
