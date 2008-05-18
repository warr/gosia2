 
C----------------------------------------------------------------------
C SUBROUTINE PODZIEL
C
C Called by: APRAM
C
C Purpose: subdivide matrix operators if the summation doesn't converge.
C
C Uses global variables:
C      IDIVE  - number of subdivisions
C      LP2    - maximum number of matrix elements (1500)
C      QAPR   - approximate Coulomb amplitudes
C
C We use the identity: exp(A) \bar{a} = exp(A/2) exp(A/2)\bar{a}.
C
C Formal parameters:
C      I      - flag (I=1,2,3) I=3 means initialise
C      J      - experiment number
 
      SUBROUTINE PODZIEL(I,J)
      IMPLICIT NONE
      INTEGER*4 I , IAPR , IDIVE , ISEX , J , k , l , l1 , l2 , LERF , 
     &          LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , 
     &          LP4 , LP6
      INTEGER*4 LP7 , LP8 , LP9
      REAL*8 QAPR
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(75)
      COMMON /APRX  / LERF , IDIVE(50,2)

      IF ( I.NE.3 ) THEN
         IF ( I.EQ.1 ) THEN
            l1 = IDIVE(J,1)
            IDIVE(J,1) = l1 + 1
            GOTO 100
         ELSE
            l1 = IDIVE(J,2)
            IDIVE(J,2) = l1 + 1
         ENDIF
      ENDIF

      l2 = IDIVE(J,2)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2 ! For each matrix element
         DO l = 1 , 7
            QAPR(k,2,l) = QAPR(k,2,l)*l1/l2
         ENDDO
      ENDDO

      IF ( I.EQ.2 ) WRITE (22,99001) J , IDIVE(J,1) , l2
      IF ( I.NE.3 ) RETURN
      
 100  l2 = IDIVE(J,1)
      IF ( I.EQ.3 ) l1 = 1
      DO k = 1 , LP2 ! For each matrix element
         DO l = 1 , 7
            QAPR(k,1,l) = QAPR(k,1,l)*l1/l2
         ENDDO
      ENDDO
       
      IF ( I.EQ.1 ) WRITE (22,99001) J , l2 , IDIVE(J,2)
      RETURN
      
99001 FORMAT (5X,'*****',1X,'EXP(A) EXPANSION FAILURE!',1X,'*****'/5X,
     &        'EXPERIMENT',1X,1I2,3X,'NEW SUBDIVISION',1X,'(',1I1,',',
     &        1I1,')')
      END
