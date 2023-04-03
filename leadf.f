
C----------------------------------------------------------------------
C FUNCTION LEADF
C
C Called by: LAIAMP, LSLOOP, NEWLV, SEQ
C
C Uses global variables:
C      LDNUM  - number of matrix elements with each multipolarity populating levels
C      LEAD   - pair of levels involved in each matrix element
C      MULTI  - number of matrix elements with a given multipolarity
C
C Formal parameters:
C      N1     - index of initial level
C      N2     - index of matrix element for given level and multipolarity
C      N3     - multipolarity
C
C Purpose: calculate the level number for the final level associated with the
C      matrix element index N2, initial level index N1 and multipolarity N3.

      INTEGER*4 FUNCTION LEADF(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMDA , LAMMAX , LDNUM , LEAD , lsum , MULTI , N1 ,
     &          n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,75) , LAMMAX ,
     &                MULTI(8)

      lsum = 0
      n3m = N3 - 1
      IF ( n3m.NE.0 ) THEN
         DO k = 1 , n3m ! Loop over multipolarities lower than one required
            lsum = lsum + MULTI(k)
         ENDDO
      ENDIF

C     lsum now points to start of the multipolarity N3

      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m ! Loop over levels below the selected one
            lsum = lsum + LDNUM(N3,k)
         ENDDO
      ENDIF

C     lsum now points to start of level N1 for multipolarity N3

      n1m = lsum + N2

C     n1m now points to the appropriate matrix element N2 for initial level N1
C     and multipolarity N3

      LEADF = LEAD(2,n1m) ! Get the final level for this matrix element
      END
