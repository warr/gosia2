 
C----------------------------------------------------------------------
C FUNCTION MEM
C
C Called by: SEQ, NEWLV
C
C Purpose: calculates an index to a matrix element given two level indices
C      and the multipolarity.
C
C Uses global variables:
C      LDNUM  - number of matrix elements with each multipolarity populating level
C      LEAD   - pair of levels involved in each matrix element
C      MULTI  - number of matrix elements having a given multipolarity
C
C Formal parameters:
C      N1     - level number for first level
C      N2     - level number for second level
C      N3     - multipolarity
C
C Return value:
C      Index of matrix element
 
      INTEGER*4 FUNCTION MEM(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMDA , LAMMAX , LDNUM , LEAD , msum , MULTI , N1 , 
     &          n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)

      msum = 0
      IF ( N3.NE.1 ) THEN
         n3m = N3 - 1
         DO k = 1 , n3m ! For each multipolarity up to one below the one we want
            msum = msum + MULTI(k) ! Add the number of matrix elements for that multipolarity
         ENDDO
      ENDIF

C     msum is now an index to the start of the matrix elements for the chosen multipolarity

      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m ! For each level up to one below the one we want
            msum = msum + LDNUM(N3,k) ! Add the number of matrix elements for that level and multipolarity
         ENDDO
      ENDIF

C     msum is now an index to the start of the matrix elements for the appropriate multipolarity and level

      n1m = msum + 1
      n3m = n1m + LDNUM(N3,N1)
      DO k = n1m , n3m ! Loop over matrix elements associated with that level and multipolarity
         msum = msum + 1
         IF ( LEAD(2,k).EQ.N2 ) GOTO 100 ! If it is the right one goto 100
      ENDDO

 100  MEM = msum ! MEM is now the index to the matrix element we want
      END
