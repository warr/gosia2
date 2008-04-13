      FUNCTION MEM(N1,N2,N3)
      IMPLICIT NONE
      INTEGER*4 k , LAMDA , LAMMAX , LDNUM , LEAD , MEM , msum , MULTI , 
     &          N1 , n1m , N2 , N3 , n3m
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      msum = 0
      IF ( N3.NE.1 ) THEN
         n3m = N3 - 1
         DO k = 1 , n3m
            msum = msum + MULTI(k)
         ENDDO
      ENDIF
      n1m = N1 - 1
      IF ( n1m.NE.0 ) THEN
         DO k = 1 , n1m
            msum = msum + LDNUM(N3,k)
         ENDDO
      ENDIF
      n1m = msum + 1
      n3m = n1m + LDNUM(N3,N1)
      DO k = n1m , n3m
         msum = msum + 1
         IF ( LEAD(2,k).EQ.N2 ) GOTO 100
      ENDDO
 100  MEM = msum
      END
