      SUBROUTINE FAKP
      IMPLICIT NONE
      INTEGER*4 i , IP , IPI , k , KF , l
      REAL*8 PILOG , x
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      DO i = 1 , 26
         x = DBLE(IP(i))
         PILOG(i) = LOG(x)
      ENDDO
      DO l = 1 , 26
         KF(1,l) = 0
         KF(2,l) = 0
      ENDDO
      DO k = 3 , 101
         CALL PRIM(k-1)
         DO i = 1 , 26
            KF(k,i) = KF(k-1,i) + IPI(i)
         ENDDO
      ENDDO
      END
