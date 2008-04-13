      SUBROUTINE PRIM(N)
      IMPLICIT NONE
      INTEGER*4 i , IP , IPI , KF , N , nni , nnk
      REAL*8 PILOG
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)
      nnk = N
      DO i = 1 , 26
         nni = nnk
         IPI(i) = 0
 50      nni = nni/IP(i)
         IF ( IP(i)*nni.EQ.nnk ) THEN
            IPI(i) = IPI(i) + 1
            nnk = nni
            GOTO 50
         ENDIF
      ENDDO
      END
