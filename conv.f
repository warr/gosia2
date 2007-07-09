 
C----------------------------------------------------------------------
 
      REAL*8 FUNCTION CONV(Ega,N)
      IMPLICIT NONE
      REAL*8 AGEli , CC , cpo , cpo1 , cv , EG , Ega , Q
      INTEGER*4 j , N , n1 , NANg , nen , NICc
      DIMENSION cpo(51) , cpo1(51)
      COMMON /CCC   / EG(50) , CC(50,5) , AGEli(50,200,2) , Q(3,200,8) , 
     &                NICc , NANg(200)
      IF ( N.EQ.0 ) THEN
         CONV = 0.0
      ELSEIF ( ABS(CC(1,N)).LT.1.E-9 ) THEN
         CONV = 0.0
      ELSE
         nen = 4
         DO j = 1 , NICc
            IF ( Ega.LE.EG(j) ) GOTO 50
         ENDDO
 50      n1 = j - 2
         IF ( n1.LT.1 ) n1 = 1
         IF ( (j+1).GT.NICc ) n1 = n1 - 1
         IF ( NICc.LE.4 ) THEN
            n1 = 1
            nen = NICc
         ENDIF
         DO j = 1 , nen
            cpo(j) = CC(n1+j-1,N)
            cpo1(j) = EG(n1+j-1)
         ENDDO
         CALL LAGRAN(cpo1,cpo,4,1,Ega,cv,2,1)
         CONV = cv
         RETURN
      ENDIF
      END
