 
C----------------------------------------------------------------------
C FUNCTION CONV
C
C Called by: BRANR, PTICC, SEQ
C Calls:     LAGRAN
C
C Purpose: calculate the conversion coefficient at a particular energy by
C interpolating over the values provided by the user.
C
C Uses global variables:
C      CC     - conversion coefficients
C      EG     - energies for conversion coefficients
C      NICC   - number of conversion coefficients
C
C Formal parameters:
C      Ega    - gamma energy
C      N      -
C
C Return value:
C      conversion coefficient interpolated to energy Ega

      REAL*8 FUNCTION CONV(Ega,N)
      IMPLICIT NONE
      REAL*8 AGELI , CC , cpo , cpo1 , cv , EG , Ega , Q
      INTEGER*4 j , N , n1 , NANG , nen , NICC
      DIMENSION cpo(51) , cpo1(51)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)

      IF ( N.EQ.0 ) THEN
         CONV = 0.0
      ELSEIF ( ABS(CC(1,N)).LT.1.E-9 ) THEN
         CONV = 0.0
      ELSE
         nen = 4
         DO j = 1 , NICC ! Loop over coefficients provided by user
            IF ( Ega.LE.EG(j) ) GOTO 50
         ENDDO
 50      n1 = j - 2
         IF ( n1.LT.1 ) n1 = 1
         IF ( (j+1).GT.NICC ) n1 = n1 - 1
         IF ( NICC.LE.4 ) THEN
            n1 = 1
            nen = NICC
         ENDIF
         DO j = 1 , nen
            cpo(j) = CC(n1+j-1,N)
            cpo1(j) = EG(n1+j-1)
         ENDDO
C        Interpolate 
         CALL LAGRAN(cpo1,cpo,4,1,Ega,cv,2,1)
         CONV = cv
         RETURN
      ENDIF
      END
