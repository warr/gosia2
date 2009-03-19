 
C----------------------------------------------------------------------
C FUNCTION CONV
C
C Called by: BRANR, PTICC, SEQ
C Calls:     LAGRAN, NEWCNV, SPLNER
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
C      N      - multipolarity N=1,2,3 = E1,2,3 and N=4,5 = M1,2 (not as elsewhere!)
C
C Return value:
C      conversion coefficient interpolated to energy Ega

      REAL*8 FUNCTION CONV(Ega,N)
      IMPLICIT NONE
      REAL*8 cpo , cpo1 , cv , Ega , NEWCNV
      INTEGER*4 j , N , n1 , nen
      DIMENSION cpo(101) , cpo1(101)
      INCLUDE 'ccc.inc'

C     If the number of conversion coefficients entered by the user is negative
C     then use read the conversion coefficients from a file on unit 29.
      IF ( NICC.LE.0 ) THEN
         CONV=NEWCNV(Ega,N)
         RETURN
      ENDIF

      IF ( N.EQ.0 ) THEN ! If no multipolarity defined
         CONV = 0.0
      ELSEIF ( ABS(CC(1,N)).LT.1.E-9 ) THEN ! If no conversion coefficients given for this multipolarity
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
         IF ( ISPL.EQ. 0 ) CALL LAGRAN(cpo1,cpo,4,1,Ega,cv,2,1)
         IF ( ISPL.EQ. 1 ) CALL SPLNER(cpo1,cpo,4,Ega,cv,2)
         CONV = cv
         RETURN
      ENDIF
      END
