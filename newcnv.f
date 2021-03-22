C----------------------------------------------------------------------
C FUNCTION NEWCNV
C
C Called by: CONV
C
C Purpose: find the conversion coefficient calculated by OP,BRIC and stored
C          in a file
C
C Formal parameters:
C      Ega    - gamma energy
C      N      - multipolarity N=1,2,3 = E1,2,3 and N=4,5 = M1,2 (not as elsewhere!)
C
C Uses global variables:
C      JZB    - unit to read from
C Return value:
C      conversion coefficient interpolated to energy Ega

      REAL*8 FUNCTION NEWCNV(Ega,N)
      IMPLICIT NONE

      INTEGER*4 isfirst(2), i, j, N, nenergies(2)
      REAL*8 energies(2, 1500), bricc(2, 1500, 5), Ega
      SAVE energies, bricc, isfirst, nenergies
      DATA isfirst/1,1/
      INCLUDE 'switch.inc'

C     The first time, we need to read the data
      IF ( isfirst(JZB-24).eq.1 ) THEN
        rewind(JZB+3)
        isfirst(JZB-24)= 0
        DO i = 1, 1500
          nenergies(JZB-24) = i - 1
          READ(JZB+3,*,END=100) energies(JZB-24,i),
     &      (bricc(JZB-24,i,j),j=1,5)
        ENDDO
      ENDIF

C     Check multipolarity is valid
 100  IF ( N.LT.1.OR.N.GT.5 ) THEN
         NEWCNV = 0.0
         RETURN
      ENDIF

C     Search for the energy in the list

      DO i = 1, nenergies(JZB-24)
        IF (ABS(Ega - energies(JZB-24,i)) .LT. 1E-3) THEN
           NEWCNV = bricc(JZB-24,i,N)
           return
        ENDIF
      ENDDO

C     We get here if the energy isn't in the list, so stop with an error
C     message
      WRITE (*,'(A,F7.3,A)')
     & 'Unable to find conversion coefficients for ',
     &  Ega, ' MeV'
      STOP 'Missing conversion coefficients'

      END
