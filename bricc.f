C----------------------------------------------------------------------
C SUBROUTINE BRICC
C
C Called by: GOSIA
C Calls:     CCLKUP
C
C Purpose: evaluate internal conversion coefficients using the BrIcc
C          database for each transition energy that gosia needs. The
C          results are stored in the file on unit 29, which is read
C          the first time CONV is called.
C
C Uses global variables:
C      EN     - energy of level
C      IZ     - Z of investigated nucleus
C      LEAD   - pair of levels involved in each matrix element
C      MEMAX  - number of matrix elements
C
      SUBROUTINE BRICC
      IMPLICIT NONE
      REAL*8 temp , egamma
      INTEGER*4 i , j , ngamma
      DIMENSION egamma(1500)
      INCLUDE 'coex.inc'
      INCLUDE 'cx.inc'
      INCLUDE 'cexc.inc'
      INCLUDE 'clcom.inc'

      INTEGER*4 n1, n2
      CHARACTER*1024 idx_name, icc_name
      REAL*8 mycc(5), CCLKUP

C     Read the names of the files
      read (*,'(A)') idx_name
      read (*,'(A)') icc_name

C     Write to output
      write(22,'(/,3A)') 'OP,BRIC interpolation of conversion ',
     &  'coefficients from the BrIcc database, which will be used by',
     &  ' gosia.'
      write(22,'(2A)') 'Please cite T. Kibedi et al. NIM A589 (2008)',
     &  ' 202-229 for the conversion coefficients!'
      write(22,'(2A)') 'Energy [MeV]   E1           E2           E3',
     &  '           M1           M2'

C     Make sure we are at start of file that we want to write
      rewind(29)
      
C     Open the BrIcc database files
      OPEN (UNIT=30,FILE=idx_name,ACCESS='direct',RECL=2048,ERR=999)
      OPEN (UNIT=31,FILE=icc_name, ACCESS='direct',RECL=44,ERR=999,
     &      FORM='UNFORMATTED')

      ngamma = 0
      DO i = 1 , MEMAX ! For each matrix element
         n1 = LEAD(2,i) ! Upper level
         n2 = LEAD(1,i) ! Lower level
         IF ( n1.EQ.n2 ) GOTO 100 ! Ignore diagonal matrix elements

         temp = EN(n1) - EN(n2) ! Energy of transition

C        Now look to see if we have it already
         DO j = 1, ngamma
            IF ( ABS(temp - egamma(j)).LT.1E-6 ) GOTO 100
         ENDDO

C        We get here if we don't have it, so add it to the list
         ngamma = ngamma + 1
         egamma(ngamma) = temp
         mycc(1) = CCLKUP(IZ, temp * 1E3, 1)
         mycc(2) = CCLKUP(IZ, temp * 1E3, 2)
         mycc(3) = CCLKUP(IZ, temp * 1E3, 3)
         mycc(4) = CCLKUP(IZ, temp * 1E3, 6)
         mycc(5) = CCLKUP(IZ, temp * 1E3, 7)
         WRITE(22,'(F7.4,3X,1P,5E13.3)') temp, (mycc(j),j=1,5)
         WRITE(29,'(F7.4,3X,1P,5E13.3)') temp, (mycc(j),j=1,5)
 100  ENDDO

C     Close BrIcc database files
      CLOSE (30)
      CLOSE (31)
      RETURN

 999  STOP 'Unable to open BrIcc database files'
      END
