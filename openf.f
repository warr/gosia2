 
C----------------------------------------------------------------------
C SUBROUTINE OPENF
C
C Called by: GOSIA
C
C Purpose: open files to specified units.
C
C Uses global variables:
C      JZB    - unit to read from
C
C The function reads three integers, the first of which is the unit to use for
C the open statement. The second is 1 if the file is required to exist already,
C 2 if it is required not to exist and 3 if it does not matter. The third is 1
C if the file is formatted and 2 if it is unformatted. A second line is read,
C which gives the name of the file to associate with that unit. If the unit is
C zero, the function returns. It keeps looping until a unit zero is reached.

      SUBROUTINE OPENF
      IMPLICIT NONE
      INTEGER*4 i , IBPS , IUNIT3 , j , JZB , k
      CHARACTER name*60 , opt1*20 , opt2*20
      COMMON /SWITCH/ JZB , IBPS , IUNIT3

 100  READ (JZB,*) i , j , k ! unit, old/new/unknown, formatted/unformatted
      IF ( i.EQ.0 ) RETURN
      IF ( j.EQ.1 ) opt1 = 'OLD'
      IF ( j.EQ.2 ) opt1 = 'NEW'
      IF ( j.EQ.3 ) opt1 = 'UNKNOWN'
      IF ( k.EQ.1 ) opt2 = 'FORMATTED'
      IF ( k.EQ.2 ) opt2 = 'UNFORMATTED'
      READ (JZB,99001) name ! name of file
99001 FORMAT (A)

C     Now open the file
      OPEN (i,IOSTAT=k,FILE=name,STATUS=opt1,FORM=opt2)
c      IF ( k.EQ.0 ) WRITE (6,99002) 'OPENED ' , name
c99002 FORMAT (1X,2A)
c      WRITE (6,99003) ' IO-num = ' , i , opt1 , opt2
c99003 FORMAT (1X,A,I4,2(1x,A))
      IF ( k.EQ.0 ) GOTO 100
c      WRITE (6,99004) 'PROBLEMS OPENING ' , name , k
c99004 FORMAT (A,A,I6)
      END
