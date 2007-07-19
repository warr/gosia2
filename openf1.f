 
C----------------------------------------------------------------------
C SUBROUTINE OPENF1
C
C Called by: GOSIA
C
C Purpose: open files to specified units.
C
      SUBROUTINE OPENF1
      IMPLICIT NONE
      INTEGER*4 i , j , k
      CHARACTER name*60 , opt1*20 , opt2*20

 100  READ * , i , j , k ! unit, old/new/unknown, formatted/unformatted
      IF ( i.EQ.0 ) RETURN
      IF ( j.EQ.1 ) opt1 = 'OLD'
      IF ( j.EQ.2 ) opt1 = 'NEW'
      IF ( j.EQ.3 ) opt1 = 'UNKNOWN'
      IF ( k.EQ.1 ) opt2 = 'FORMATTED'
      IF ( k.EQ.2 ) opt2 = 'UNFORMATTED'
      READ 99001 , name
99001 FORMAT (A)

      OPEN (i,IOSTAT=k,FILE=name,STATUS=opt1,FORM=opt2)
c      IF (K.EQ.0) WRITE(6,1030) 'OPENED ',NAME
c 1030 FORMAT (1X,2A)
c      WRITE(6,1010) ' IO-num = ',I,OPT1,OPT2
c 1010 FORMAT (1X,A,I4,2(1x,A))
      IF ( k.EQ.0 ) GOTO 100
c      WRITE (6,1020) 'PROBLEMS OPENING ',NAME,K
c 1020 FORMAT(A,A,I6)
      END
