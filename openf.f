      SUBROUTINE OPENF
      IMPLICIT NONE
      INTEGER*4 i , IBPS , j , JZB , k
      CHARACTER name*60 , opt1*20 , opt2*20
      COMMON /SWITCH/ JZB , IBPS
 100  READ (JZB,*) i , j , k
      IF ( i.EQ.0 ) RETURN
      IF ( j.EQ.1 ) opt1 = 'OLD'
      IF ( j.EQ.2 ) opt1 = 'NEW'
      IF ( j.EQ.3 ) opt1 = 'UNKNOWN'
      IF ( k.EQ.1 ) opt2 = 'FORMATTED'
      IF ( k.EQ.2 ) opt2 = 'UNFORMATTED'
      READ (JZB,99001) name
99001 FORMAT (A)
      IF ( i.NE.25 .AND. i.NE.26 )
     &     OPEN (i,IOSTAT=k,FILE=name,STATUS=opt1,FORM=opt2)
      IF ( i.EQ.25 .OR. i.EQ.26 ) k = 0
c      IF (K.EQ.0) WRITE(6,1030) 'OPENED ',NAME
c 1030 FORMAT (1X,2A)
c      WRITE(6,1010) ' IO-num = ',I,OPT1,OPT2
c 1010 FORMAT (1X,A,I4,2(1x,A))
      IF ( k.EQ.0 ) GOTO 100
c      WRITE (6,1020) 'PROBLEMS OPENING ',NAME,K
c 1020 FORMAT(A,A,I6)
      END
