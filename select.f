 
C----------------------------------------------------------------------
C SUBROUTINE SELECT
C
C Called by: GOSIA
C
C Purpose: integrate the functionality of the program SELECT into gosia as
C          OP,SELE
C
C PJN April 2008
C
      SUBROUTINE SELECT
      IMPLICIT NONE
      REAL*8 a , al , am , y
      INTEGER*4 i , ie , iexp , indx , ixf , j , l , lm , lu , lum , 
     &        lx , memax
      DIMENSION lm(1500) , y(175,1500) , a(1500,1500)

      ixf = 0
      REWIND 18
      READ (18,*) memax
      DO i = 1 , memax
         DO j = 1 , memax
           a(i,j) = 0.d0
         ENDDO
      ENDDO
      DO i = 1 , 1000
         READ (18,*) l , j
         IF ( l.EQ.0 ) GOTO 100
         IF ( j.NE.0 ) THEN
            a(l,j) = -1.d0
            a(j,l) = -1.d0
         ENDIF
      ENDDO
 100  ie = 1
 200  lum = 0
      DO i = 1 , 175
         DO j = 1 , memax
            y(i,j) = 0.
         ENDDO
      ENDDO
      DO i = 1 , 10000
         READ (18,*) lu , indx , iexp , al
         IF ( iexp.NE.ie ) GOTO 300
         lum = MAX0(lum,lu)
         y(lu,indx) = al
      ENDDO
 300  BACKSPACE 18
      ie = iexp
      IF ( ie.EQ.0 ) ixf = 1
      DO i = 1 , memax
         DO j = i , memax
           IF ( a(i,j).NE.-1.d0 ) THEN
               DO l = 1 , lum
                  a(i,j) = a(i,j) + y(l,i)*y(l,j)
               ENDDO
               a(j,i) = a(i,j)
            ENDIF
         ENDDO
      ENDDO
      IF ( ixf.NE.1 ) GOTO 200
      DO i = 1 , memax
         am = 0.
         DO j = 1 , memax
            am = MAX(a(i,j),am)
         ENDDO
         IF ( am.EQ.0.d0 ) am = 1.
         DO j = 1 , memax
           IF ( a(i,j).NE.-1.d0 ) THEN
               a(i,j) = a(i,j)/am
               IF ( a(i,j).LT..1d0 ) THEN
                 a(i,j) = 0.d0
                  GOTO 350
               ENDIF
            ENDIF
            a(i,j) = 1.d0
 350        CONTINUE
         ENDDO
         WRITE (10,*) (INT(a(i,j)),j=1,memax)
      ENDDO
      DO i = 1 , memax
         WRITE (22,99001) i
99001    FORMAT (10X,'ME=',1I3,3X,'PREDICTED CORRELATION'/)
         lx = 0
         DO j = 1 , memax
           IF ( a(i,j).NE.0.d0 ) THEN
               lx = lx + 1
               lm(lx) = j
            ENDIF
         ENDDO
         WRITE (22,*) (lm(j),j=1,lx)
      ENDDO
      END
