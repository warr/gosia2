      SUBROUTINE LAGRAN(X,Y,Ndata,Ipc,Xx,Yy,Iscal,Irc)
      IMPLICIT NONE
      REAL*8 arh , FUNC , FUNC1 , t , w , X , Xx , Y , y1 , Yy
      INTEGER*4 i , Ipc , Irc , Iscal , j , Ndata
      DIMENSION X(51) , Y(51) , w(51) , arh(51,51)
      IF ( Irc.EQ.2 ) THEN
      ELSEIF ( Irc.EQ.3 ) THEN
         DO i = 1 , Ndata
            t = 1.
            DO j = 1 , Ndata
               IF ( i.NE.j ) t = t*(Xx-X(j))/(X(i)-X(j))
            ENDDO
            arh(Ipc,i) = t
         ENDDO
         GOTO 100
      ELSEIF ( Irc.EQ.4 ) THEN
         GOTO 100
      ELSE
         DO i = 1 , Ndata
            w(i) = 1.
            DO j = 1 , Ndata
               IF ( i.NE.j ) w(i) = w(i)*(Xx-X(j))/(X(i)-X(j))
            ENDDO
         ENDDO
      ENDIF
      Yy = 0.
      DO j = 1 , Ndata
         y1 = Y(j)
         Yy = Yy + w(j)*FUNC(y1,Iscal)
      ENDDO
      Yy = FUNC1(Yy,Iscal)
      RETURN
 100  Yy = 0.
      DO j = 1 , Ndata
         y1 = Y(j)
         Yy = Yy + arh(Ipc,j)*FUNC(y1,Iscal)
      ENDDO
      Yy = FUNC1(Yy,Iscal)
      END
