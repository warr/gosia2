C----------------------------------------------------------------------
C SUBROUTINE SPLNER
C
C Called by: CCLKUP
C Calls:     FUNC, FUNC1, SPLINE, SPLINT
C
C Purpose: interpolates using a cubic spline
C
C Formal parameters:
C      X      - x-coordinate of input data
C      Yr     - y-coordinate of input data
C      Ndata  - number of data points
C      Ipc    - index for storing results
C      Xx     - value for which to interpolate
C      Yy     - result of interpolation
C      Iscal  - mode: 1 = linear, 2 = exponential, 3 = square root
C      Irc    - weighting mode
C
C Note that the effect of FUNC and FUNC1 depends on Iscal:
C Iscal = 1   FUNC(y) = y        FUNC1(y) = y
C Iscal = 2   FUNC(y) = ln(y)    FUNC1(y) = exp(y)
C Iscal = 3   FUNC(y) = sqrt(y)  FUNC1(y) = y^2
C
      SUBROUTINE SPLNER(X,Yr,N,Ipc,Xx,Yy,Iscal,Irc)
      IMPLICIT NONE
      REAL*8 FUNC , FUNC1
      INTEGER N
      REAL*8 X(*) , Yr(*) , w(1500) , arh(1500,1500)
      REAL*8 yp1 , ypn , y(1500) , ys
      INTEGER*4 i , Iscal , Irc , Ipc
      REAL*8 Xx , Yy

      IF ( Irc.EQ.2 ) THEN
      ELSEIF ( Irc.EQ.3 ) THEN

         DO i = 1 , N
            y(i) = FUNC(Yr(i),Iscal)
c      print*,iscal,x(i),yr(i),y(i),irc,ipc
         ENDDO

         yp1 = (y(2)-y(1))/(X(2)-X(1))
         ypn = (y(N)-y(N-1))/(X(N)-X(N-1))
         CALL SPLINE(X,y,N,yp1,ypn,w)
         DO i = 1 , N
            arh(Ipc,i) = w(i)
         ENDDO
         CALL SPLINT(X,y,w,N,Xx,ys)
         Yy = FUNC1(ys,Iscal)
c         PRINT * , 'Spline' , Xx , ys , Yy , Iscal , Irc , Ipc
         RETURN
      ELSEIF ( Irc.EQ.4 ) THEN

         DO i = 1 , N
            w(i) = arh(Ipc,i)
         ENDDO

         CALL SPLINT(X,y,w,N,Xx,ys)
         Yy = FUNC1(ys,Iscal)
         PRINT * , 'Spline' , Xx , ys , Yy , Iscal , Irc , Ipc
         GOTO 99999
      ELSE


         DO i = 1 , N
            y(i) = FUNC(Yr(i),Iscal)
c      print*,iscal,x(i),yr(i),y(i),irc
         ENDDO

         yp1 = (y(2)-y(1))/(X(2)-X(1))
         ypn = (y(N)-y(N-1))/(X(N)-X(N-1))
         CALL SPLINE(X,y,N,yp1,ypn,w)
      ENDIF
      CALL SPLINT(X,y,w,N,Xx,ys)
      Yy = FUNC1(ys,Iscal)
c      PRINT * , 'Spline' , Xx , ys , Yy , Iscal , Irc , Ipc
      RETURN
99999 END

