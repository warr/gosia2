C----------------------------------------------------------------------
C SUBROUTINE SPLNER
C
C Called by: CCLKUP, EFFIX , GOSIA
C Calls:     FUNC, FUNC1, SPLINE, SPLINT
C
C Purpose: interpolates using a cubic spline
C
C Formal parameters:
C      X      - x-coordinate of input data
C      Yr     - y-coordinate of input data
C      N      - number of data points
C      Xx     - value for which to interpolate
C      Yy     - result of interpolation
C      Iscal  - mode: 1 = linear, 2 = logarithmic, 3 = square root
C
C Note that the effect of FUNC and FUNC1 depends on Iscal:
C Iscal = 1   FUNC(y) = y        FUNC1(y) = y
C Iscal = 2   FUNC(y) = ln(y)    FUNC1(y) = exp(y)
C Iscal = 3   FUNC(y) = sqrt(y)  FUNC1(y) = y^2
C
      SUBROUTINE SPLNER(X,Yr,N,Xx,Yy,Iscal)
      IMPLICIT NONE
      REAL*8 FUNC , FUNC1
      INTEGER*4 N
      REAL*8 X(*) , Yr(*) , w(1500)
      REAL*8 yp1 , ypn , y(1500) , ys
      INTEGER*4 i , Iscal
      REAL*8 Xx , Yy

C     We need at least three points, so if we don't have them, use the
C     Lagrangian method instead
      IF ( N.LE. 3 ) THEN
        CALL LAGRAN(X,Yr,N,1,Xx,Yy,Iscal,1)
        RETURN
      ENDIF 

C     Apply the scaling function
      DO i = 1 , N
        y(i) = FUNC(Yr(i),Iscal)
      ENDDO

C     Get the slope at each end
      yp1 = (y(2)-y(1))/(X(2)-X(1))
      ypn = (y(N)-y(N-1))/(X(N)-X(N-1))

C     Fit a spline
      CALL SPLINE(X,y,N,yp1,ypn,w)

C     Evaluate the spline at the desired point
      CALL SPLINT(X,y,w,N,Xx,ys)

C     Apply the inverse scaling function
      Yy = FUNC1(ys,Iscal)
      RETURN
      END
 
