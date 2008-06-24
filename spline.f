C----------------------------------------------------------------------
C SUBROUTINE SPLINE
C
C Called by: SPLNER
C
C Purpose: fit a spline to data points
C
C Formal parameters:
C      X      - x-coordinate of input data
C      Y      - y-coordinate of input data
C      N      - number of data points
C      Yp1    - slope of first two data points
C      Yyn    - slope of last two data points
C      Y2     - second derivative vector
C
      SUBROUTINE SPLINE(X,Y,N,Yp1,Ypn,Y2)
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 Yp1 , Ypn , X(*) , Y(*) , Y2(*)
      INTEGER*4 i , k
      REAL*8 p , qn , sig , un , u(1500)
 
      IF ( Yp1.GT..99E30 ) THEN
         Y2(1) = 0.
         u(1) = 0.
      ELSE
         Y2(1) = -.5
         u(1) = (3./(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-Yp1)
      ENDIF
      DO i = 2 , N - 1
         sig = (X(i)-X(i-1))/(X(i+1)-X(i-1))
         p = sig*Y2(i-1) + 2.
         Y2(i) = (sig-1.)/p
         u(i) = (6.*((Y(i+1)-Y(i))/(X(i+1)-X(i))-(Y(i)-Y(i-1))/(X(i)-X(i
     &          -1)))/(X(i+1)-X(i-1))-sig*u(i-1))/p
      ENDDO
      IF ( Ypn.GT..99E30 ) THEN
         qn = 0.
         un = 0.
      ELSE
         qn = .5
         un = (3./(X(N)-X(N-1)))*(Ypn-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N) = (un-qn*u(N-1))/(qn*Y2(N-1)+1.)
      DO k = N - 1 , 1 , -1
         Y2(k) = Y2(k)*Y2(k+1) + u(k)
      ENDDO
      END
