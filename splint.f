C----------------------------------------------------------------------
C SUBROUTINE SPLINT
C
C Called by: SPLNER
C
C Purpose: evaluate spline
C
C Formal parameters:
C      Xa     - x-coordinate of input data
C      Ya     - y-coordinate of input data
C      Ya2    - second derivative vector
C      N      - number of data points
C      Xx     - value for which to evaluate spline
C      Yy     - result
C
      SUBROUTINE SPLINT(Xa,Ya,Y2a,N,Xx,Yy)
c
c xa,ya - tabulated function
c n - number of tabulated points
c y2a - second derivatives vector
c xx - interpolated point
c yy - intetpoalted value
c
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 Xx , Yy , Xa(*) , Y2a(*) , Ya(*)
      INTEGER*4 k , khi , klo
      REAL*8 a , b , h
      DATA k/0/

      klo = 1
      khi = N
 100  IF ( khi-klo.GT.1 ) THEN
         k = (khi+klo)/2
         IF ( Xa(k).GT.Xx ) THEN
            khi = k
         ELSE
            klo = k
         ENDIF
         GOTO 100
      ENDIF
      h = Xa(khi) - Xa(klo)
C
C  Linear extrapolation
C
      IF ( ABS(h).LT.1E-9 ) THEN
         IF ( Xx.LT.Xa(1) ) h = 1
         IF ( Xx.GT.Xa(N) ) k = N - 1
         a = (Ya(k)-Ya(k+1))/(Xa(k)-Xa(k+1))
         b = (Ya(k)+Ya(k+1)-a*(Xa(k)+Xa(k+1)))*.5
         Yy = a*Xx + b
c         PRINT * , 'splint' , Xx , Yy , 'extrapolation'
         RETURN
      ENDIF
      a = (Xa(khi)-Xx)/h
      b = (Xx-Xa(klo))/h
      Yy = a*Ya(klo) + b*Ya(khi) + ((a**3-a)*Y2a(klo)+(b**3-b)*Y2a(khi))
     &     *(h**2)/6.
      END
