 
C----------------------------------------------------------------------
C SUBROUTINE GCF
C
C Called by: GAMATT
C
C Purpose:
C
C Formal parameters:
C      Tau    -
C      Thing  -
C      Q      -
 
      SUBROUTINE GCF(Tau,Thing,Q)
      IMPLICIT NONE
      REAL*8 A , b , D , dl , ev , ex , f , fint , od , ODL , Q , R , 
     &       Tau , Thing , XL , xm , yl , yu
      INTEGER*4 i , j , k , m
      COMMON /DIMX  / A , R , XL , D , ODL(200)
      DIMENSION f(101) , b(4) , Q(9)

      b(1) = ATAN2(A,D+XL)
      b(2) = ATAN2(A,D)
      b(3) = ATAN2(R,D+XL)
      b(4) = ATAN2(R,D)
      DO k = 1 , 9
         Q(k) = 0.0
         DO j = 1 , 3
            yl = b(j)
            yu = b(j+1)
            dl = (yu-yl)/100.
            DO m = 1 , 101
               xm = yl + dl*(m-1)
               IF ( j.EQ.2 ) THEN
                  ex = -Tau*XL/COS(xm)
               ELSEIF ( j.EQ.3 ) THEN
                  ex = Tau*(D*TAN(xm)-R)/SIN(xm)
               ELSE
                  ex = Tau*(A-(D+XL)*TAN(xm))/SIN(xm)
               ENDIF
               f(m) = SIN(xm)*(1-EXP(ex))*EXP(Thing/COS(xm))
               IF ( j.EQ.1 ) f(m) = f(m)*EXP(-Tau*(A/SIN(xm)-D/COS(xm)))
               IF ( k.EQ.1 ) THEN
               ELSEIF ( k.EQ.3 ) THEN
                  f(m) = f(m)*(1.5*COS(xm)**2-0.5)
               ELSEIF ( k.EQ.4 ) THEN
                  f(m) = f(m)*(2.5*COS(xm)**3-1.5*COS(xm))
               ELSEIF ( k.EQ.5 ) THEN
                  f(m) = f(m)*(4.375*COS(xm)**4-3.75*COS(xm)**2+.375)
               ELSEIF ( k.EQ.6 ) THEN
                  f(m) = f(m)*((63.*COS(xm)**5-70.*COS(xm)**3+15.)/8.)
               ELSEIF ( k.EQ.7 ) THEN
                  f(m) = f(m)
     &                   *((21.*COS(xm)**2*(11.*COS(xm)**4-15.*COS(xm)
     &                   **2+5.)-5.)/16.)
               ELSEIF ( k.EQ.8 ) THEN
                  f(m) = f(m)
     &                   *(429.*COS(xm)**7-693.*COS(xm)**5+315.*COS(xm)
     &                   **3-35.*COS(xm))/16.
               ELSEIF ( k.EQ.9 ) THEN
                  f(m) = f(m)
     &                   *(6435.*COS(xm)**8-12012.*COS(xm)**6+6930.*COS
     &                   (xm)**4-1260.*COS(xm)**2+35.)/128.
               ELSE
                  f(m) = f(m)*COS(xm)
               ENDIF
            ENDDO
            ev = 0.0
            od = 0.0
            DO m = 2 , 98 , 2
               ev = ev + f(m)
               od = od + f(m+1)
            ENDDO
            fint = dl/3.*(f(1)+4.*(ev+f(100))+2.*od+f(101))
            Q(k) = Q(k) + fint
         ENDDO
      ENDDO
      DO i = 1 , 8
         Q(i+1) = Q(i+1)/Q(1)
      ENDDO
      Q(1) = Q(1)/2.
      END