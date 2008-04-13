      FUNCTION SIMIN(Np,H,Y)
      IMPLICIT NONE
      REAL*8 ee , H , SIMIN , sm , Y
      INTEGER*4 ik , in , Np
      DIMENSION Y(101)
      IF ( Np.GE.3 ) THEN
         ik = Np - 2
         sm = Y(1) + Y(Np)
         DO in = 1 , ik
            ee = in/2.
            sm = sm + 2.*Y(in+1)/(1.+INT(ee)-ee)
         ENDDO
         SIMIN = sm*H/3.
         RETURN
      ELSEIF ( Np.EQ.1 ) THEN
         SIMIN = Y(1)
         GOTO 99999
      ENDIF
      SIMIN = (Y(1)+Y(2))*H/2.
99999 END
