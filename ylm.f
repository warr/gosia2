
C----------------------------------------------------------------------
C SUBROUTINE YLM
C
C Called by: ANGULA
C
C Purpose: evaluate the even spherical harmonics.
C
C Uses global variables:
C      IEXP   - experiment number
C      IAXS   - axial symmetry flag
C
C Formal parameters:
C      Theta  - theta for which to evaluate (read only)
C      Ylmr   - return value for that theta (write only)
C
C Ylmr(l,m) = 1 / \sqrt{4 \pi} Y_{2l}^{m - 1}
C
C Note the factor of 1 / \sqrt{4 \pi} compared to the orthonormal spherical
C harmonics.
C
C 0.0889703179  = sqrt(5) / (8 pi)
C 0.0298415518  = 3 / (32 pi)
C 0.0179325408  = sqrt(13) / (64 pi)
C 0.1089659406  = sqrt(30) / (16 pi)
C -0.2179318812 = -1 * sqrt(30) / (8 pi)
C 0.1248361677  = 3 * sqrt(70) / (64 pi)
C -0.3530900028 = -3 * sqrt(140) / (32 pi)
C 0.0943672726  = 3 * sqrt(10) / (32 pi)
C -0.1334554768 = -3 * sqrt(20) / (32 pi)
C etc.

      SUBROUTINE YLM(Theta,Ylmr)
      IMPLICIT NONE
      REAL*8 ct , ctsq , st , Theta , Ylmr
      INTEGER*4 i , j , l , lf , m
      INCLUDE 'kin.inc'
      INCLUDE 'fconst.inc'
      DIMENSION Ylmr(9,9) , st(7)

      ct = COS(Theta)
      ctsq = ct*ct
      IF ( IAXS(IEXP).EQ.0 ) THEN ! If axially symmetric
         Ylmr(1,1) = SQRT(5.D0)/8.D0/pi*(3.D0*ctsq-1.D0)
         Ylmr(2,1) = 3.D0/32.D0/pi*((35.D0*ctsq-30.D0)*ctsq+3.D0)
         Ylmr(3,1) = SQRT(13.D0)/64.D0/pi*(((231.D0*ctsq-315.D0)*ctsq+
     &               105.D0)*ctsq-5.D0)
         GOTO 99999
      ENDIF
      st(1) = SIN(Theta)
      DO i = 2 , 7
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      Ylmr(1,3) = SQRT(30.D0)/16.D0/pi
      Ylmr(1,2) = -SQRT(30.D0)/8.D0/pi*ct
      Ylmr(1,1) = SQRT(5.D0)/8.D0/pi*(3.D0*ctsq-1.D0)
      Ylmr(2,5) = 3.D0*SQRT(70.D0)/64.D0/pi
      Ylmr(2,4) = -3.D0*SQRT(140.D0)/32.D0/pi*ct
      Ylmr(2,3) = 3.D0*SQRT(10.D0)/32.D0/pi*(7.D0*ctsq-1.D0)
      Ylmr(2,2) = -3.D0*SQRT(20.D0)/32.D0/pi*ct*(7.D0*ctsq-3.D0)
      Ylmr(2,1) = 3.0D0/32.D0/pi*((35.D0*ctsq-30.D0)*ctsq+3.D0)
      Ylmr(3,7) = SQRT(3003.D0)/128.D0/pi
      Ylmr(3,6) = -SQRT(9009.D0)/64.D0/pi*ct
      Ylmr(3,5) = 3.D0*SQRT(182.D0)/128.D0/pi*(11.D0*ctsq-1.D0)
      Ylmr(3,4) = -SQRT(1365.D0)/64.D0/pi*ct*(11.D0*ctsq-3.D0)
      Ylmr(3,3) = SQRT(1365.D0)/128.D0/pi*((33.D0*ctsq-18.D0)*ctsq+1.D0)
      Ylmr(3,2) = -SQRT(546.D0)/64.D0/pi*ct*((33.D0*ctsq-30.D0)*ctsq+
     &            5.D0)
      Ylmr(3,1) = SQRT(13.D0)/64.D0/pi*(((231.D0*ctsq-315.D0)*ctsq+
     &            105.D0)*ctsq-5.D0)
      DO l = 1 , 3
         lf = 2*l + 1
         DO m = 2 , lf
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)
         ENDDO
      ENDDO
99999 END
