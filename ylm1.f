 
C----------------------------------------------------------------------
C SUBROUTINE YLM1
C
C Called by: ANGULA
C
C Purpose: evaluate the odd spherical harmonics.
C
C Formal parameters:
C      Theta  - theta for which to evaluate
C      Ylmr   - return value for that theta
C
C Ylmr(l,m) = 1 / \sqrt{4 \pi} Y_{l - 1}^{m - 1}
C
C Note the factor of 1 / \sqrt{4 \pi} compared to the orthonormal spherical
C harmonics.
C
C Note also that YLM1 and YLM have some values in common.
C e.g. YLM1(5,3) = YLM(2,3)
 
      SUBROUTINE YLM1(Theta,Ylmr)
      IMPLICIT NONE
      REAL*8 ct , ctsq , st , Theta , Ylmr
      INTEGER*4 i , j , l , m
      DIMENSION Ylmr(9,9) , st(9)
      INCLUDE 'fconst.inc'
      
      ct = COS(Theta)
      ctsq = ct*ct
      st(1) = SIN(Theta)
      DO i = 2 , 9
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      DO l = 2 , 9
         DO m = 1 , 9
            Ylmr(l,m) = 0.0
         ENDDO
      ENDDO
      Ylmr(2,2) = -SQRT(6.D0)/2.D0
      Ylmr(2,1) = SQRT(3.D0)*ct
      Ylmr(3,3) = SQRT(30.D0)/4.D0
      Ylmr(3,2) = -(SQRT(30.D0)/2.D0)*ct
      Ylmr(3,1) = (SQRT(5.D0)/2.D0)*(3.D0*ctsq-1.D0)
      Ylmr(4,4) = -SQRT(35.D0)/4.D0
      Ylmr(4,3) = (SQRT(210.D0)/4.D0)*ct
      Ylmr(4,2) = -(SQRT(21.D0)/4.D0)*(5.D0*ctsq-1.D0)
      Ylmr(4,1) = (SQRT(7.D0)/2.D0)*ct*(5.D0*ctsq-3.D0)
      Ylmr(5,5) = 3.D0*SQRT(70.D0)/16.D0
      Ylmr(5,4) = -(3.D0*SQRT(35.D0)/4.D0)*ct
      Ylmr(5,3) = (3.D0*SQRT(10.D0)/8.D0)*(7.D0*ctsq-1.D0)
      Ylmr(5,2) = -(3.D0*SQRT(5.D0)/4.D0)*ct*(7.D0*ctsq-3.D0)
      Ylmr(5,1) = (3.D0/8.D0)*((35.D0*ctsq-30.D0)*ctsq+3.D0)
      Ylmr(6,6) = -3.D0*SQRT(77.D0)/16.D0
      Ylmr(6,5) = (3.D0*SQRT(770.D0)/16.D0)*ct
      Ylmr(6,4) = -(SQRT(385.D0)/16.D0)*(9.D0*ctsq-1.D0)
      Ylmr(6,3) = (SQRT(2310.D0)/8.D0)*ct*(3.D0*ctsq-1.D0)
      Ylmr(6,2) = -(SQRT(330.D0)/16.D0)*((21.D0*ctsq-14.D0)*ctsq+1.D0)
      Ylmr(6,1) = (SQRT(11.D0)/8.D0)*ct*((63.D0*ctsq-70.D0)*ctsq+15.D0)
      Ylmr(7,7) = SQRT(3003.D0)/32.D0
      Ylmr(7,6) = -(3.D0*SQRT(1001.D0)/16.D0)*ct
      Ylmr(7,5) = (3.D0*SQRT(182.D0)/32.D0)*(11.D0*ctsq-1.D0)
      Ylmr(7,4) = -(SQRT(1365.D0)/16.D0)*ct*(11.D0*ctsq-3.D0)
      Ylmr(7,3) = (SQRT(1365.D0)/32.D0)*((33.D0*ctsq-18.D0)*ctsq+1.D0)
      Ylmr(7,2) = -(SQRT(546.D0)/16.D0)*ct*((33.D0*ctsq-30.D0)*
     &            ctsq+5.D0)
      Ylmr(7,1) = (SQRT(13.D0)/16.D0)*(((231.D0*ctsq-315.D0)*ctsq+
     &            105.D0)*ctsq-5.D0)
      Ylmr(8,8) = -3.D0*SQRT(1430.D0)/64.D0
      Ylmr(8,7) = (3.D0*SQRT(5005.D0)/32.D0)*ct
      Ylmr(8,6) = -(3.D0*SQRT(770.D0)/64.D0)*(13.D0*ctsq-1.D0)
      Ylmr(8,5) = (3.D0*SQRT(770.D0)/32.D0)*(13.D0*ctsq-3.D0)*ct
      Ylmr(8,4) = -(3.D0*SQRT(70.D0)/64.D0)*((143.D0*ctsq-66.D0)*
     &            ctsq+3.D0)
      Ylmr(8,3) = (3.D0*SQRT(35.D0)/32.D0)*((143.D0*ctsq-110.D0)*ctsq+
     &            15.D0)*ct
      Ylmr(8,2) = -(SQRT(210.D0)/64.D0)
     &            *(((429.D0*ctsq-495.D0)*ctsq+135.D0)*ctsq-5.D0)
      Ylmr(8,1) = (SQRT(15.D0)/16.D0)
     &            *(((429.D0*ctsq-693.D0)*ctsq+315.D0)*ctsq-35.D0)*ct
      Ylmr(9,9) = 3.D0*SQRT(24310.D0)/256.D0
      Ylmr(9,8) = -(3.D0*SQRT(24310.D0)/64.D0)*ct
      Ylmr(9,7) = (SQRT(7293.D0)/64.D0)*(15.D0*ctsq-1.D0)
      Ylmr(9,6) = -(3.D0*SQRT(34034.D0)/64.D0)*(5.D0*ctsq-1.D0)*ct
      Ylmr(9,5) = (3.D0*SQRT(2618.D0)/128.D0)*((65.D0*ctsq-26.D0)*
     &            ctsq+1.D0)
      Ylmr(9,4) = -(SQRT(39270.D0)/64.D0)*((39.D0*ctsq-26.D0)*ctsq+
     &            3.D0)*ct
      Ylmr(9,3) = (3.D0*SQRT(595.D0)/64.D0)
     &            *(((143.D0*ctsq-143.D0)*ctsq+33.D0)*ctsq-1.D0)
      Ylmr(9,2) = -(3.D0*SQRT(34.D0)/64.D0)
     &            *(((715.D0*ctsq-1001.D0)*ctsq+385.D0)*ctsq-35.D0)*ct
      Ylmr(9,1) = (SQRT(17.D0)/128.D0)
     &            *((((6435.D0*ctsq-12012.D0)*ctsq+6930.D0)*ctsq-
     &            1260.D0)*ctsq+35.D0)
      DO l = 2 , 9
         Ylmr(l,1) = Ylmr(l,1)/4.D0/pi
         DO m = 2 , l
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)/4.D0/pi
         ENDDO
      ENDDO
      END
