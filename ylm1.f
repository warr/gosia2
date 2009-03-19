 
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
      Ylmr(2,2) = -SQRT(6.)/2.
      Ylmr(2,1) = SQRT(3.)*ct
      Ylmr(3,3) = SQRT(30.)/4.
      Ylmr(3,2) = -(SQRT(30.)/2.)*ct
      Ylmr(3,1) = (SQRT(5.)/2.)*(3.*ctsq-1.)
      Ylmr(4,4) = -SQRT(35.)/4.
      Ylmr(4,3) = (SQRT(210.)/4.)*ct
      Ylmr(4,2) = -(SQRT(21.)/4.)*(5.*ctsq-1.)
      Ylmr(4,1) = (SQRT(7.)/2.)*ct*(5.*ctsq-3.)
      Ylmr(5,5) = 3.*SQRT(70.)/16.
      Ylmr(5,4) = -(3.*SQRT(35.)/4.)*ct
      Ylmr(5,3) = (3.*SQRT(10.)/8.)*(7.*ctsq-1.)
      Ylmr(5,2) = -(3.*SQRT(5.)/4.)*ct*(7.*ctsq-3.)
      Ylmr(5,1) = (3./8.)*((35.*ctsq-30.)*ctsq+3.)
      Ylmr(6,6) = -3.*SQRT(77.)/16.
      Ylmr(6,5) = (3.*SQRT(770.)/16.)*ct
      Ylmr(6,4) = -(SQRT(385.)/16.)*(9.*ctsq-1.)
      Ylmr(6,3) = (SQRT(2310.)/8.)*ct*(3.*ctsq-1.)
      Ylmr(6,2) = -(SQRT(330.)/16.)*((21.*ctsq-14.)*ctsq+1.)
      Ylmr(6,1) = (SQRT(11.)/8.)*ct*((63.*ctsq-70.)*ctsq+15.)
      Ylmr(7,7) = SQRT(3003.)/32.
      Ylmr(7,6) = -(3.*SQRT(1001.)/16.)*ct
      Ylmr(7,5) = (3.*SQRT(182.)/32.)*(11.*ctsq-1.)
      Ylmr(7,4) = -(SQRT(1365.)/16.)*ct*(11.*ctsq-3.)
      Ylmr(7,3) = (SQRT(1365.)/32.)*((33.*ctsq-18.)*ctsq+1.)
      Ylmr(7,2) = -(SQRT(546.)/16.)*ct*((33.*ctsq-30.)*ctsq+5.)
      Ylmr(7,1) = (SQRT(13.)/16.)*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
      Ylmr(8,8) = -3.*SQRT(1430.)/64.
      Ylmr(8,7) = (3.*SQRT(5005.)/32.)*ct
      Ylmr(8,6) = -(3.*SQRT(770.)/64.)*(13.*ctsq-1.)
      Ylmr(8,5) = (3.*SQRT(770.)/32.)*(13.*ctsq-3.)*ct
      Ylmr(8,4) = -(3.*SQRT(70.)/64.)*((143.*ctsq-66.)*ctsq+3.)
      Ylmr(8,3) = (3.*SQRT(35.)/32.)*((143.*ctsq-110.)*ctsq+15.)*ct
      Ylmr(8,2) = -(SQRT(210.)/64.)
     &            *(((429.*ctsq-495.)*ctsq+135.)*ctsq-5.)
      Ylmr(8,1) = (SQRT(15.)/16.)
     &            *(((429.*ctsq-693.)*ctsq+315.)*ctsq-35.)*ct
      Ylmr(9,9) = 3.*SQRT(24310.)/256.
      Ylmr(9,8) = -(3.*SQRT(24310.)/64.)*ct
      Ylmr(9,7) = (SQRT(7293.)/64.)*(15.*ctsq-1.)
      Ylmr(9,6) = -(3.*SQRT(34034.)/64.)*(5.*ctsq-1.)*ct
      Ylmr(9,5) = (3.*SQRT(2618.)/128.)*((65.*ctsq-26.)*ctsq+1.)
      Ylmr(9,4) = -(SQRT(39270.)/64.)*((39.*ctsq-26.)*ctsq+3.)*ct
      Ylmr(9,3) = (3.*SQRT(595.)/64.)
     &            *(((143.*ctsq-143.)*ctsq+33.)*ctsq-1.)
      Ylmr(9,2) = -(3.*SQRT(34.)/64.)
     &            *(((715.*ctsq-1001.)*ctsq+385.)*ctsq-35.)*ct
      Ylmr(9,1) = (SQRT(17.)/128.)
     &            *((((6435.*ctsq-12012.)*ctsq+6930.)*ctsq-1260.)
     &            *ctsq+35.)
      DO l = 2 , 9
         Ylmr(l,1) = Ylmr(l,1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
         DO m = 2 , l
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
         ENDDO
      ENDDO
      END
