
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
      DIMENSION Ylmr(9,9) , st(7)

      ct = COS(Theta)
      ctsq = ct*ct
      IF ( IAXS(IEXP).EQ.0 ) THEN ! If axially symmetric
         Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
         Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
         Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
         GOTO 99999
      ENDIF
      st(1) = SIN(Theta)
      DO i = 2 , 7
         j = i - 1
         st(i) = st(j)*st(1)
      ENDDO
      Ylmr(1,3) = .1089659406
      Ylmr(1,2) = -.2179318812*ct
      Ylmr(1,1) = .0889703179*(3.*ctsq-1.)
      Ylmr(2,5) = .1248361677
      Ylmr(2,4) = -.3530900028*ct
      Ylmr(2,3) = .0943672726*(7.*ctsq-1.)
      Ylmr(2,2) = -.1334554768*ct*(7.*ctsq-3.)
      Ylmr(2,1) = .0298415518*((35.*ctsq-30.)*ctsq+3.)
      Ylmr(3,7) = .1362755124
      Ylmr(3,6) = -.4720722226*ct
      Ylmr(3,5) = .100646136*(11.*ctsq-1.)
      Ylmr(3,4) = -.1837538634*ct*(11.*ctsq-3.)
      Ylmr(3,3) = .0918769316*((33.*ctsq-18.)*ctsq+1.)
      Ylmr(3,2) = -.1162161475*ct*((33.*ctsq-30.)*ctsq+5.)
      Ylmr(3,1) = .0179325408*(((231.*ctsq-315.)*ctsq+105.)*ctsq-5.)
      DO l = 1 , 3
         lf = 2*l + 1
         DO m = 2 , lf
            Ylmr(l,m) = Ylmr(l,m)*st(m-1)
         ENDDO
      ENDDO
99999 END
