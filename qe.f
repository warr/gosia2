 
C----------------------------------------------------------------------
C SUBROUTINE QE
C
C Called by: SNAKE
C
C Purpose: calculate Qe values
C
C Formal parameters:
C      C     - cosh(omega) + epsilon
C      D     - sqrt(epsilon^2 - 1) * sinh(omega)
C      B2    - B^2 = (epsilon * cosh(omega) + 1)^2
C      C2    - C^2
C      D2    - D^2
C      B4    - B^4
C      B6    - B^6
C      D3    - D^3
C      B8    - B^8
C      C4    - C^4
C      D4    - D^4
C      B10   - B^10
C      D5    - D^5
C      B12   - B^12
C      D6    - D^6
C      Lmda  - lambda
C      Pol   - E1 polarisation factor to allow for GDR excitation
C      Cq    - array where the results are returned
C
C We used different formulae depending on lambda (see the table of electric
C collision functions in the gosia manual).
C
C Note that we multiply by Pol, which is the E1 dipole correction factor for
C the GDR in the case of the quadrupole. This is correct. It isn't really a
C quadrupole, but the operator has the same form as the E2 operator, so it
C is easiest to consider the effect as a correction to this operator.
C
C Lmda = lambda (1 = E1, 2 = E2... 6 = E6)
 
      SUBROUTINE QE(C,D,B2,C2,D2,B4,B6,D3,B8,C4,D4,B10,D5,B12,D6,Lmda,
     &              Pol,Cq)
      IMPLICIT NONE
      REAL*8 B10 , B12 , B2 , B4 , B6 , B8 , C , C2 , C4 , Cq , D , D2 ,
     &       D3 , D4 , D5 , D6 , Pol
      INTEGER*4 Lmda
      DIMENSION Cq(7)

      IF ( Lmda.EQ.2 ) THEN ! E2
         Cq(1) = 0.75*(2.0*C2-D2)/B4*Pol
         Cq(2) = -1.83711730*C*D/B4*Pol
         Cq(3) = -0.91855865*D2/B4*Pol
         RETURN
      ELSEIF ( Lmda.EQ.3 ) THEN ! E3
         Cq(1) = 1.875*C*(2.0*C2-3.0*D2)/B6
         Cq(2) = -1.62379763*(4.0*C2-D2)*D/B6
         Cq(3) = -5.13489890*C*D2/B6
         Cq(4) = 2.09631373*D3/B6
         RETURN
      ELSEIF ( Lmda.EQ.4 ) THEN ! E4
         Cq(1) = 1.09375000*(8.0*C4-24.0*C2*D2+3.0*D4)/B8
         Cq(2) = -4.89139867*C*(4.0*C2-3.0*D2)*D/B8
         Cq(3) = -3.45874113*(6.0*C2-D2)*D2/B8
         Cq(4) = 12.9414244*C*D3/B8
         Cq(5) = 4.57548440*D4/B8
         RETURN
      ELSEIF ( Lmda.EQ.5 ) THEN ! E5
         Cq(1) = 1.230468*C*(-14.*C2*(9.*D2+B2)+30.*B4)/B10
         Cq(2) = -1.347911*D*(35.*C2*(-3.*D2+B2)+5.*B4)/B10
         Cq(3) = -35.662372*D2*C*(-3.*D2+2.*B2)/B10
         Cq(4) = 7.279552*D3*(9.*C2-B2)/B10
         Cq(5) = 30.884521*D4*C/B10
         Cq(6) = -9.766543*D5/B10
         RETURN
      ELSEIF ( Lmda.EQ.6 ) THEN ! E6
         Cq(1) = 2.707031*(21.*C2*(-C2*(11.*D2+4.*B2)+5.*B4)-5.*B6)/B12
         Cq(2) = -17.543567*D*C*(3.*C2*(-11.*D2+B2)+5.*B4)/B12
         Cq(3) = -13.869408*D2*(3.*C2*(-11.*D2+5.*B2)+B4)/B12
         Cq(4) = 27.738815*D3*C*(-11.*D2+8.*B2)/B12
         Cq(5) = 15.193177*D4*(11.*C2-B2)/B12
         Cq(6) = -71.262308*D5*C/B12
         Cq(7) = -20.571656*D6/B12
         GOTO 99999
      ENDIF ! E1
      Cq(1) = 0.5*C/B2
      Cq(2) = -0.35355339*D/B2
99999 END
