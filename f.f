 
C----------------------------------------------------------------------
C FUNCTION F
C
C Called by: SEQ
C Calls:     WSIXJ, WTHREJ
C
C Purpose: evaluates the F coefficients.
C
C Formal coefficients:
C      K      - K
C      Sji    - initial spin
C      Sjf    - final spin
C      L1     - lambda
C      L2     - lambda'
C
C Return value:
C      F-coefficient
C
C We evaluate:
C F_k(\lambda \lambda^\prime I_2 I1) = (-1)^{I_1 + I_2 -l} *
C        \sqrt{(2 k + 1) (2 I_1 + 1) (2 \lambda + 1) (2 \lambda^\prime + 1) *
C        \threej{\lambda \lambda^\prime k 1 -1 0} *
C        \sixj{\lambda \lambda^\prime k I_1 I_1 I_2}
C
C Here \lambda = L1, \lambda^\prime = L2, I_1 = Sji, I_2 = Sjf, k = K
C
C Note that the code actually evaluates:
C \sixj{I_1 I_1 k \lambda^\prime \lambda I_2} which is equal to
C \sixj{\lambda \lambda^\prime k I_1 I_1 I_2} by the symmetry rules for 6-j
C symbols.
C
C Note also that both WTHREJ and WSIXJ need to have parameters which are twice
C the values to calculate, so that they can handle half-integers correctly.

      REAL*8 FUNCTION F(K,Sji,Sjf,L1,L2)
      IMPLICIT NONE
      INTEGER*4 ix , jfz , jiz , K , kz , l , L1 , l1z , L2 , l2z
      REAL*8 phase , Sjf , Sji , WSIXJ , WTHREJ
      
      F = 0.
      IF ( (L1*L2).EQ.0 ) RETURN
      ix = INT(Sji+Sjf+.0001)
      l = ix - 1
      phase = 1.
      IF ( l/2*2.NE.l ) phase = -1.
      kz = K*2
      jiz = INT(Sji*2)
      jfz = INT(Sjf*2)
      l1z = L1*2
      l2z = L2*2
      F = phase*SQRT((l1z+1.)*(l2z+1.)*(jiz+1.)*(kz+1.))
     &    *WTHREJ(l1z,l2z,kz,2,-2,0)*WSIXJ(jiz,jiz,kz,l2z,l1z,jfz)
      END
