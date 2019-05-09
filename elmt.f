
C----------------------------------------------------------------------
C FUNCTION ELMT
C
C Called by: GOSIA
C Calls:     WTHREJ
C
C Purpose: collective model matrix elements (OP,THEO)
C
C Formal parameters:
C      Xi1    - spin of initial level
C      Xi2    - spin of final level
C      Lam    - multipolarity
C      Nb1    - band number of initial level
C      Nb2    - band number of final level
C      Xk1    - initial level
C      Xk2    - final level
C      Xm1    - intrinsic moment Q1
C      Xm2    - intrinsic moment Q2
C      Xm3    - intrinsic moment Q3
C
C Return value:
C      Collective model matrix element
C
C Note that the parameters to WTHREJ are doubled to allow it to handle half
C integers.

      REAL*8 FUNCTION ELMT(Xi1,Xi2,Lam,Nb1,Nb2,Xk1,Xk2,Xm1,Xm2,Xm3)
      IMPLICIT NONE
      REAL*8 addt , fac , fct , pha1 , pha2 , s1 , s2 , WTHREJ , Xi1 ,
     &       Xi2 , Xk1 , Xk2 , xlam , Xm1 , Xm2 , Xm3 , xn
      INTEGER*4 i1 , i2 , ipha , k1 , k2 , l , la , Lam , llam , n ,
     &          Nb1 , Nb2

      la = Lam
      IF ( la.GT.6 ) la = la - 6
      xlam = DBLE(la)
      i1 = INT(2.*Xi1)
      i2 = INT(2.*Xi2)
      llam = 2*la
      k1 = INT(2.*Xk1)
      k2 = INT(2.*Xk2)
      fac = SQRT(2.*Xi1+1.)*SQRT(2.*Xi2+1.)
C-----In-band matrix element

      IF ( Nb1.NE.Nb2 ) THEN
C-----Interband, K-allowed
C-----One K=0
         IF ( ABS(k1-k2).GE.llam ) THEN
C-----Forbidden and K1-K2=lambda, Mikhailov formula
            addt = 0.
            IF ( k1.EQ.1 ) addt = (-1.)**((i1+1)/2)*(i1+1)/2.*Xm3
            xn = ABS(Xk1-Xk2) - xlam
            n = INT(xn+.1)
            IF ( n.EQ.0 ) THEN
               fct = 1.
            ELSEIF ( n.EQ.1 ) THEN
               fct = SQRT((Xi1-Xk1)*(Xi1+Xk1+1.))
            ELSE
               s1 = Xi1 - Xk1
               s2 = Xi1 + Xk1 + 1.
               DO l = 1 , n
                  s1 = s1*(Xi1-Xk1-DBLE(l))
                  s2 = s2*(Xi1+Xk2+1.+DBLE(l))
               ENDDO
               fct = SQRT(s1*s2)
            ENDIF
            pha1 = (-1.)**INT((Xi1-xlam+Xk2)+.1)
            ELMT = fac*pha1*fct*WTHREJ(i1,llam,i2,k2-llam,llam,-k2)
     &             *(Xm1+Xm2*(Xi2*(Xi2+1.)-Xi1*(Xi1+1.))+addt)
         ELSEIF ( k1.NE.0 .AND. k2.NE.0 ) THEN
C-----Both K's non-zero
            pha1 = (-1.)**((i1-llam+k2)/2)
            pha2 = (-1.)**((i1+k1)/2)*pha1
            ELMT = fac*(pha1*WTHREJ(i1,llam,i2,k1,k2-k1,-k2)
     &             *Xm1+pha2*WTHREJ(i1,llam,i2,-k1,k1+k2,-k2)*Xm2)
            RETURN
         ELSE
            ipha = (i1-llam+k2)/2
            IF ( k2.EQ.0 ) ipha = ((i2-llam+k1)/2)
            pha1 = (-1.)**ipha
            ELMT = fac*pha1*WTHREJ(i1,llam,i2,0,k2,-k2)*Xm1
            IF ( k2.EQ.0 ) ELMT = fac*pha1*WTHREJ(i2,llam,i1,0,k1,-k1)
     &                            *Xm1
            IF ( k1.NE.0 .OR. k2.NE.0 ) ELMT = ELMT*SQRT(2.)
            RETURN
         ENDIF
C-----K=0
      ELSEIF ( k1.NE.0 ) THEN
C-----In band, K.ne.0
         pha1 = (-1.)**((i1-llam+k1)/2)
         pha2 = (-1.)**((k1+i1)/2+1)*pha1
         ELMT = fac*(pha1*WTHREJ(i1,llam,i2,k1,0,-k1)
     &          *Xm1+pha2*WTHREJ(i1,llam,i2,-k1,2*k1,-k1)*Xm2)
         RETURN
      ELSE
         ELMT = fac*WTHREJ(i1,llam,i2,0,0,0)*Xm1
         RETURN
      ENDIF
      END
