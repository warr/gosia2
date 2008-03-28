 
C----------------------------------------------------------------------
C SUBROUTINE GAMATT
C
C Called by: QFIT
C Calls:     GCF
C
C Purpose: calculate gamma attenuation in absorbers
C
C Formal parameters:
C      Qui    - attenuation
C      Tau1   - table of absorption coefficients
C      Tau2   - table of absorption coefficients
C      Xl1    - thickness of each kind of absorber
C      Nl     - number of kinds of absorber, we can treat (7)
C
C Note the absorbers are: Al, C, Fe, Cu, Ag/Cd/Sn, Ta and Pb, respectively.
 
      SUBROUTINE GAMATT(Qui,Tau1,Tau2,Xl1,Nl)
      IMPLICIT NONE
      INTEGER*4 i , i1 , k , Nl
      REAL*8 q , Qui , tau , Tau1 , Tau2 , thing , thing1 , thing3 , Xl1
      DIMENSION Tau1(10) , Tau2(10,7) , Xl1(7) , thing3(10) , q(9) , 
     &          Qui(8,10)

      DO i = 1 , 10 ! Loop over energies
         i1 = 1
         thing3(i) = 0.
 50      thing1 = -Tau2(i,i1)*Xl1(i1) + thing3(i)
         i1 = i1 + 1
         thing3(i) = thing1
         IF ( i1.LE.Nl ) GOTO 50 ! Loop over Nl absorbers
      ENDDO

      DO i = 1 , 10 ! Loop over energies
         tau = Tau1(i)
         thing = thing3(i)
         CALL GCF(tau,thing,q)
         DO k = 2 , 9
            Qui(k-1,i) = q(k)
         ENDDO
      ENDDO
      END
