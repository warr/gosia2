
C----------------------------------------------------------------------
C SUBROUTINE QFIT
C
C Called by: GOSIA
C Calls:     GAMATT
C
C Purpose: for OP,GDET, fit attenuation by absorbers
C
C Formal parameters:
C      Qui    - attenuation coefficients
C      Tau1   - absorption coefficients' table
C      Tau2   - absorption coefficients' table
C      Eng    - gamma energy
C      Xl1    - thickness of absorbers
C      Cf     - coefficients of fit
C      Nl     - number of types of absorber (7)
C      Ind    - type of absorber
C
C Note the absorbers are: Al, C, Fe, Cu, Ag/Cd/Sn, Ta and Pb, respectively.

      SUBROUTINE QFIT(Qui,Tau1,Tau2,Eng,Xl1,Cf,Nl,Ind)
      IMPLICIT NONE
      REAL*8 ca , cb , Cf , cm , cn , co , d , d1 , d2 , Eng , Qui ,
     &       Tau1 , Tau2 , Xl1
      INTEGER*4 Ind , ind1 , k , Nl
      DIMENSION Tau1(10) , Eng(10) , Tau2(10,7) , Xl1(7) , Qui(8,10) ,
     &          Cf(8,2)

      CALL GAMATT(Qui,Tau1,Tau2,Xl1,Nl)

      ind1 = 5
      IF ( Ind.EQ.4 ) ind1 = 6
      IF ( Ind.EQ.5 ) ind1 = 7
      DO k = 1 , 8
         co = Qui(k,Ind)
         cn = Qui(k,10)
         cm = Qui(k,ind1)
         ca = (Eng(ind1)-Eng(Ind))**2
         cb = (Eng(10)-Eng(Ind))**2
         d = ca*(co-cn) - cb*(co-cm)
         d1 = ca*cm*(co-cn) - cb*cn*(co-cm)
         d2 = ca*cb*(cn-cm)
         Cf(k,1) = d1/d
         Cf(k,2) = d2/d
      ENDDO
      END
