
C----------------------------------------------------------------------
C SUBROUTINE ROTATE
C
C Called by: RECOIL
C Calls:     GOSIA, ROTATE, TENS
C
C Purpose: rotate the frame of reference
C
C Formal parameters:
c      Alab   - matrix in lab frame
C      Attl   - matrix rotated
C      Theta  - angle to rotate by
C      K2     - maximum dimension to rotate
C      Kd     - step over dimension
C
C We use:
C \rho_{k \xi} = (-1)^\xi *
C  \sum_\xi^\prime{i^\xi\prime d^k_{xi^\prime \xi}
C     ({\pi + \theta \over 2}) \rho_{k \xi^\prime}}

      SUBROUTINE ROTATE(Alab,Attl,Theta,K2,Kd)
      IMPLICIT NONE
      REAL*8 Alab , Attl , djarg , DJMM , dkkk , sum , Theta
      INTEGER*4 idj , idm , idmp , j , k , K2 , ka , kappa , kapri , Kd
      DIMENSION Alab(9,9) , Attl(9,9)

      IF ( ABS(Theta).GT..01 ) THEN
         djarg = Theta
         DO ka = 1 , K2 , Kd
            idj = ka - 1
            DO kappa = 1 , ka
               idmp = kappa - 1
               sum = 0.0
               DO kapri = 1 , ka
                  idm = kapri - 1
                  dkkk = DJMM(djarg,idj,idm,idmp)
                  sum = sum + dkkk*Alab(ka,kapri)
               ENDDO
               IF ( ka.NE.1 ) THEN
                  DO kapri = 2 , ka
                     idm = -kapri + 1
                     dkkk = DJMM(djarg,idj,idm,idmp)
                     sum = sum + dkkk*Alab(ka,kapri)*(-1.0)**(kapri-1)
                  ENDDO
               ENDIF
               Attl(ka,kappa) = sum
            ENDDO
         ENDDO
         GOTO 99999
      ENDIF
      DO j = 1 , 9
         DO k = 1 , 9
            Attl(j,k) = Alab(j,k)
         ENDDO
      ENDDO
99999 END
