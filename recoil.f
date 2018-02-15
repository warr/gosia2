 
C----------------------------------------------------------------------
C SUBROUTINE RECOIL
C
C Called by: ANGULA
C Calls:     ROTATE
C
C Purpose: correct for relativistic effects of recoiling nucleus.
C
C Formal parameters:
C      Alab   - matrix of (l,m) pairs in lab frame
C      Attl   - matrix of (l,m) pairs in rotated frame
C      Beta   - beta of recoil
C      Theta  - angle to rotate
C
C We transform into the frame of the recoiling nucleus, correct according to
C the method of Lesser and then rotate back to the laboratory frame.
 
      SUBROUTINE RECOIL(Alab,Attl,Beta,Theta)
      IMPLICIT NONE
      REAL*8 Alab , atemp , Attl , Beta , betasq , dum , hold , test , 
     &       Theta
      INTEGER*4 i , i1 , j , l , m
      DIMENSION Alab(9,9) , Attl(9,9) , atemp(16)
      
      hold = Alab(1,1)
      IF ( ABS(hold).LT.1.E-9 ) RETURN

C     Rotate into frame of recoiling nucleus
      CALL ROTATE(Alab,Attl,-Theta,7,2)

C     Correct for relativistic effects
      Attl(2,1) = (2./SQRT(15.))*(SQRT(5.)*Attl(1,1)-Attl(3,1))
      Attl(2,2) = -Attl(3,2)/SQRT(5.)
      Attl(4,1) = (4./SQRT(35.))*(3.*Attl(3,1)-SQRT(5.)*Attl(5,1))
      Attl(4,2) = (8.*SQRT(2.)*Attl(3,2)-5.*SQRT(3.)*Attl(5,2))
     &            /SQRT(35.)
      Attl(4,3) = (2./SQRT(7.))*(2.*Attl(3,3)-SQRT(3.)*Attl(5,3))
      Attl(4,4) = -Attl(5,4)
      Attl(6,1) = (10./SQRT(11.))*(Attl(5,1)-(3./SQRT(13.))*Attl(7,1))
      Attl(6,2) = (1./SQRT(11.))
     &            *(4.*SQRT(6.)*Attl(5,2)-5.*SQRT(35./13.)*Attl(7,2))
      Attl(6,3) = SQRT(4./11.)
     &            *(SQRT(21.)*Attl(5,3)-10.*SQRT(2./13.)*Attl(7,3))
      Attl(6,4) = SQRT(1./11.)*(8.*Attl(5,4)-15.*SQRT(3./13.)*Attl(7,4))
      Attl(6,5) = SQRT(4./11.)*(3.*Attl(5,5)-5.*SQRT(5./13.)*Attl(7,5))
      Attl(6,6) = -Attl(7,6)*SQRT(25./13.)
      Attl(8,1) = (56./SQRT(195.))*Attl(7,1)
      Attl(8,2) = (32./SQRT(65.))*Attl(7,2)
      Attl(8,3) = (8.*SQRT(3./13.))*Attl(7,3)
      Attl(8,4) = (16.*SQRT(2./39.))*Attl(7,4)
      Attl(8,5) = (8.*SQRT(11./65.))*Attl(7,5)
      Attl(8,6) = (16.*SQRT(2./65.))*Attl(7,6)
      Attl(8,7) = (8./SQRT(15.))*Attl(7,7)
      DO l = 2 , 8 , 2
         DO m = 1 , l
            Attl(l,m) = Beta*Attl(l,m)
         ENDDO
      ENDDO
      betasq = Beta*Beta
      IF ( betasq.GE.1.0E-10 ) THEN
         i1 = 0
         DO i = 1 , 7 , 2
            DO j = 1 , i
               i1 = i1 + 1
               atemp(i1) = Attl(i,j)
            ENDDO
         ENDDO
         dum = (2./5.)*SQRT(5.)*atemp(1) - (10./7.)*atemp(2) + (12./35.)
     &         *SQRT(5.)*atemp(5)
         Attl(3,1) = atemp(2) + betasq*dum
         dum = -(17./14.)*atemp(3) + (2./7.)*SQRT(6.)*atemp(6)
         Attl(3,2) = atemp(3) + betasq*dum
         dum = -(4./7.)*atemp(4) + (2./7.)*SQRT(3.)*atemp(7)
         Attl(3,3) = atemp(4) + betasq*dum
         dum = (8./7.)*SQRT(5.)*atemp(2) - (380./77.)*atemp(5)
     &         + (100./11.)*SQRT(1./13.)*atemp(10)
         Attl(5,1) = atemp(5) + betasq*dum
         dum = (20./21.)*SQRT(6.)*atemp(3) - (723./154.)*atemp(6)
     &         + (20./11.)*SQRT(70./39.)*atemp(11)
         Attl(5,2) = atemp(6) + betasq*dum
         dum = (20./21.)*SQRT(3.)*atemp(4) - (306./77.)*atemp(7)
     &         + (40./11.)*SQRT(14./39.)*atemp(12)
         Attl(5,3) = atemp(7) + betasq*dum
         dum = -(61./22.)*atemp(8) + (40./11.)*SQRT(3./13.)*atemp(13)
         Attl(5,4) = atemp(8) + betasq*dum
         dum = -(12./11.)*atemp(9) + (20./11.)*SQRT(5./13.)*atemp(14)
         Attl(5,5) = atemp(9) + betasq*dum
         dum = (210./11.)*SQRT(1./13.)*atemp(5) - (574./55.)*atemp(10)
         Attl(7,1) = atemp(10) + betasq*dum
         dum = (14./11.)*SQRT(210./13.)*atemp(6) - (1121./110.)
     &         *atemp(11)
         Attl(7,2) = atemp(11) + betasq*dum
         dum = (28./11.)*SQRT(42./13.)*atemp(7) - (104./11.)*atemp(12)
         Attl(7,3) = atemp(12) + betasq*dum
         dum = (84./11.)*SQRT(3./13.)*atemp(8) - (181./22.)*atemp(13)
         Attl(7,4) = atemp(13) + betasq*dum
         dum = (42./11.)*SQRT(5./13.)*atemp(9) - (358./55.)*atemp(14)
         Attl(7,5) = atemp(14) + betasq*dum
         Attl(7,6) = atemp(15)*(1.-(43./10.)*betasq)
         Attl(7,7) = atemp(16)*(1.-(8./5.)*betasq)
         Attl(9,1) = (672./5.)*SQRT(1./221.)*atemp(10)*betasq
         Attl(9,2) = (144./5.)*SQRT(21./221.)*atemp(11)*betasq
         Attl(9,3) = 36.*SQRT(12./221.)*atemp(12)*betasq
         Attl(9,4) = 24.*SQRT(22./221.)*atemp(13)*betasq
         Attl(9,5) = (144./5.)*SQRT(11./221.)*atemp(14)*betasq
         Attl(9,6) = (72./5.)*SQRT(2./17.)*atemp(15)*betasq
         Attl(9,7) = (24./5.)*SQRT(7./17.)*atemp(16)*betasq
      ENDIF

C     Rotate back into laboratory frame
      CALL ROTATE(Attl,Alab,Theta,9,1)
      test = ABS(1.0-Alab(1,1)/hold)
      IF ( test.GT.1.0E-07 ) THEN
         WRITE (22,99001) test
99001    FORMAT (' ERROR IN ROTATION',1X,1E10.3/)
      ENDIF
      END
