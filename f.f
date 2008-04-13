      FUNCTION F(K,Sji,Sjf,L1,L2)
      IMPLICIT NONE
      REAL*8 F , phase , Sjf , Sji , WSIXJ , WTHREJ
      INTEGER*4 ix , jfz , jiz , K , kz , l , L1 , l1z , L2 , l2z
      F = 0.
      IF ( (L1*L2).EQ.0 ) RETURN
      ix = INT(Sji+Sjf+.0001)
      l = ix - 1
      phase = 1.
      IF ( l/2*2.NE.l ) phase = -1.
      kz = K*2
      jiz = Sji*2
      jfz = Sjf*2
      l1z = L1*2
      l2z = L2*2
      F = phase*SQRT((l1z+1.)*(l2z+1.)*(jiz+1.)*(kz+1.))
     &    *WTHREJ(l1z,l2z,kz,2,-2,0)*WSIXJ(jiz,jiz,kz,l2z,l1z,jfz)
      END
