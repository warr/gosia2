      FUNCTION GF(K,Sji,Sjf,L)
      IMPLICIT NONE
      REAL*8 GF , phase , Sjf , Sji , WSIXJ
      INTEGER*4 i , ix , jfz , jiz , K , kz , L , lz
      GF = 0.
      IF ( L.EQ.0 ) RETURN
      ix = INT(Sji+Sjf+.0001)
      i = ix + L + K
      phase = 1.
      IF ( i/2*2.NE.i ) phase = -1.
      kz = K*2
      jiz = Sji*2
      jfz = Sjf*2
      lz = L*2
      GF = phase*SQRT((jiz+1.)*(jfz+1.))*WSIXJ(jiz,jiz,kz,jfz,jfz,lz)
      END