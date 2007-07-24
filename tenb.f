 
C----------------------------------------------------------------------
C SUBROUTINE TENB
C
C Called by: FTBM, GOSIA
C Calls:     WTHREJ
C
C Purpose: calculate the state of polarization of the decaying level
C
C Uses global variables:
C      ARM    - reduced matrix elements
C      CAT    -
C      NMAX   - number of levels
C      NSTART -
C      NSTOP  -
C      SPIN   - spin of level
C
C Formal parameters:
C      Icl    -
C      Bten   -
C      Lmax   -
C
C Note that the parameters to WTHREJ are all doubled, so that this routine
C can cope with half-integers.

      SUBROUTINE TENB(Icl,Bten,Lmax)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , Bten , CAT , ce , DIPOL , EN , fc , si , 
     &       SPIN , WTHREJ , x , ZPOL
      INTEGER*4 i , Icl , iha , ila , ilg , ind , isi , ISMAX , ISO , 
     &          ite , jm , jmp , k , kk , kp , l , ll , Lmax , lp , m
      INTEGER*4 mm , mp , ms , msp , NDIM , NMAX , NMAX1 , NSTART , 
     &          NSTOP
      COMPLEX*16 ARM
      DIMENSION Bten(1200)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /AZ    / ARM(600,7)
      
      iha = (-1)**INT(2.*SPIN(1)+.01)
      IF ( Icl.EQ.1 ) THEN
         ms = 16*(NMAX-1)
         DO i = 1 , ms
            Bten(i) = 0.
         ENDDO
      ENDIF

      DO i = 2 , NMAX
         ms = NSTART(i)
         IF ( ms.NE.0 ) THEN
            msp = NSTOP(i)
            si = SPIN(i)
            isi = INT(2.*si+.01)
            ce = SQRT(2.*si+1.)
            DO kp = 1 , 7 , 2
               k = kp - 1
               kk = 2*k
               IF ( isi.GE.k ) THEN
                  ila = -1
                  DO lp = 1 , kp
                     ila = -ila
                     l = lp - 1
                     ll = 2*l
                     ind = k*k/4 + lp + (i-2)*16
                     DO m = ms , msp
                        mm = m
                        mp = m + l
                        jm = INT(2.01*CAT(mm,3))
                        IF ( mp.GT.NSTOP(i) ) GOTO 4
                        ilg = (-1)**INT(si-CAT(mp,3))
                        jmp = -INT(2.01*CAT(mp,3))
                        fc = WTHREJ(isi,kk,isi,jmp,ll,jm)
                        ite = 1
 2                      IF ( ila.EQ.1 ) x = DBLE(ARM(mp,5))
     &                       *DBLE(ARM(mm,5)) + IMAG(ARM(mp,5))
     &                       *IMAG(ARM(mm,5))
                        IF ( ila.NE.1 ) x = DBLE(ARM(mp,5))
     &                       *IMAG(ARM(mm,5)) - DBLE(ARM(mm,5))
     &                       *IMAG(ARM(mp,5))
                        Bten(ind) = Bten(ind) + x*fc*ilg
                        IF ( ite.EQ.2 ) GOTO 6
 4                      IF ( iha.NE.1 .OR. Icl.NE.Lmax ) THEN
                           ite = 2
                           mp = mp - 2*l
                           IF ( mp.GE.NSTART(i) ) THEN
                              jmp = INT(2.01*CAT(mp,3))
                              jm = -jm
                              fc = WTHREJ(isi,kk,isi,jmp,ll,jm)
                              ilg = (-1)**INT(si+CAT(mp,3))
                              GOTO 2
                           ENDIF
                        ENDIF
 6                   ENDDO
                     IF ( Icl.EQ.Lmax ) Bten(ind) = Bten(ind)
     &                    *ce/(2.*SPIN(1)+1.)
                  ENDDO ! Loop over lp
               ENDIF ! If isi.GE.k
            ENDDO ! Loop over kp
         ENDIF ! If ms.NE.0
      ENDDO ! Loop over i
      END
