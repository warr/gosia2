 
C----------------------------------------------------------------------
C SUBROUTINE TENS
C
C Called by: FTBM, GOSIA
C Calls:     DJMM
C
C Purpose:
C
C Uses global variables:
C      IAXS   - axial symmetry flag
C      IEXP   - experiment number
C      NMAX   - number of levels
C      TETACM - theta of particle detector in center of mass frame
C      ZETA   - various coefficients
C
C Formal parameter:
C      Bten   - 
 
      SUBROUTINE TENS(Bten)
      IMPLICIT NONE
      REAL*8 arg , Bten , DJMM
      INTEGER*4 i , ind , inz , iph , ix , k , k1 , kp , l , lp , 
     &          lpp , lx , lxx
      DIMENSION Bten(1200)
      INCLUDE 'kin.inc'
      INCLUDE 'tcm.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'coex2.inc'

      ix = NMAX*28
      arg = 1.570796327 + TETACM(IEXP)/2.
      DO i = 1 , ix
         ZETA(i) = 0.
      ENDDO

      DO i = 2 , NMAX ! For each level except ground state
         DO kp = 1 , 7 , 2
            k = kp - 1
            k1 = INT(DBLE(k)/2.+.01)
            IF ( k.EQ.0 ) THEN
               ind = (i-2)*16 + 1
               inz = (i-1)*28 + 1
               ZETA(inz) = Bten(ind)
            ELSE
               DO lp = 1 , kp
                  IF ( IAXS(IEXP).NE.0 .OR. lp.EQ.1 ) THEN
                     inz = (i-1)*28 + k1*7 + lp
                     l = lp - 1
                     DO lpp = 1 , kp
                        ind = k*k/4 + lpp + (i-2)*16
                        lx = lpp - 1
                        lxx = lx
 2                      iph = (-1)**(l+INT(DBLE(lxx)/2.))
                        ZETA(inz) = ZETA(inz) + Bten(ind)
     &                              *iph*DJMM(arg,k,lx,l)
                        IF ( lpp.NE.1 ) THEN
                           IF ( lx.GE.0 ) THEN
                              lx = -lx
                              lxx = lx - 1
                              GOTO 2
                           ENDIF ! if lx .ge. 0
                        ENDIF ! if lpp .ne. 1
                     ENDDO ! Loop over lpp
                  ENDIF ! if iaxs .ne.0 .or. lp.eq.1
               ENDDO ! Loop over lp
            ENDIF ! If k .eq. 0
         ENDDO ! Loop over kp
      ENDDO ! Loop over level i
      END
