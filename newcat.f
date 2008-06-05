 
C----------------------------------------------------------------------
C SUBROUTINE NEWCAT
C
C Called by: APRAM
C Calls:     FXIS1, FXIS2
C
C Purpose: create a new catalog of matrix elements
C
C Uses global variables:
C      IAPR   - index of initial and final levels for each matrix element
C      MAGA   - number of magnetic substates in approximate calculation
C      MULTI  - number of matrix elements having given multipolarity
C      NMAX   - number of levels
C      PARX   -
C      PARXM  -
C      QAPR   - approximate Coulomb amplitudes
C      XI     - xi coupling coefficients
C      XIR    - [for maps]
C
C Formal parameters:
C     Iexp    - experiment number
C     Jidim   - 
      
      SUBROUTINE NEWCAT(Iexp,Jidim)
      IMPLICIT NONE
      REAL*8 a , b , FXIS1 , FXIS2 , q1 , q2 , QAPR , wg , wl , 
     &       xp , xx , zt
      INTEGER*4 IAPR , Iexp , ISEX , ist , istop , Jidim , k , 
     &          kk , n , ng , nl
      INCLUDE 'map.inc'
      INCLUDE 'cxi.inc'
      INCLUDE 'clcom.inc'
      INCLUDE 'coex2.inc'
      INCLUDE 'pth.inc'
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(75)

      Jidim = NMAX + 1
      IF ( MAGA(Iexp).NE.0 ) Jidim = 3*NMAX + 1
      ist = 1
      DO kk = 1 , 6 ! For each multipolarity E1 to E6
         IF ( MULTI(kk).NE.0 ) THEN ! If there are matrix elements for this multipolarity
            istop = MULTI(kk) - 1 + ist ! Last matrix element for this multipolarity
            DO k = ist , istop ! For each matrix element for this multipolarity
               xx = ABS(XI(k))
               xx = xx/XIR(kk,Iexp)

               DO n = 1 , 7 , 3 ! For n in 1, 4 and 7
                  IF ( MAGA(Iexp).NE.0 .OR. n.EQ.4 ) THEN
                     zt = QAPR(k,1,n)
                     zt = ABS(zt)
                     xp = 9.*xx
                     nl = INT(xp) + 1
                     wg = xp - DBLE(nl-1)
                     ng = nl + 1
                     wl = DBLE(nl) - xp
                     a = wg*PARXM(Iexp,1,ng,kk) + wl*PARXM(Iexp,1,nl,kk)
                     b = wg*PARXM(Iexp,2,ng,kk) + wl*PARXM(Iexp,2,nl,kk)
                     q1 = a*zt + b
                     a = wg*PARXM(Iexp,3,ng,kk) + wl*PARXM(Iexp,3,nl,kk)
                     b = wg*PARXM(Iexp,4,ng,kk) + wl*PARXM(Iexp,4,nl,kk)
                     q2 = a*zt + b
                     QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n) ! FXIS2 = 1
                     QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n) ! FXIS1 = sign(xi)
                     IF ( IAPR(k,1).EQ.IAPR(k,2) ) THEN
                        QAPR(k,1,n) = 0.
                        QAPR(k,2,n) = QAPR(k,2,n)/2.
                     ENDIF
                  ENDIF
               ENDDO ! Loop over n

               IF ( MAGA(Iexp).NE.0 ) THEN
                  DO n = 2 , 6
                     IF ( n.NE.4 ) THEN ! For N in 2, 3, 5, 6
                        zt = QAPR(k,1,n)
                        zt = ABS(zt)
                        xp = 4.*xx
                        nl = INT(xp) + 1
                        wg = xp - DBLE(nl-1)
                        ng = nl + 1
                        wl = DBLE(nl) - xp
                        q1 = wg*PARX(Iexp,2*kk-1,ng)
     &                       + wl*PARX(Iexp,2*kk-1,nl)
                        q2 = wg*PARX(Iexp,2*kk,ng)
     &                       + wl*PARX(Iexp,2*kk,nl)
                        QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n) ! FXIS2 = sign(xi)
                        QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n) ! FXIS1 = 1
                     ENDIF
                  ENDDO ! Loop over n
               ENDIF
            ENDDO ! Loop over matrix elements for this mult. k
            ist = istop + 1
         ENDIF
      ENDDO ! Loop over electric multipolarities kk
      END
