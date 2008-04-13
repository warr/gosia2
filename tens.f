      SUBROUTINE TENS(Bten)
      IMPLICIT NONE
      REAL*8 arg , Bten , DJMM , DSIGS , EPS , EROOT , FIEX , TETACM , 
     &       TREP , ZETA
      INTEGER*4 i , IAXS , IEXP , ind , inz , iph , ix , k , k1 , kp , 
     &          l , lp , lpp , lx , lxx , LZETA , NDIM , NMAX , NMAX1
      DIMENSION Bten(1200)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /TCM   / TETACM(50) , TREP(50) , DSIGS(50)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      ix = NMAX*28
      arg = 1.570796327 + TETACM(IEXP)/2.
      DO i = 1 , ix
         ZETA(i) = 0.
      ENDDO
      DO i = 2 , NMAX
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
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ENDDO
      END
