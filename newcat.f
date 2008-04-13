      SUBROUTINE NEWCAT(Iexp,Jidim)
      IMPLICIT NONE
      REAL*8 a , b , FXIS1 , FXIS2 , PARX , PARXM , q1 , q2 , QAPR , 
     &       wg , wl , XI , XIR , xp , xx , zt
      INTEGER*4 IAPR , Iexp , IPATH , ISEX , ist , istop , Jidim , k , 
     &          kk , LAMDA , LAMMAX , LDNUM , LEAD , MAGA , MULTI , n , 
     &          NDIM , ng , nl , NMAX
      INTEGER*4 NMAX1
      COMMON /MAP   / PARX(50,12,5) , PARXM(50,4,10,6) , XIR(6,50)
      COMMON /CXI   / XI(500)
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      Jidim = NMAX + 1
      IF ( MAGA(Iexp).NE.0 ) Jidim = 3*NMAX + 1
      ist = 1
      DO kk = 1 , 6
         IF ( MULTI(kk).NE.0 ) THEN
            istop = MULTI(kk) - 1 + ist
            DO k = ist , istop
               xx = ABS(XI(k))
               xx = xx/XIR(kk,Iexp)
               DO n = 1 , 7 , 3
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
                     QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n)
                     QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n)
                     IF ( IAPR(k,1).EQ.IAPR(k,2) ) THEN
                        QAPR(k,1,n) = 0.
                        QAPR(k,2,n) = QAPR(k,2,n)/2.
                     ENDIF
                  ENDIF
               ENDDO
               IF ( MAGA(Iexp).NE.0 ) THEN
                  DO n = 2 , 6
                     IF ( n.NE.4 ) THEN
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
                        QAPR(k,2,n) = QAPR(k,1,n)*q2*FXIS2(k,n)
                        QAPR(k,1,n) = QAPR(k,1,n)*q1*FXIS1(k,n)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            ist = istop + 1
         ENDIF
      ENDDO
      END
