      SUBROUTINE TAPMA(Lx,Iske,Isko,Iskf,Nflr,Idr,Nco,Nft,Enb)
      IMPLICIT NONE
      REAL*8 DS , DSE , DSG , emn , emx , en0 , Enb , tmn , tmx , tta , 
     &       XV , YGN , YGP , YV , ZETA , ZV
      INTEGER*4 Idr , IFMO , Iske , Iskf , Isko , j , jf , jj , js , k , 
     &          Lx , lx1 , LZETA , na , Nco , ne , nfil , nfilt , Nflr , 
     &          Nft
      INTEGER*4 ng , ng1 , ntt
      COMMON /VLIN  / XV(51) , YV(51) , ZV(20) , DSG(20) , DSE(20) , DS
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /YTEOR / YGN(500) , YGP(500) , IFMO
      Nft = 0
      nfilt = 0
      REWIND 14
      IF ( Iske.NE.0 ) THEN
 50      READ (14,*) ne , ntt , emn , emx , tmn , tmx , na , tmx , tmx , 
     &               tmx
         nfil = ne*ntt*na
         nfilt = nfilt + nfil
         DO j = 1 , nfil
            READ (14,*) lx1 , Enb , tta , ng , DS , (YGN(k),k=1,Idr)
         ENDDO
         IF ( nfilt.NE.Iske ) GOTO 50
      ENDIF
      IF ( Nco.EQ.0 ) RETURN
      READ (14,*) ne , ntt , emn , emx , tmn , tmx , na , tmx , tmx , 
     &            tmx
      IF ( Isko.NE.0 ) THEN
         DO j = 1 , Isko
            READ (14,*) lx1 , Enb , tta , ng , DS , (YGN(k),k=1,Idr)
         ENDDO
      ENDIF
      DO j = 1 , Nflr
         js = (j-1)*Idr + 1
         jf = js + Idr - 1
         READ (14,*) lx1 , Enb , tta , ng1 , DS , (ZETA(k),k=js,jf)
         IF ( lx1.NE.Lx ) Nft = 1
         IF ( Nft.EQ.1 ) GOTO 100
         XV(j) = tta/57.2957795
         IF ( Iskf.NE.0 .AND. j.NE.Nflr ) THEN
            DO jj = 1 , Iskf
               READ (14,*) lx1 , en0 , tta , ng , DS , (YGN(k),k=1,Idr)
            ENDDO
         ENDIF
      ENDDO
      RETURN
 100  WRITE (22,99001)
99001 FORMAT (10X///10X,'TAPE READ ERROR'/10X,'JOB ABORTED')
      END
