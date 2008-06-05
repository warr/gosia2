 
C----------------------------------------------------------------------
C SUBROUTINE TAPMA
C
C Called by: GOSIA
C
C Purpose: read parameters for sensitivity maps
C
C Uses global variables:
C      DS     -
C      XV     - energy meshpoints (sometimes theta meshpoints) where we calculate exact Coulex
C      YGN    - gamma yield calculated without correction to angular distribution from finite recoil distance
C      ZETA   - various coefficients
C
C Formal parameters:
C      Lx     -
C      Iske   -
C      Isko   -
C      Iskf   -
C      Nflr   -
C      Idr    -
C      Nco    -
C      Nft    -
C      Enb    -
C
C Note that unit 14 is used internally for the purpose of sensitivity
C maps.
 
      SUBROUTINE TAPMA(Lx,Iske,Isko,Iskf,Nflr,Idr,Nco,Nft,Enb)
      IMPLICIT NONE
      REAL*8 emn , emx , en0 , Enb , tmn , tmx , tta , 
     &       ZETA
      INTEGER*4 Idr , Iske , Iskf , Isko , j , jf , jj , js , k , 
     &          Lx , lx1 , LZETA , na , Nco , ne , nfil , nfilt , Nflr , 
     &          Nft
      INTEGER*4 ng , ng1 , ntt
      INCLUDE 'vlin.inc'
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      INCLUDE 'yteor.inc'

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
