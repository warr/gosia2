      SUBROUTINE PATH(Irld)
      IMPLICIT NONE
      REAL*8 CAT , spm , vl
      INTEGER*4 i , IPATH , Irld , ISMAX , isp , ist , j , MAGA , NDIM , 
     &          NMAX , NMAX1 , NSTART , NSTOP
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      spm = CAT(Irld,3)
      DO i = 2 , NMAX
         IPATH(i) = 0
         ist = NSTART(i)
         IF ( ist.NE.0 ) THEN
            isp = NSTOP(i)
            DO j = ist , isp
               vl = CAT(j,3)
               IF ( ABS(vl-spm).LT.1.E-6 ) GOTO 50
            ENDDO
         ENDIF
         GOTO 100
 50      IPATH(i) = j
 100  ENDDO
      IPATH(1) = Irld
      END
