      SUBROUTINE LIMITS
      IMPLICIT NONE
      REAL*8 ELM , ELML , ELMU , SA
      INTEGER*4 IVAR , j , LMAXE , MAGEXC , MEMAX , MEMX6
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      DO j = 1 , MEMAX
         IF ( IVAR(j).NE.0 ) THEN
            IF ( ELM(j).GT.ELMU(j) .OR. ELM(j).LT.ELML(j) ) THEN
               IF ( ELM(j).GT.ELMU(j) ) THEN
                  ELM(j) = ELMU(j)
                  WRITE (22,99001) j , ELM(j)
               ELSE
                  ELM(j) = ELML(j)
                  WRITE (22,99001) j , ELM(j)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
99001 FORMAT (2X,'Warning - matrix element ',1I3,' reset to ',1F10.6)
      END
