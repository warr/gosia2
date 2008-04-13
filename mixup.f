      SUBROUTINE MIXUP
      IMPLICIT NONE
      REAL*8 ELM , ELML , ELMU , RNDM , SA , SE
      INTEGER*4 IVAR , k , k1 , LMAXE , MAGEXC , MEMAX , MEMX6
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /XRA   / SE
      DO k = 1 , MEMAX
         IF ( IVAR(k).NE.0 .AND. IVAR(k).LE.999 ) ELM(k) = ELML(k)
     &        + RNDM(SE)*(ELMU(k)-ELML(k))
      ENDDO
      DO k = 1 , MEMAX
         IF ( IVAR(k).GE.999 ) THEN
            k1 = IVAR(k) - 1000
            IF ( ABS(ELMU(k1)).LT.1.E-9 ) THEN
               ELM(k) = 0.
            ELSE
               ELM(k) = ELM(k1)*SA(k)
            ENDIF
         ENDIF
      ENDDO
      END
