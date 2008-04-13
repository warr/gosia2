      SUBROUTINE CHMEM(Nw,Chi,Chilo)
      IMPLICIT NONE
      REAL*8 Chi , Chilo , di , EAMX , ELM , ELML , ELMU , SA
      INTEGER*4 ia , IAMX , IAMY , ib , NAMX , Nw
      COMMON /ME2D  / EAMX(100,2) , NAMX , IAMX(100) , IAMY(100,2)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      IF ( NAMX.EQ.0 ) RETURN
      Nw = Nw + NAMX
      DO ia = 1 , NAMX
         ib = IAMX(ia)
         IF ( IAMY(ia,1).NE.IAMY(ia,2) ) THEN
            di = (ELM(ib)-EAMX(ia,1))/EAMX(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &              /EAMX(ia,2))**2
            Chi = Chi + di*di
         ELSE
            di = (ELM(ib)-EAMX(ia,1))/EAMX(ia,2)
            Chilo = Chilo + 
     &              (LOG(ABS(ELM(ib)/EAMX(ia,1)))*ABS(EAMX(ia,1))
     &              /EAMX(ia,2))**2
            Chi = Chi + di*di
         ENDIF
      ENDDO
      END
