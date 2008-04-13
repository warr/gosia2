      SUBROUTINE APRAM(Iexp,Inc,Indx,Irld,Acca)
      IMPLICIT NONE
      REAL*8 Acca , accah , ELM , ELML , ELMU , QAPR , SA , uwa
      INTEGER*4 i1 , i56 , i7 , IAPR , IDIVE , Iexp , img , Inc , Indx , 
     &          IPATH , Irld , ISEX , itm , IVAR , j , jidim , jj , k , 
     &          ktoto , l
      INTEGER*4 l1 , l2 , l3 , LERF , LMAXE , m , MAGA , MAGEXC , 
     &          MEMAX , MEMX6
      COMPLEX*16 ARM
      COMMON /AZ    / ARM(600,7)
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /APRX  / LERF , IDIVE(50,2)
      LERF = 0
      accah = Acca
 100  i7 = 7
      itm = -1
      img = 3
      i1 = 1
      IF ( MAGA(Iexp).EQ.0 ) THEN
         i7 = 4
         i1 = 4
         img = 1
      ENDIF
      IF ( Inc.EQ.0 ) GOTO 300
      IF ( LERF.EQ.0 ) CALL NEWCAT(Iexp,jidim)
      IF ( LERF.EQ.0 ) CALL PODZIEL(3,Iexp)
      i56 = 5
      DO k = 1 , jidim
         ARM(k,2) = (0.,0.)
         ARM(k,5) = (0.,0.)
      ENDDO
      ARM(Irld+1,5) = (1.,0.)
 200  ktoto = 0
      LERF = 0
      l1 = IDIVE(Iexp,1)
      DO l3 = 1 , l1
         Acca = accah*l3/l1
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim)
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(1,Iexp)
            GOTO 100
         ENDIF
      ENDDO
      l2 = IDIVE(Iexp,2)
      DO l3 = 1 , l2
         Acca = accah + accah*l3/l2
         CALL POMNOZ(Acca,2,i56,ktoto,img,jidim)
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(2,Iexp)
            GOTO 100
         ENDIF
      ENDDO
      DO l = 1 , MEMX6
         DO m = i1 , i7
            QAPR(l,1,m) = -QAPR(l,1,m)
         ENDDO
      ENDDO
      DO l3 = 1 , l1
         Acca = accah*2. + accah*l3/l1
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim)
      ENDDO
      Acca = accah
      DO l = 1 , MEMX6
         DO m = i1 , i7
            QAPR(l,1,m) = -QAPR(l,1,m)
         ENDDO
      ENDDO
      IF ( Inc.NE.0 .OR. itm.NE.0 ) THEN
         IF ( Inc.EQ.0 ) THEN
            DO l = 1 , jidim
               ARM(l,6) = ARM(l,6) - ARM(l,7)
               ARM(l,6) = 50.*ARM(l,6)/ELM(Indx)
            ENDDO
            DO l = 1 , 2
               DO j = i1 , i7
                  QAPR(Indx,l,j) = QAPR(Indx,l,j)/.99
               ENDDO
            ENDDO
            DO jj = 2 , jidim
               ARM(jj-1,6) = ARM(jj,6)
            ENDDO
            GOTO 99999
         ELSE
            DO jj = 2 , jidim
               ARM(jj-1,5) = ARM(jj,5)
            ENDDO
            RETURN
         ENDIF
      ENDIF
 300  itm = itm + 1
      i56 = itm + 6
      DO k = 1 , jidim
         ARM(k,i56) = (0.,0.)
      ENDDO
      ARM(Irld+1,i56) = (1.,0.)
      uwa = -itm*.0298019802 + 1.01
      DO l = 1 , 2
         DO j = i1 , i7
            QAPR(Indx,l,j) = QAPR(Indx,l,j)*uwa
         ENDDO
      ENDDO
      DO j = 1 , jidim
         ARM(j,2) = (0.,0.)
      ENDDO
      GOTO 200
99999 END
