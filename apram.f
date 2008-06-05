 
C----------------------------------------------------------------------
C SUBROUTINE APRAM
C
C Called by: FTBM
C Calls:     NEWCAT, PODZIEL, POMNOZ
C
C Purpose: calculate approximate value of the Coulomb excitation amplitudes.
C
C Uses global parameters:
C      ARM    - excitation amplitudes of substates.
C      ELM    - matrix elements
C      IDIVE  - number of subdivisions
C      LERF   - error flag for expansion in POMNOZ
C      MAGA   - number of magnetic substates in approximate calculation
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      QAPR   - approximate Coulomb amplitudes
C
C Formal parameters:
C      Iexp   - experiment number
C      Inc    - flag: first time we call after LOAD, Inc=0, afterwards Inc=1
C      Indx   - index of matrix element
C      Irld   - index into ARM array
C      Acca   - accuracy required

      SUBROUTINE APRAM(Iexp,Inc,Indx,Irld,Acca)
      IMPLICIT NONE
      REAL*8 Acca , accah , ELM , ELML , ELMU , SA , uwa
      INTEGER*4 i1 , i56 , i7 , IDIVE , Iexp , img , Inc , Indx , 
     &          Irld , itm , IVAR , j , jidim , jj , k , 
     &          ktoto , l
      INTEGER*4 l1 , l2 , l3 , LERF , LMAXE , m , MAGEXC , 
     &          MEMAX , MEMX6
      INCLUDE 'az.inc'
      INCLUDE 'aprcat.inc'
      INCLUDE 'pth.inc'
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      COMMON /COMME / ELM(1500) , ELMU(1500) , ELML(1500) , SA(1500)
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
      IF ( LERF.EQ.0 ) CALL PODZIEL(3,Iexp) ! Subdivide
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
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim) ! Expansion for L=1
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(1,Iexp) ! Subdivide
            GOTO 100
         ENDIF
      ENDDO

      l2 = IDIVE(Iexp,2)
      DO l3 = 1 , l2
         Acca = accah + accah*l3/l2
         CALL POMNOZ(Acca,2,i56,ktoto,img,jidim) ! Expansion for L=2
         IF ( LERF.NE.0 ) THEN
            CALL PODZIEL(2,Iexp) ! Subdivide
            GOTO 100
         ENDIF
      ENDDO

      DO l = 1 , MEMX6 ! Matrix elements for E1...6
         DO m = i1 , i7
            QAPR(l,1,m) = -QAPR(l,1,m)
         ENDDO
      ENDDO

      DO l3 = 1 , l1
         Acca = accah*2. + accah*l3/l1
         CALL POMNOZ(Acca,1,i56,ktoto,img,jidim) ! Expansion for L=1
      ENDDO

      Acca = accah
      DO l = 1 , MEMX6 ! Matrix elements for E1...6
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

C     Initialise (Inc = 0)
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
