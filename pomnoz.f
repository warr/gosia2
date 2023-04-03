
C----------------------------------------------------------------------
C SUBROUTINE POMNOZ
C
C Called by: APRAM
C Calls:     TCABS
C
C Purpose: perform the expansion to calculate the approximate Coulomb
C          amplitudes
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      IAPR   - index of initial and final levels for each matrix element
C      INHB   - inhibit error flag setting (LERF)
C      IPATH  - index of substate in level with same m as substate Irld
C      ISEX   -
C      LERF   - error flag which is set here and used in APRAM
C      MEMX6  - number of matrix elements with E1...6 multipolarity
C      QAPR   - approximate Coulomb amplitudes
C
C Formal parameters:
C      Acca   - accuracy required
C      L
C      Iw
C      Img
C      Jidim
C      Ktoto  - number of iterations needed

      SUBROUTINE POMNOZ(Acca,L,Iw,Ktoto,Img,Jidim)
      IMPLICIT NONE
      REAL*8 Acca , QAPR , sig , TCABS , test , u
      INTEGER*4 IAPR , IDIVE , Img , INHB , IPATH , ISEX , IVAR , Iw ,
     &          Jidim , k , kk , Ktoto , L , LERF , LMAXE , m , MAGA ,
     &          MAGEXC , mc , mc1
      INTEGER*4 MEMAX , MEMX6 , mw , mw1
      COMPLEX*16 ARM , ci
      COMMON /INHI  / INHB
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(1500)
      COMMON /AZ    / ARM(600,7)
      COMMON /APRX  / LERF , IDIVE(50,2)
      DATA ci/(0.,-1.)/ ! -sqrt(-1)

      sig = 1.
      IF ( L.NE.2 ) sig = -1.
      DO kk = 1 , Jidim
         ARM(kk,1) = ARM(kk,Iw)
      ENDDO

      DO k = 1 , 100 ! Perform up to 100 iterations
         Ktoto = Ktoto + 1
         DO m = 1 , MEMX6 ! Matrix elements for E1...6
            mw1 = IAPR(m,1)
            mc1 = IAPR(m,2)
            IF ( IPATH(mw1).NE.0 .AND. IPATH(mc1).NE.0 ) THEN
               mw = IPATH(mw1) + 1
               mc = IPATH(mc1) + 1
               IF ( Ktoto.GE.ISEX(mc1) ) THEN
                  IF ( Img.EQ.1 ) THEN
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,4)*ARM(mw,1)
                  ELSE
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,4)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,4)*ARM(mw,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPR(m,L,2)*ARM(mc,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,2)*ARM(mw-1,1)
                     ARM(mw-1,2) = ARM(mw-1,2) + QAPR(m,L,1)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPR(m,L,1)
     &                             *ARM(mw-1,1)
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,3)*ARM(mc-1,1)
                     ARM(mc-1,2) = ARM(mc-1,2) + sig*QAPR(m,L,3)
     &                             *ARM(mw,1)
                     ARM(mw,2) = ARM(mw,2) + QAPR(m,L,5)*ARM(mc+1,1)
                     ARM(mc,2) = ARM(mc,2) + sig*QAPR(m,L,6)*ARM(mw+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPR(m,L,5)
     &                             *ARM(mw,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPR(m,L,6)*ARM(mc,1)
                     ARM(mw+1,2) = ARM(mw+1,2) + QAPR(m,L,7)*ARM(mc+1,1)
                     ARM(mc+1,2) = ARM(mc+1,2) + sig*QAPR(m,L,7)
     &                             *ARM(mw+1,1)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO

C        Calculate accuracy we have achieved
         test = 0.
         DO m = 1 , Jidim
            ARM(m,1) = ARM(m,2)/k
            ARM(m,2) = (0.,0.)
            IF ( L.NE.1 ) ARM(m,1) = ARM(m,1)*ci
            ARM(m,Iw) = ARM(m,Iw) + ARM(m,1)
            IF ( k.GT.5 ) THEN
               u = TCABS(ARM(m,Iw))
               test = test + u*u
            ENDIF
         ENDDO
C        Test to see if we have achieved required accuracy
         IF ( ABS(test-1.).LT.Acca ) GOTO 99999 ! Accuracy OK, so end
      ENDDO ! Iteration loop

      IF ( INHB.NE.1 ) LERF = 1 ! Accuracy not achieved, so set error flag
99999 END
