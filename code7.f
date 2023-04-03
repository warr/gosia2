
C----------------------------------------------------------------------
C SUBROUTINE CODE7
C
C Called by: LSLOOP
C
C Purpose:
C
C Uses global variables:
C      IAPR   - index of initial and final levels for each matrix element
C      IPATH  - index of substate in level with same m as substate Irld
C
C Formal parameters:
C      Ir     - index of initial substate
C      Is     - index of final substate
C      N      - index of initial level
C      Mt     - index of final level
C      Inqa   - result of operation
C      Indx   - Index of matrix element

      SUBROUTINE CODE7(Ir,Is,N,Mt,Inqa,Indx)
      IMPLICIT NONE
      INTEGER*4 IAPR , idm , idn , Indx , Inqa , IPATH , Ir , Is ,
     &          ISEX , ism , MAGA , Mt , N
      REAL*8 QAPR
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /APRCAT/ QAPR(1500,2,7) , IAPR(1500,2) , ISEX(75)

      IAPR(Indx,1) = N  ! Index of initial level
      IAPR(Indx,2) = Mt ! Index of final level
      IF ( IPATH(N).EQ.0 .OR. IPATH(Mt).EQ.0 ) THEN
         Inqa = -1
         GOTO 99999
      ELSE
         idn = Ir - IPATH(N)
         idm = Is - IPATH(Mt)
         ism = idn + idm + 3
         IF ( ism.EQ.2 ) THEN
            Inqa = 2
            IF ( idn.GT.idm ) Inqa = 3
            RETURN
         ELSEIF ( ism.EQ.3 ) THEN
            Inqa = 4
            RETURN
         ELSEIF ( ism.EQ.4 ) THEN
            Inqa = 5
            IF ( idn.GT.idm ) Inqa = 6
            RETURN
         ELSEIF ( ism.NE.5 ) THEN
            Inqa = 1
            RETURN
         ENDIF
      ENDIF
      Inqa = 7
99999 END
