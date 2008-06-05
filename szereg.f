 
C----------------------------------------------------------------------
C SUBROUTINE SZEREG
C
C Called by: READY
C
C Purpose: sort out the decay sequence
C
C Uses global variables:
C      DYEX   - error on experimental yield
C      IY     - index for yields
C      YEXP   - experimental yield
C
C Formal parameters:
C      Lst    - first yield in set
C      Ls     - last yield in set
C      L      - number of dataset
 
      SUBROUTINE SZEREG(Lst,Ls,L)
      IMPLICIT NONE
      REAL*8 dyh , yh
      INTEGER*4 ia , ib , ih , inx , k , L , Ls , lsp , Lst , lst1
      INCLUDE 'yexpt.inc'

      IF ( Lst.EQ.Ls ) RETURN
      lst1 = Lst
      lsp = Ls - 1

 100  ia = IY(lst1,L)
      IF ( ia.GT.1000 ) ia = ia/1000

      inx = lst1
      DO k = lst1 , lsp
         ib = IY(k+1,L)
         IF ( ib.GT.1000 ) ib = ib/1000
         ia = MIN(ia,ib)
         IF ( ia.EQ.ib ) inx = k + 1
      ENDDO

C     Swap them
      IF ( inx.NE.lst1 ) THEN
         ih = IY(lst1,L)
         IY(lst1,L) = IY(inx,L)
         IY(inx,L) = ih
         yh = YEXP(L,lst1)
         dyh = DYEX(L,lst1)
         YEXP(L,lst1) = YEXP(L,inx)
         DYEX(L,lst1) = DYEX(L,inx)
         YEXP(L,inx) = yh
         DYEX(L,inx) = dyh
      ENDIF

      lst1 = lst1 + 1
      IF ( lst1.GT.lsp ) RETURN
      GOTO 100
      END
