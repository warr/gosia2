      SUBROUTINE SZEREG(Lst,Ls,L)
      IMPLICIT NONE
      REAL*8 CORF , DYEX , dyh , UPL , YEXP , yh , YNRM
      INTEGER*4 ia , ib , IDRN , ih , ILE , inx , IY , k , L , Ls , 
     &          lsp , Lst , lst1 , NYLDE
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
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
