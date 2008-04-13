      SUBROUTINE PTICC(Idr)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , cone1 , cone2 , conm1 , CONV , DIPOL , EN , 
     &       enet , SPIN , TAU , ZPOL
      INTEGER*4 Idr , iinx , ISO , KSEQ , l , LAMDA , LAMMAX , LDNUM , 
     &          LEAD , MULTI , nf , ni
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      WRITE (22,99001)
99001 FORMAT (1X//20X,'CALCULATED INTERNAL CONVERSION ',
     &        'COEFFICIENTS FOR E1,E2 AND M1'//5X,'NI',5X,'NF',7X,'II',
     &        8X,'IF',9X,'ENERGY(MEV)',6X,'ICC(E1)',8X,'ICC(E2)',8X,
     &        'ICC(M1)')
      DO l = 1 , Idr
         iinx = KSEQ(l,1)
         ni = KSEQ(l,3)
         nf = KSEQ(l,4)
         enet = EN(ni) - EN(nf)
         cone2 = CONV(enet,2)
         IF ( ABS(SPIN(ni)-SPIN(nf)).GT.2. ) cone2 = 0.
         conm1 = 0.
         cone1 = 0.
         IF ( iinx.LE.MULTI(1) ) cone1 = CONV(enet,1)
         IF ( ABS(SPIN(ni)-SPIN(nf)).LT.2. ) conm1 = CONV(enet,4)
         WRITE (22,99002) ni , nf , SPIN(ni) , SPIN(nf) , enet , cone1 , 
     &                    cone2 , conm1
99002    FORMAT (5X,I2,5X,I2,7X,F4.1,6X,F4.1,9X,F6.4,8X,E9.4,6X,E9.4,6X,
     &           E9.4)
      ENDDO
      END
