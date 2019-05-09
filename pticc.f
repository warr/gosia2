
C----------------------------------------------------------------------
C SUBROUTINE PTICC
C
C Called by: GOSIA
C Calls:     CONV
C
C Purpose: print the conversion coefficients
C
C Uses global variables:
C      EN     - energy of level
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      MULTI  - number of matrix elements having a given multipolarity
C      SPIN   - spin of level
C
C Formal parameters:
C      Idr    - number of decays

      SUBROUTINE PTICC(Idr)
      IMPLICIT NONE
      REAL*8 cone1 , cone2 , conm1 , CONV , enet
      INTEGER*4 Idr , iinx , l , nf , ni
      INCLUDE 'clcom.inc'
      INCLUDE 'coex.inc'
      INCLUDE 'lev.inc'

      WRITE (22,99001)
99001 FORMAT (1X//20X,'CALCULATED INTERNAL CONVERSION ',
     &        'COEFFICIENTS FOR E1,E2 AND M1'//5X,'NI',5X,'NF',7X,'II',
     &        8X,'IF',9X,'ENERGY(MEV)',6X,'ICC(E1)',8X,'ICC(E2)',8X,
     &        'ICC(M1)')
      DO l = 1 , Idr
         iinx = KSEQ(l,1) ! Index of l'th decay
         ni = KSEQ(l,3) ! Initial level of l'th decay
         nf = KSEQ(l,4) ! Final level of l'th decay
         enet = EN(ni) - EN(nf)
         cone2 = CONV(enet,2)
         IF ( ABS(SPIN(ni)-SPIN(nf)).GT.2. ) cone2 = 0.
         conm1 = 0.
         cone1 = 0.
         IF ( iinx.LE.MULTI(1) ) cone1 = CONV(enet,1)
         IF ( ABS(SPIN(ni)-SPIN(nf)).LT.2. ) conm1 = CONV(enet,4)
         WRITE (22,99002) ni , nf , SPIN(ni) , SPIN(nf) , enet , cone1 ,
     &                    cone2 , conm1
99002    FORMAT (4X,I3,4X,I3,7X,F4.1,6X,F4.1,9X,F6.4,8X,E9.4,6X,E9.4,6X,
     &           E9.4)
      ENDDO
      END
