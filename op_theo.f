
C----------------------------------------------------------------------
C SUBROUTINE OP_THEO
C
C Called by: GOSIA
C
C Purpose: Calculate collective model matrix elements
C
C Uses global variables:
C      ELM    - matrix elements
C      JZB    - unit to read from
C      LEAD   - pair of levels involved in each matrix element
C      LP1    - maximum number of experiments (50)
C      MEMAX  - number of matrix elements
C      SPIN   - spin of level
C
C This function handles OP,THEO, which generates the matrix elements specified
C by the user and writes them into the file read by OP,REST

      SUBROUTINE OP_THEO
      IMPLICIT NONE
      INCLUDE 'clcom.inc'
      INCLUDE 'cexc.inc'
      INCLUDE 'coex.inc'
      INCLUDE 'comme.inc'
      INCLUDE 'mgn.inc'
      INCLUDE 'switch.inc'
      REAL*8 xlevb(50,2), bm(8,20,20,3), xm1, xm2, xm3
      REAL*8 xk1, xk2, xi1, xi2, ELMT
      INTEGER*4 irix, ibaf, ib, jb, lb, nbands, levl(50), ilevls
      INTEGER*4 jb1, jb2, jl, kb, kl, lamd, nnl, j
      INTEGER*4 nb1, nb2, nl, ind1, ind2, inva, bk

C     Initialise
      irix = 12
      REWIND (irix)
      ibaf = 1
      DO jb = 1 , LP1 ! LP1 = 50 (maximum number of experiments)
        DO lb = 1 , 2
          xlevb(jb,lb) = 0
        ENDDO
      ENDDO

C     Read the number of bands and initialise internal array bm
      READ (JZB,*) nbands ! Number of bands
      IF ( nbands.LE.0 ) ibaf = 0
      nbands = ABS(nbands)
      DO nl = 1 , 8
        DO jb = 1 , nbands
          DO jl = 1 , nbands
            DO kl = 1 , 3
              bm(nl,jb,jl,kl) = 0.
            ENDDO
          ENDDO
        ENDDO
      ENDDO

C     For each band, read in the band definition and list of levels in
C     that band
      DO jb = 1 , nbands
        READ (JZB,*) bk , ilevls ! K of band, number of levels in band
        READ (JZB,*) (levl(ib),ib=1,ilevls) ! Level list for band
        DO kb = 1 , ilevls
          inva = levl(kb)
          xlevb(inva,2) = bk
          xlevb(inva,1) = DBLE(jb)
        ENDDO
      ENDDO

C     For each multipolarity, read the band indices and the intrinsic
C     moments
      DO nl = 1 , 8
        READ (JZB,*) nnl ! Multipolarity
 126    IF ( nnl.LE.0 ) GOTO 130
        READ (JZB,*) jb1 , jb2 ! band indices
        IF ( jb1.NE.0 ) THEN
          READ (JZB,*) (bm(nnl,jb1,jb2,j),j=1,3) ! intrinsic moments
          DO j = 1 , 3
            bm(nnl,jb2,jb1,j) = bm(nnl,jb1,jb2,j)
          ENDDO
          GOTO 126
        ENDIF
      ENDDO

C     Now do the calculation loop. The work is done by the ELMT function
 130  DO kb = 1 , MEMAX
        IF ( ibaf.NE.0 ) THEN
          ind1 = LEAD(1,kb)
          ind2 = LEAD(2,kb)
          xi1 = SPIN(ind1)
          xi2 = SPIN(ind2)
          lamd = mlt(kb)
          nb1 = INT(xlevb(ind1,1)+.1)
          nb2 = INT(xlevb(ind2,1)+.1)
          xk1 = xlevb(ind1,2)
          xk2 = xlevb(ind2,2)
          xm1 = bm(lamd,nb1,nb2,1)
          xm2 = bm(lamd,nb1,nb2,2)
          xm3 = bm(lamd,nb1,nb2,3)
          ELM(kb) = ELMT(xi1,xi2,lamd,nb1,nb2,xk1,xk2,xm1,xm2,xm3)
          IF ( ABS(ELM(kb)).LT.1D-6 ) ELM(kb) = 1.D-6
          irix = 12
          IF ( IBPS.EQ.1 ) irix = 32 ! unit 12 for particle 1, 32 for particle 2 (only gosia2)
          WRITE (irix,*) ELM(kb)
        ENDIF
      ENDDO
      END
