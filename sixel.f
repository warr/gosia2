 
C----------------------------------------------------------------------
C SUBROUTINE SIXEL
C
C Called by: CEGRY
C
C Purpose:
C
C Uses global variables:
C      ARM    - excitation amplitudes of substates.
C      DEV    -
C      IEXP   - experiment number
C      ITS    - create tape 18 file (OP,CONT switch SEL,)
C      KVAR   -
C
C Formal parameters:
C      Rik    - DEV + YEXP
C      Rv     - difference between experimental and calculated yields
C      Em     - matrix element
C      Jk     -
C      Kk     -
C      Indx   - index of matrix element
C      Lu     -
 
      SUBROUTINE SIXEL(Rik,Rv,Em,Jk,Kk,Indx,Lu)
      IMPLICIT NONE
      REAL*8 a1 , al , al1 , c1 , c2 , Em , Rik , rn , Rv , rx
      INTEGER*4 Indx , j , j1 , Jk , Kk , kk6 , l , l1 , Lu
      INCLUDE 'az.inc'
      INCLUDE 'odch.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'trb.inc'
      INCLUDE 'sel.inc'
      
      kk6 = Kk + 5
      rn = DEV(Lu)
      al = (Rv-rn)*20./Rik
      IF ( ITS.EQ.1 .AND. KVAR(Indx).NE.0 ) WRITE (18,*) Lu , Indx , 
     &     IEXP , al/Em
      al1 = ABS(al)
      IF ( ITS.EQ.2 ) WRITE (18,*) Lu , Indx , IEXP , al1
      IF ( al1.LE.ABS(IMAG(ARM(kk6,Jk))) ) RETURN

      DO j = Kk , kk6
         a1 = ABS(IMAG(ARM(j,Jk)))
         IF ( al1.GT.a1 ) THEN
            j1 = j + 1
            DO l = j1 , kk6
               l1 = kk6 + j1 - l
               c1 = DBLE(ARM(l1-1,Jk))
               c2 = IMAG(ARM(l1-1,Jk))
               ARM(l1,Jk) = DCMPLX(c1,c2)
            ENDDO
            rx = DBLE(Indx)
            ARM(j,Jk) = DCMPLX(rx,al)
            GOTO 99999
         ENDIF
      ENDDO

99999 END
