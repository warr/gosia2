
C----------------------------------------------------------------------
C SUBROUTINE DECAY
C
C Called by: CEGRY, GOSIA
C Calls:     GKVAC
C
C Purpose: Calculate the gamma decay following excitation.
C
C Uses global variables:
C      DELLA  - products of matrix elements: e1^2, e2^2, e1*e2
C      DELTA  - \delta_\lambda: index 1 = electric^2, 2 = magnetic^2, 3 = cross term
C      GKP    - Gk * DELTA^2
C      IAXS   - axial symmetry flag
C      IBYP   - flag to indicate whether we calculate <\alpha_k>
C      IEXP   - experiment number
C      KLEC   - number of decays for each level
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LIFCT  - index of level for lifetimes
C      NMAX   - number of levels
C      NMAX1  - number of levels with decays
C      TAU    - lifetime in picoseconds
C      TIMEL  - lifetimes and their errors
C      VACDP  - G_k for each level
C      ZETA   - various coefficients
C
C Formal parameters:
C      Chisq  - chi squared
C      Nlift  - number of lifetimes
C      Chilo  - chi squared of logs

      SUBROUTINE DECAY(Chisq,Nlift,Chilo)
      IMPLICIT NONE
      REAL*8 bsum , Chilo , Chisq , df , el1 , emt , emt1 , gk , vcd
      INTEGER*4 i , ibra , idr , idrh , ifn , il , inx , inx1 , iu ,
     &          j , jlt , k , kl , kq , l , l1 , lc1 , lc2 , n1 , n2
      INTEGER*4 Nlift
      INCLUDE 'tra.inc'
      INCLUDE 'life1.inc'
      INCLUDE 'vac.inc'
      INCLUDE 'ccoup.inc'
      INCLUDE 'lev.inc'
      INCLUDE 'coex2.inc'
      INCLUDE 'comme.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'catlf.inc'
      INCLUDE 'lcdl.inc'
      DIMENSION gk(4)
      DATA emt1/0./

      idr = 1
      DO il = 1 , NMAX1 ! For each level with decays
         l = KSEQ(idr,3) ! Initial level of idr'th decay
         n1 = 28*(l-1)
         ibra = KLEC(l) ! Number of decays from level l
         bsum = 0.
         idrh = idr
         DO j = 1 , ibra ! For each decay from level l
            inx = KSEQ(idr,1) ! Index to matrix element of idr'th decay
            inx1 = KSEQ(idr,2) ! Index 2 of idr'th decay
            el1 = 0.
            IF ( inx.NE.0 ) el1 = ELM(inx)
            emt = el1*el1
            DELLA(idr,1) = emt
            IF ( inx1.NE.0 ) emt1 = ELM(inx1)*ELM(inx1)
            bsum = bsum + DELTA(idr,1)*emt
            IF ( inx1.NE.0 ) THEN
               DELLA(idr,3) = el1*ELM(inx1)
               DELLA(idr,2) = emt1
               bsum = bsum + DELTA(idr,2)*emt1
            ENDIF
            idr = idr + 1
         ENDDO ! Loop on j

         idr = idrh
         TAU(l) = 1./bsum
         CALL GKVAC(l) ! Evaluate G_k

         DO j = 1 , ibra ! For each decay from level l
            l1 = KSEQ(idr,4) ! Final energy of idr'th decay
            n2 = 28*(l1-1)
            inx1 = KSEQ(idr,2) ! Index 2 of idr'th decay
            DO i = 1 , 4
               gk(i) = GKP(i,idr,1)*DELLA(idr,1)
            ENDDO
            IF ( inx1.NE.0 ) THEN
               DO i = 1 , 4
                  gk(i) = gk(i) + GKP(i,idr,2)*DELLA(idr,2)
               ENDDO
            ENDIF
            DO i = 1 , 4
               vcd = 1.
               IF ( i.NE.1 ) vcd = VACDP(i-1,l)
               gk(i) = gk(i)*TAU(l)
               ifn = 2*i - 1
               iu = (i-1)*7
               IF ( IAXS(IEXP).EQ.0 ) ifn = 1
               DO kq = 1 , ifn
                  lc1 = n1 + iu + kq
                  lc2 = n2 + iu + kq
                  ZETA(lc2) = ZETA(lc2) + gk(i)*vcd*ZETA(lc1)
               ENDDO
            ENDDO
            idr = idr + 1
         ENDDO ! Loop on j
      ENDDO ! Loop on l

      IBYP = 1 ! Set flag to indicate we have calculated <\alpha_k>
      IF ( Nlift.NE.0 .AND. IEXP.EQ.1 ) THEN
         DO jlt = 1 , Nlift ! For each lifetime
            kl = LIFCT(jlt) ! Get level for this lifetime
            df = TAU(kl)-TIMEL(1,jlt) ! TIMEL(1,X) is lifetime
            IF ( df .LT. 0 ) THEN
              df = df / TIMEL(2,jlt) ! TIMEL(2,X) is lower limit
            Chilo = Chilo + (LOG(TAU(kl)/TIMEL(1,jlt))*TIMEL(1,jlt)
     &              /TIMEL(2,jlt))**2 ! Log chisqr
            ELSE
              df = df / TIMEL(3,jlt) ! TIMEL(3,X) is upper limit
            Chilo = Chilo + (LOG(TAU(kl)/TIMEL(1,jlt))*TIMEL(1,jlt)
     &              /TIMEL(3,jlt))**2 ! Log chisqr
            ENDIF
            Chisq = Chisq + df*df ! Chisqr
         ENDDO
      ENDIF

      DO l = 2 , NMAX ! For each level except the ground state
         IF ( KLEC(l).NE.0 ) THEN ! If there are decays from this level
            n1 = 28*(l-1)
            DO j = 1 , 4
               vcd = 1.
               IF ( j.NE.1 ) vcd = VACDP(j-1,l) ! G_k for each level
               ifn = 2*j - 1
               iu = (j-1)*7
               DO k = 1 , ifn
                  lc1 = n1 + iu + k
                  ZETA(lc1) = ZETA(lc1)*vcd
               ENDDO
            ENDDO
         ENDIF
      ENDDO ! Loop on levels
      END
