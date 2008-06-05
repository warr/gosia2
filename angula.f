 
C----------------------------------------------------------------------
C SUBROUTINE ANGULA
C
C Called by: GOSIA, CEGRY
C Calls:     FIINT, FIINT1, RECOIL, YLM, YLM1
C
C Purpose: calculate angular distribution of emitted gamma rays
C
C Uses global variables:
C      BETAR  - recoil beta
C      DELLA  - products of matrix elements: e1^2, e2^2, e1*e2
C      ENDEC  - energy difference for each matrix element
C      ENZ    - something to do with the absorption
C      FP     - F coefficient * DELTA^2
C      IAXS   - axial symmetry flag
C      IEXP   - experiment number
C      ITMA   - identify detectors according to OP,GDET
C      ITTE   - thick target experiment flag
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      TAU    - lifetime in picoseconds
C      ZETA   - various coefficients
C
C Formal parameters:
C      Ygn    - Gamma-ray yield
C      Idr    - number of decays
C      Iful   - flag to select full basis or not
C      Fi0    - phi_0
C      Fi1    - phi_1
C      Trec   - Theta of recoiling nucleus
C      Gth    - Theta of gamma
C      Figl   - Phi of gamma
C      Ngl    - detector number
C      Op2    - The part after the OP, for the option we are processing (gosia2)
      
      SUBROUTINE ANGULA(Ygn,Idr,Iful,Fi0,Fi1,Trec,Gth,Figl,Ngl,Op2)
      IMPLICIT NONE
      REAL*8 alab , arg , at , attl , bt , f , Fi0 , fi01 , 
     &       Fi1 , fi11
      REAL*8 Figl , Gth , qv , sm , Trec , 
     &       Ygn , ylmr
      INTEGER*4 Idr , ifn , Iful , ig , il , inat , inx1 , 
     &          ipd , is , ITTE , iu , ixs , j , ji , jj , jm , k
      INTEGER*4 kq , l , lf , lf1 , mind , Ngl , nlv
      REAL*8 dsig, ttx
      CHARACTER*4 Op2 ! For gosia2
      DIMENSION f(4) , ylmr(9,9) , at(28) , alab(9,9) , attl(9,9) , 
     &          Ygn(*)
      INCLUDE 'ccoup.inc'
      INCLUDE 'tra.inc'
      INCLUDE 'lev.inc'
      INCLUDE 'ccc.inc'
      INCLUDE 'kin.inc'
      INCLUDE 'lcdl.inc'
      INCLUDE 'catlf.inc'
      INCLUDE 'brec.inc'
      COMMON /THTAR / ITTE(50)
      INCLUDE 'tcm.inc' ! For gosia2
      INCLUDE 'cx.inc' ! For gosia2
      
      DO l = 1 , Idr ! For each decay

         nlv = KSEQ(l,3) ! Level number of l'th decay
         il = (nlv-1)*28
         inx1 = KSEQ(l,2) ! Index of l'th decay

         DO j = 1 , 4
            f(j) = FP(j,l,1)*DELLA(l,1)
         ENDDO

         IF ( inx1.NE.0 ) THEN
            DO j = 1 , 4
               f(j) = f(j) + 2.*FP(j,l,3)*DELLA(l,3) + FP(j,l,2)
     &                *DELLA(l,2)
            ENDDO
         ENDIF

         DO j = 1 , 4
            f(j) = f(j)*TAU(nlv)
            iu = (j-1)*7
            ifn = 2*j - 1
            IF ( IAXS(IEXP).EQ.0 ) ifn = 1
            DO kq = 1 , ifn
               is = iu + kq
               ig = is + il
               at(is) = ZETA(ig)*f(j)
            ENDDO
         ENDDO

         IF ( Iful.EQ.1 ) THEN
            DO j = 1 , 9
               DO k = 1 , 9
                  alab(j,k) = 0.
                  attl(j,k) = 0.
               ENDDO
            ENDDO
            DO j = 1 , 4
               lf = 2*j - 1
               lf1 = lf
               IF ( IAXS(IEXP).EQ.0 ) lf1 = 1
               DO k = 1 , lf1
                  inat = (j-1)*7 + k
                  alab(lf,k) = at(inat)
               ENDDO
            ENDDO
            bt = BETAR(IEXP) ! Get beta
            IF ( ITTE(IEXP).NE.1 ) CALL RECOIL(alab,attl,bt,Trec)
            IF ( l.EQ.1 ) CALL YLM1(Gth,ylmr)
            ixs = IAXS(IEXP) ! Get axial symmetry flag
            fi01 = Fi0 - Figl ! Get lower phi limit
            fi11 = Fi1 - Figl ! Get upper phi limit
            CALL FIINT1(fi01,fi11,alab,ixs) ! Integrate over phi in lab frame
            Ygn(l) = alab(1,1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
            DO j = 2 , 9
               sm = ylmr(j,1)*alab(j,1)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  DO k = 2 , j
                     sm = sm + 2.*ylmr(j,k)*alab(j,k)
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl) ! Detector ID
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,j-1)*Q(2,ipd,j-1)+Q(1,ipd,j-1)*arg)
     &              /(Q(2,ipd,j-1)+arg)
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ELSE
            ixs = IAXS(IEXP) ! Get axial symmetry flag
            fi01 = Fi0 - Figl ! Get lower phi limit
            fi11 = Fi1 - Figl ! Get upper phi limit
            CALL FIINT(fi01,fi11,at,ixs) ! Integrate over phi in recoiling nucleus frame, result in at
            IF ( l.EQ.1 ) CALL YLM(Gth,ylmr)
            Ygn(l) = at(1)*.0795774715 ! 0.0795774715 = 1 / (4 pi)
            DO jj = 1 , 3
               ji = jj*7 + 1
               sm = ylmr(jj,1)*at(ji)
               IF ( IAXS(IEXP).NE.0 ) THEN
                  mind = 2*jj + 1
                  DO jm = 2 , mind
                     ji = ji + 1
                     sm = ylmr(jj,jm)*at(ji)*2. + sm
                  ENDDO
               ENDIF
               ipd = ITMA(IEXP,Ngl) ! Detector ID
               arg = (ENDEC(l)-ENZ(ipd))**2
               qv = (Q(3,ipd,2*jj)*Q(2,ipd,2*jj)+Q(1,ipd,2*jj)*arg)
     &              /(Q(2,ipd,2*jj)+arg) ! solid angle attenuation coefficients
               Ygn(l) = Ygn(l) + sm*qv
            ENDDO
         ENDIF
      ENDDO ! Loop over decays

C     Added for gosia2
      IF ( Op2.EQ.'INTG' ) RETURN
      dsig = DSIGS(IEXP)
      ttx = TLBDG(IEXP)/57.2957795 ! Theta in lab frame in radians
      DO j = 1 , Idr  ! For each decay
         Ygn(j) = Ygn(j)*dsig*SIN(ttx)
      ENDDO ! Loop on decays Idr
C     End of addition for gosia2
      END
