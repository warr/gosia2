 
C----------------------------------------------------------------------
C SUBROUTINE SNAKE
C
C Called by: FTBM, GOSIA
C Calls:     QE, QM, QRANGE
C
C Purpose: evaluate and store the dimensionless collision functions Qe and Qm.
C
C Uses global variables:
C      CH     - table of cosh values
C      EPS    - epsilon
C      EROOT  - sqrt(epsilon^2 -1)
C      LOCQ   - location of collision function in ZETA array
C      LP7    - start of collision functions in ZETA (45100)
C      SH     - table of sinh values
C      ZETA   - various coefficients (here the collision functions)
C
C Formal parameters:
C      Nexp   - experiment number
C      Zpol   -
C
C The function QE is used to calculate Qe and QM to calculate Qm, but first
C we call QRANGE to determine the range over which we need to calculate them.
C
C The results are stored in the ZETA array, but not starting from the
C beginning, which is where zeta itself is written, but from ZETA(LP7).
C
C LOCQ (in ALLC) is used as an index to these values.
C
C EROOT is set in CMLAB to \sqrt(\epsilon^2 - 1).
C
C Note that when we call QE and QM that lmda = 1...6 for E1...6 and 7,8 for
C M1, M2.
 
      SUBROUTINE SNAKE(Nexp,Zpol)
      IMPLICIT NONE
      REAL*8 b10 , b12 , b2 , b4 , b6 , b8 , c , c2 , c4 , c6 , CH , 
     &       chi , cq , d , d2 , d3 , d4 , d5 , d6 , EPS
      REAL*8 EROOT , ert , FIEX , pol , SH , shi , ZETA , Zpol
      INTEGER*4 IAXS , ibm , icm , icnt , idm , IEXP , irl , j , k , 
     &          lloc , lmd , lmda , LOCQ , LP1 , LP10 , LP11 , LP12 , 
     &          LP13 , LP14 , LP2
      INTEGER*4 LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , LZETA , mimx , 
     &          Nexp , nind , nlm
      DIMENSION lloc(8) , cq(7) , irl(8)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /HIPER / SH(365) , CH(365)
      
      icnt = 0
 100  icnt = icnt + 1

C     Calculate range over which we will want Qe and Qm
      CALL QRANGE(icnt,nlm,lloc,ibm,icm,idm,irl)
      IF ( nlm.EQ.0 ) RETURN

C     Calculate some parameters, which we will pass to QE or QM
      chi = CH(icnt) ! \cosh(\omega)
      shi = SH(icnt) ! \sinh(\omega)
      b2 = EPS(Nexp)*chi + 1.
      pol = 1. - Zpol/b2
      b2 = b2*b2 ! b^2 = (\epsilon \cosh(\omega) + 1)^2
      IF ( ibm.NE.2 ) THEN
         b4 = b2*b2
         IF ( ibm.NE.4 ) THEN
            b6 = b4*b2
            IF ( ibm.NE.6 ) THEN
               b8 = b4*b4
               IF ( ibm.NE.8 ) THEN
                  b10 = b6*b4
                  IF ( ibm.NE.10 ) b12 = b6*b6
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      IF ( icm.NE.0 ) THEN
         c = chi + EPS(Nexp) ! c = \cosh(\omega + \epsilon)
         IF ( icm.NE.1 ) THEN
            c2 = c*c
            IF ( icm.NE.2 ) THEN
               c4 = c2*c2
               IF ( icm.NE.4 ) c6 = c2*c4
            ENDIF
         ENDIF
      ENDIF
      IF ( idm.NE.0 ) THEN
         d = EROOT(Nexp)*shi ! d = \sinh(\omega) * \sqrt(epsilon^2 - 1)
         IF ( idm.NE.1 ) THEN
            d2 = d*d
            IF ( idm.NE.2 ) THEN
               d3 = d*d2
               IF ( idm.NE.3 ) THEN
                  d4 = d2*d2
                  IF ( idm.NE.4 ) THEN
                     d5 = d3*d2
                     IF ( idm.NE.5 ) d6 = d3*d3
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      DO j = 1 , nlm
         lmda = lloc(j)
         IF ( lmda.GT.6 ) THEN
            lmd = lmda
            lmda = lmda - 6
            ert = EROOT(Nexp)
            CALL QM(c,d,b2,b4,ert,lmda,cq)
            mimx = lmda
            DO k = 1 , mimx
               nind = LOCQ(lmd,k) + icnt
               ZETA(nind+LP7) = cq(k) ! These are the collision functions
            ENDDO
         ELSE
            CALL QE(c,d,b2,c2,d2,b4,b6,d3,b8,c4,d4,b10,d5,b12,d6,lmda,
     &              pol,cq)
            mimx = lmda + 1
            DO k = 1 , mimx
               nind = LOCQ(lmda,k) + icnt
               ZETA(nind+LP7) = cq(k) ! These are the collision functions
            ENDDO
         ENDIF
      ENDDO
      GOTO 100
      END