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
      CALL QRANGE(icnt,nlm,lloc,ibm,icm,idm,irl)
      IF ( nlm.EQ.0 ) RETURN
      chi = CH(icnt)
      shi = SH(icnt)
      b2 = EPS(Nexp)*chi + 1.
      pol = 1. - Zpol/b2
      b2 = b2*b2
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
         c = chi + EPS(Nexp)
         IF ( icm.NE.1 ) THEN
            c2 = c*c
            IF ( icm.NE.2 ) THEN
               c4 = c2*c2
               IF ( icm.NE.4 ) c6 = c2*c4
            ENDIF
         ENDIF
      ENDIF
      IF ( idm.NE.0 ) THEN
         d = EROOT(Nexp)*shi
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
               ZETA(nind+LP7) = cq(k)
            ENDDO
         ELSE
            CALL QE(c,d,b2,c2,d2,b4,b6,d3,b8,c4,d4,b10,d5,b12,d6,lmda,
     &              pol,cq)
            mimx = lmda + 1
            DO k = 1 , mimx
               nind = LOCQ(lmda,k) + icnt
               ZETA(nind+LP7) = cq(k)
            ENDDO
         ENDIF
      ENDDO
      GOTO 100
      END
