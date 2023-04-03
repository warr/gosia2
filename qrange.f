
C----------------------------------------------------------------------
C SUBROUTINE QRANGE
C
C Called by: SNAKE
C
C Purpose: determine the range for which we will need Qe and Qm values.
C
C Uses global variables:
C      IRA    - range to integrate over omega
C      MAXLA  - multipolarity to calculate
C      MULTI  - number of matrix elements having given multipolarity
C
C Formal parameters:
C      Icnt   - index of omega to calculate
C      Nlm    - returns the number of l,m values
C      Lloc   -
C      Ibm    -
C      Icm    -
C      Idm    -
C      Irl    -

      SUBROUTINE QRANGE(Icnt,Nlm,Lloc,Ibm,Icm,Idm,Irl)
      IMPLICIT NONE
      INTEGER*4 Ibm , Icm , Icnt , Idm , IRA , Irl , is , k , ke , km ,
     &          l , LAMDA , LAMMAX , ld , LDNUM , LEAD , Lloc , ls ,
     &          MAXLA , MULTI
      INTEGER*4 nlend , Nlm
      DIMENSION Lloc(8) , Irl(8)
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,75) , LAMMAX ,
     &                MULTI(8)

      IF ( Icnt.EQ.1 ) THEN
         Nlm = 0
         DO l = 1 , 8
            Lloc(l) = 0
            Irl(l) = 0
         ENDDO
         DO k = 1 , 6
            ke = 7 - k
            km = 13 - k
            IF ( km.LE.8 ) THEN
               IF ( MULTI(km).NE.0 ) THEN
                  Nlm = Nlm + 1
                  Lloc(Nlm) = km
                  Irl(Nlm) = IRA(km)
               ENDIF
            ENDIF
            IF ( MULTI(ke).NE.0 ) THEN
               Nlm = Nlm + 1
               Lloc(Nlm) = ke
               Irl(Nlm) = IRA(ke)
            ENDIF
         ENDDO
         nlend = INT((DBLE(Nlm)+1.1)/2.)
         DO k = 1 , nlend
            ke = Nlm - k + 1
            ls = Lloc(ke)
            is = Irl(ke)
            Lloc(ke) = Lloc(k)
            Irl(ke) = Irl(k)
            Lloc(k) = ls
            Irl(k) = is
         ENDDO
         l = 0
         DO k = 1 , 6
            IF ( MULTI(k).NE.0 ) l = k
         ENDDO
         Icm = MIN(4,l)
         Ibm = 2*l
         Idm = l
         l = 0
         DO k = 7 , 8
            ke = k - 6
            IF ( MULTI(k).NE.0 ) l = ke
         ENDDO
         Ibm = MAX(Ibm,2*l)
         Idm = MAX(Idm,l)
         IF ( Icm.EQ.1 .AND. l.GT.1 ) Icm = 2
         MAXLA = Lloc(1)
         RETURN
      ELSE
         IF ( Irl(Nlm).GE.Icnt ) RETURN
         ld = Lloc(Nlm)
         Lloc(Nlm) = 0
         Nlm = Nlm - 1
         IF ( Nlm.EQ.0 ) RETURN
         IF ( ld.GT.6 ) RETURN
         l = Lloc(Nlm)
         IF ( l.GT.6 ) l = l - 6
         Icm = MIN(2,l)
         Ibm = 2*l
         Idm = l
      ENDIF
      END
