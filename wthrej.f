
C----------------------------------------------------------------------
C FUNCTION WTHREJ
C
C Called by: ELMT, F, GOSIA, LSLOOP, TENB
C
C Purpose: evaluates a Wigner 3-j symbol.
C
C Uses global variables:
C      IP     - table of prime numbers
C      KF     - sum of factors of primes
C      PILOG  - table of natural logs of primes
C
C Formal parameters:
C      J1     - twice the value of J1
C      J2     - twice the value of J2
C      J3     - twice the value of J3
C      M1     - twice the value of M1
C      M2     - twice the value of M2
C      M3     - twice the value of M3
C
C Return value:
C      The value of the 3-j symbol
C
C Note that the values of the parameters are doubled, so that this function
C can handle half-integers. In other words if you want to evaluate
C \threej(J1 J2 J3 M1 M2 M3) you need to use call the function as:
C WTHREJ(2 * J1, 2 * J2, 2 * J3, 2 * M1, 2 * M2, 2 * M3).

      REAL*8 FUNCTION WTHREJ(J1,J2,J3,M1,M2,M3)
      IMPLICIT NONE
      INTEGER*4 IP , IPI , iz , iza , izb , izc , izd , ize , izexp ,
     &          izf , izmax , izmin , J1 , J2 , J3 , jabc , jabm ,
     &          jbma , jj1 , jj2
      INTEGER*4 jj3 , jjha , jjhb , jjhc , jjhd , jlp , jma , jmax ,
     &          jmb , jmc , jmd , jme , jmf , jta , jtb , jtc , jvo ,
     &          jvora , KF , M1
      INTEGER*4 M2 , M3 , mm1 , mm2 , mm3 , n , nmax
      REAL*8 PILOG , qsumlo , sumlo , vorz , wthrep , zuthre
      DIMENSION jvora(26)
      COMMON /FAKUL / IP(26) , IPI(26) , KF(101,26) , PILOG(26)

      wthrep = 0.E+00
      jjha = (J1+J2-J3)/2 + 1
      jjhb = (J1-J2+J3)/2 + 1
      jjhc = (-J1+J2+J3)/2 + 1

      IF ( (jjha.LT.1) .OR. (jjhb.LT.1) .OR. (jjhc.LT.1) .OR.
     &     ((M1+M2+M3).NE.0) ) THEN
         WTHREJ = wthrep
         GOTO 99999
      ENDIF

      jjhd = (J1+J2+J3+4)/2
      jmax = MAX(J1,J2,J3)
      IF ( jmax.NE.J1 ) THEN
         IF ( jmax.EQ.J2 ) THEN
            jj1 = J3
            jj2 = J1
            jj3 = J2
            mm1 = M3
            mm2 = M1
            mm3 = M2
            GOTO 100
         ELSEIF ( jmax.EQ.J3 ) THEN
            jj1 = J1
            jj2 = J2
            jj3 = J3
            mm1 = M1
            mm2 = M2
            mm3 = M3
            GOTO 100
         ENDIF
      ENDIF
      jj1 = J2
      jj2 = J3
      jj3 = J1
      mm1 = M2
      mm2 = M3
      mm3 = M1

 100  jma = (jj1+mm1)/2
      jmb = (jj1-mm1)/2
      jmc = (jj2+mm2)/2
      jmd = (jj2-mm2)/2
      jme = (jj3+mm3)/2
      jmf = (jj3-mm3)/2
      jabc = (jj1+jj2-jj3)/2
      jabm = (jj2-jj3-mm1)/2
      jbma = (jj1+mm2-jj3)/2
      izmin = MAX(jabm,jbma,0)
      izmax = MIN(jabc,jmb,jmc)
      nmax = MAX(jjhd,izmax+1)
      DO n = 1 , 26
         IF ( IP(n).GE.nmax ) GOTO 200
      ENDDO
      WTHREJ = wthrep
      GOTO 99999

 200  DO jlp = 1 , n
         jta = KF(jjha,jlp) + KF(jjhb,jlp) + KF(jjhc,jlp) - KF(jjhd,jlp)
         jtb = KF(jma+1,jlp) + KF(jmb+1,jlp) + KF(jmc+1,jlp)
         jtc = KF(jmd+1,jlp) + KF(jme+1,jlp) + KF(jmf+1,jlp)
         jvora(jlp) = jta + jtb + jtc
      ENDDO

      vorz = -1.E+00
      IF ( 2*(izmin/2).EQ.izmin ) vorz = +1.E+00
      IF ( izmin.LE.izmax ) THEN
         DO iz = izmin , izmax
            qsumlo = 0.E+00
            iza = iz + 1
            izb = jabc + 1 - iz
            izc = jmb + 1 - iz
            izd = jmc + 1 - iz
            ize = iz - jabm + 1
            izf = iz - jbma + 1
            DO jlp = 1 , n
               izexp = jvora(jlp) - 2*KF(iza,jlp) - 2*KF(izb,jlp)
     &                 - 2*KF(izc,jlp) - 2*KF(izd,jlp) - 2*KF(ize,jlp)
     &                 - 2*KF(izf,jlp)
               sumlo = izexp
               qsumlo = qsumlo + sumlo*PILOG(jlp)*(.5E+00)
            ENDDO
            zuthre = vorz*EXP(qsumlo)
            wthrep = wthrep + zuthre
            vorz = -vorz
         ENDDO
         jvo = jj1 - jj2 - mm3
         IF ( 4*(jvo/4).NE.jvo ) wthrep = -wthrep
      ENDIF
      WTHREJ = wthrep
99999 END
