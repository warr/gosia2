 
C----------------------------------------------------------------------
C FUNCTION WSIXJ
C
C Called by: F, GF, GKK, GOSIA
C
C Purpose: evaluates a Wigner 6-j symbol.
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
C      L1     - twice the value of L1
C      L2     - twice the value of L2
C      L3     - twice the value of L3
C
C Return value:
C      The value of the 6-j symbol
C
C Note that the values of the parameters are doubled, so that this function
C can handle half-integers. In other words if you want to evaluate
C \sixj(J1 J2 J3 L1 L2 L3) you need to use call the function as:
C WSIXJ(2 * J1, 2 * J2, 2 * J3, 2 * L1, 2 * L2, 2 * L3).
 
      REAL*8 FUNCTION WSIXJ(J1,J2,J3,L1,L2,L3)
      IMPLICIT NONE
      INTEGER*4 irj , irl , isa , isb , isc , isumfa , iva , ivb , 
     &          ivc , ivd , ivorfa , iz , iza , izb , izc , izd , 
     &          ize , izf
      INTEGER*4 izg , izh , izmax , izmin , J1 , J2 , J3 , kqa , 
     &          kqb , kqc , kqd , kra , krb , krc , krd , ksa , ksb , 
     &          ksc , ksd
      INTEGER*4 kta , ktb , ktc , ktd , kua , kub , kuc , L1 , L2 , L3 , 
     &          n , nmax
      REAL*8 qsumfa , qsumlo , sumlo , vorz , wsixp , zusix
      DIMENSION isumfa(26) , ivorfa(26)
      INCLUDE 'fakul.inc'
      
      wsixp = 0.E+00
      IF ( ((J1+J2-J3).GE.0) .AND. ((J1-J2+J3).GE.0) .AND. 
     &     ((-J1+J2+J3).GE.0) ) THEN
         IF ( ((J1+L2-L3).GE.0) .AND. ((J1-L2+L3).GE.0) .AND. 
     &        ((-J1+L2+L3).GE.0) ) THEN
            IF ( ((L1+J2-L3).GE.0) .AND. ((L1-J2+L3).GE.0) .AND. 
     &           ((-L1+J2+L3).GE.0) ) THEN
               IF ( ((L1+L2-J3).GE.0) .AND. ((L1-L2+J3).GE.0) .AND. 
     &              ((-L1+L2+J3).GE.0) ) THEN
                  kqa = (J1+J2-J3)/2
                  kqb = (J1-J2+J3)/2
                  kqc = (J2+J3-J1)/2
                  kqd = (J1+J2+J3)/2
                  kra = (J1+L2-L3)/2
                  krb = (J1-L2+L3)/2
                  krc = (L2+L3-J1)/2
                  krd = (J1+L2+L3)/2
                  ksa = (L1+J2-L3)/2
                  ksb = (L1-J2+L3)/2
                  ksc = (J2+L3-L1)/2
                  ksd = (L1+J2+L3)/2
                  kta = (L1+L2-J3)/2
                  ktb = (L1-L2+J3)/2
                  ktc = (L2+J3-L1)/2
                  ktd = (L1+L2+J3)/2
                  izmin = MAX(kqd,krd,ksd,ktd)
                  kua = kqa + kta + J3
                  kub = ksc + ktc + L1
                  kuc = krb + ktb + L2
                  izmax = MIN(kua,kub,kuc)
                  IF ( izmin.LE.izmax ) THEN
                     nmax = MAX(izmax+2,kqd+2,krd+2,ksd+2,ktd+2)
                     DO n = 1 , 26
                        IF ( IP(n).GE.nmax ) GOTO 5
                     ENDDO
                  ENDIF
                  GOTO 100
 5                vorz = -1.E+00
                  IF ( 2*(izmin/2).EQ.izmin ) vorz = +1.E+00
                  DO irl = 1 , n
                     iva = KF(kqa+1,irl) + KF(kqb+1,irl) + KF(kqc+1,irl)
     &                     - KF(kqd+2,irl)
                     ivb = KF(kra+1,irl) + KF(krb+1,irl) + KF(krc+1,irl)
     &                     - KF(krd+2,irl)
                     ivc = KF(ksa+1,irl) + KF(ksb+1,irl) + KF(ksc+1,irl)
     &                     - KF(ksd+2,irl)
                     ivd = KF(kta+1,irl) + KF(ktb+1,irl) + KF(ktc+1,irl)
     &                     - KF(ktd+2,irl)
                     ivorfa(irl) = iva + ivb + ivc + ivd
                  ENDDO
                  DO iz = izmin , izmax
                     sumlo = 0.E+00
                     iza = iz + 2
                     izb = iz - kqd + 1
                     izc = iz - krd + 1
                     izd = iz - ksd + 1
                     ize = iz - ktd + 1
                     izf = kua - iz + 1
                     izg = kub - iz + 1
                     izh = kuc - iz + 1
                     DO irj = 1 , n
                        isa = 2*KF(iza,irj) - 2*KF(izb,irj)
     &                        - 2*KF(izc,irj)
                        isb = -2*KF(izd,irj) - 2*KF(ize,irj)
     &                        - 2*KF(izf,irj)
                        isc = ivorfa(irj) - 2*KF(izg,irj)
     &                        - 2*KF(izh,irj)
                        isumfa(irj) = isa + isb + isc
                        qsumfa = isumfa(irj)
                        sumlo = sumlo + qsumfa*PILOG(irj)
                     ENDDO
                     qsumlo = (.5E+00)*sumlo
                     zusix = EXP(qsumlo)*vorz
                     wsixp = wsixp + zusix
                     vorz = -vorz
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
      ENDIF
 100  WSIXJ = wsixp
      END
