      SUBROUTINE XSTATIC(Iz,Ido,Iup,Beta)
      IMPLICIT NONE
      REAL*8 AKS , Beta , DQ , h , QCEN , VACDP , XNOR
      INTEGER*4 IBYP , Ido , Iup , Iz , lq
      COMMON /VAC   / VACDP(3,75) , QCEN , DQ , XNOR , AKS(6,75) , IBYP
      h = 1./(1.+(Iz**.45*.012008/Beta)**1.666667)
      QCEN = Iz*h**.6
      DQ = SQRT(QCEN*(1.-h))/2.
      Iup = INT(QCEN+3.*DQ+.5)
      Ido = INT(QCEN-3.*DQ-.5)
      IF ( Iup.GT.Iz ) Iup = Iz
      IF ( Ido.LT.1 ) Ido = 1
      XNOR = 0.
      DO lq = Ido , Iup
         XNOR = XNOR + EXP(-((QCEN-DBLE(lq))/DQ)**2/2.)
      ENDDO
      END
