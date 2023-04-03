
C----------------------------------------------------------------------
C FUNCTION DJMM
C
C Called by: GOSIA, ROTATE, TENS
C
C Purpose: calculate the rotation functions D^k_{\xi \xi^\prime}
C
C Uses global variables:
C      B      - array of factorials
C      BEQ    - identifier for angle for rotations
C
C Formal parameters:
C      Beta   - v/c
C      K      - K
C      Kpp    - \xi
C      Kp     - \xi^\prime
C
C Return value:
C      Element of rotation matrix D
C
C Note that to be efficient, this function remembers values that it has
C previously calculated. For this to work, the variable djm has to be
C declared with a SAVE statement, otherwise the variable is an automatic one,
C which is created freshly each time the function is called.


      REAL*8 FUNCTION DJMM(Beta,K,Kpp,Kp)
      IMPLICIT NONE
      REAL*8 B , b1 , b2 , be , BEQ , Beta , cb , ctb , djm , f , g ,
     &       sb , sk , ul
      INTEGER*4 iczy , ifla , ifza , ill , j , ja , jb , jc , jd , K ,
     &          Kp , Kpp , lca , loc , mas , mis
      DIMENSION djm(525) , iczy(525)
      COMMON /IDENT / BEQ
      COMMON /CB    / B(20)
      SAVE djm ! Added N. Warr Jul2007

      ifza = 1
      IF ( Beta.LT.0. ) ifza = (-1)**(Kp+Kpp)
      sk = DBLE(K)
      ul = sk*((sk-1.)*(4.*sk+7)/6.+1.)
      lca = INT(ul+.1)

C     Calculate position in djm and iczy arrays
      loc = lca + (2*K+1)*Kp + Kpp + K + 1

      IF ( ABS(BEQ-ABS(Beta)).GT.1.E-6 ) THEN ! If beta doesn't match the identifier, initialise
         BEQ = ABS(Beta)
         DO ill = 1 , 525
            iczy(ill) = 0
         ENDDO
      ELSEIF ( iczy(loc).EQ.1 ) THEN ! We have already calculated it, so return
         DJMM = djm(loc)*ifza
         GOTO 99999
      ENDIF

C     We have to calculate it
      be = BEQ/2.
      cb = COS(be)
      sb = SIN(be)
      ifla = 0
      IF ( BEQ.GT..01 .AND. ABS(BEQ-6.2832).GT..01 ) ifla = 1
      IF ( ifla.NE.1 ) THEN
         IF ( Kp.EQ.Kpp ) THEN
            sb = 1.
         ELSE
            DJMM = 0.
            RETURN
         ENDIF
      ENDIF

      ctb = cb*cb/sb/sb
      ja = K + Kp + 1 ! K + \xi^\prime (+1 for indexing array)
      jb = K - Kp + 1 ! K - \xi^\prime (+1 for indexing array)
      jc = K + Kpp + 1 ! K + \xi (+1 for indexing array)
      jd = K - Kpp + 1 ! K - \xi (+1 for indexing array)
      b1 = B(ja)*B(jb)*B(jc)*B(jd) ! B array holds factorials

      ja = Kp + Kpp
      jb = 2*K - Kp - Kpp
      IF ( ABS(BEQ-3.141592654).LT..01 .AND. ja.LT.0 ) ifla = 3
      IF ( ifla.EQ.3 ) cb = 1.
      f = (-1)**(K-Kp)*(cb**ja)*(sb**jb)*SQRT(b1)
      mis = 0
      IF ( ja.LT.0 ) mis = -ja
      mas = K - Kpp
      IF ( Kpp.LT.Kp ) mas = K - Kp
      ja = Kp + Kpp + mis + 1
      jb = K - Kpp - mis + 1
      jc = K - Kp - mis + 1
      jd = mis + 1
      b2 = B(ja)*B(jb)*B(jc)*B(jd)
      IF ( ifla.NE.3 ) THEN
         g = (-ctb)**mis/b2
         DJMM = g
         ja = mis + 1
         IF ( mas.GE.ja ) THEN
            DO j = ja , mas
               g = -g*ctb*(K-Kpp-j+1)*(K-Kp-j+1)/(Kp+Kpp+j)/j
               DJMM = DJMM + g
            ENDDO
         ENDIF
         IF ( ifla.EQ.0 ) DJMM = g
         DJMM = DJMM*f*ifza
         djm(loc) = DJMM/ifza
         iczy(loc) = 1
         RETURN
      ENDIF
      DJMM = f*ifza/((-sb*sb)**mis)/b2
      djm(loc) = DJMM/ifza
      iczy(loc) = 1
99999 END
