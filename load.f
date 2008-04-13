      SUBROUTINE LOAD(Iexp,Ient,Icg,Polm,Joj)
      IMPLICIT NONE
      REAL*8 a1 , a2 , aaz2 , aaz3 , aazz , ACCA , ACCUR , ah , CAT , 
     &       cpsi , dep , DIPOL , EMMA , EN , EP , eta , etan , Polm , 
     &       pp1 , pp2
      REAL*8 ppp , PSI , QAPR , rlam , SPIN , ssqrt , szet , TLBDG , 
     &       VINF , wrt , wrtm , XA , XA1 , XI , z1 , z2 , zet , ZETA , 
     &       ZPOL , zsqa
      INTEGER*4 i , i1 , i2 , i3 , IAPR , Icg , Ient , Iexp , IPATH , 
     &          ir , is , ISEX , ISHA , ISMAX , ISO , ispi , ispo , 
     &          IVAR , IZ , IZ1
      INTEGER*4 jj , jjj , Joj , la , lam , lam1 , LAMDA , LAMMAX , ld , 
     &          LDNUM , LEAD , LMAX , LMAXE , LP1 , LP10 , LP11 , LP12 , 
     &          LP13 , LP14 , LP2
      INTEGER*4 LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , LZETA , m , m1 , 
     &          m2 , MAGA , MAGEXC , MEMAX , MEMX6 , mstop , MULTI , n , 
     &          n2 , n3 , NCM
      INTEGER*4 NDIM , NEXPT , NMAX , NMAX1 , nn , NSTART , NSTOP , nz
      LOGICAL ERR
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /PSPIN / ISHA
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /PCOM  / PSI(500)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CLM   / LMAX
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /CLCOM9/ ERR
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /CXI   / XI(500)
      COMMON /CAUX0 / EMMA(75) , NCM
      COMMON /APRCAT/ QAPR(500,2,7) , IAPR(500,2) , ISEX(75)
      COMMON /PTH   / IPATH(75) , MAGA(75)
      DIMENSION etan(75) , cpsi(8)
      LMAX = INT(SPIN(1)+1.1)
      IF ( Ient.EQ.1 ) THEN
         ISHA = 0
         ispi = INT(SPIN(1)+.51)
         ispo = INT(SPIN(1)+.49)
         IF ( ispi.NE.ispo ) ISHA = 1
         z1 = DBLE(ABS(IZ1(Iexp)))
         z2 = DBLE(IZ)
         a1 = XA1(Iexp)
         a2 = XA
         ZPOL = DIPOL*EP(Iexp)*a2/(z2*z2*(1.+a1/a2))
         IF ( IZ1(Iexp).LT.0 ) ZPOL = DIPOL*EP(Iexp)
     &                                *a1/(z1*z1*(1.+a2/a1))
         IF ( IZ1(Iexp).LE.0 ) THEN
            ah = a1
            a1 = a2
            a2 = ah
         ENDIF
         eta = z1*z2*SQRT(a1/EP(Iexp))/6.349770
         DO m = 1 , NMAX
            dep = (1.0+a1/a2)*EN(m)
            zet = dep/EP(Iexp)
            szet = SQRT(1.0-zet)
            etan(m) = eta/szet
         ENDDO
         DO n = 1 , MEMAX
            i1 = LEAD(1,n)
            i2 = LEAD(2,n)
            XI(n) = etan(i1) - etan(i2)
         ENDDO
         aazz = 1./(1.+a1/a2)/z1/z2
         cpsi(1) = 5.169286*aazz
         IF ( LMAXE.NE.1 ) THEN
            aaz2 = aazz*aazz
            cpsi(2) = 14.359366*aaz2
            IF ( LMAXE.NE.2 ) THEN
               aaz3 = aazz*aaz2
               cpsi(3) = 56.982577*aaz3
               IF ( LMAXE.NE.3 ) THEN
                  aazz = aaz2*aaz2
                  cpsi(4) = 263.812653*aazz
                  IF ( LMAXE.NE.4 ) THEN
                     aaz2 = aaz3*aaz2
                     cpsi(5) = 1332.409500*aaz2
                     IF ( LMAXE.NE.5 ) THEN
                        aazz = aaz3*aaz3
                        cpsi(6) = 7117.691577*aazz
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF ( MAGEXC.NE.0 ) THEN
            aazz = VINF(Iexp)/95.0981942
            cpsi(7) = aazz*cpsi(1)
            IF ( LAMMAX.NE.8 ) cpsi(8) = aazz*cpsi(2)
         ENDIF
         zsqa = z1*SQRT(a1)
         i3 = 1
         ppp = 1. + a1/a2
         DO i1 = 1 , LAMMAX
            lam = LAMDA(i1)
            lam1 = lam
            IF ( lam.GT.6 ) lam1 = lam - 6
            DO n2 = 1 , NMAX
               nn = LDNUM(lam,n2)
               IF ( nn.NE.0 ) THEN
                  n3 = LEAD(1,i3)
                  pp1 = EP(Iexp) - ppp*EN(n3)
                  DO m1 = 1 , nn
                     m2 = LEAD(2,i3)
                     i2 = i3
                     i3 = i3 + 1
                     pp2 = EP(Iexp) - ppp*EN(m2)
                     PSI(i2) = cpsi(lam)*zsqa*(pp1*pp2)
     &                         **((2.*DBLE(lam1)-1.)/4.)
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
         IF ( Ient.EQ.1 ) RETURN
      ENDIF
      DO n = 1 , NMAX
         NSTART(n) = 0
         NSTOP(n) = 0
      ENDDO
      is = 1
      NSTART(1) = 1
      DO n = 1 , NMAX
         wrt = Polm - EMMA(Iexp)
         wrtm = Polm + EMMA(Iexp)
         IF ( Icg.EQ.2 ) wrt = Polm - DBLE(MAGA(Iexp))
         IF ( Icg.EQ.2 ) wrtm = Polm + DBLE(MAGA(Iexp))
         IF ( wrtm.LT.-SPIN(n) ) THEN
            NSTART(n) = 0
         ELSE
            IF ( ABS(wrt).GT.SPIN(n) ) wrt = -SPIN(n)
            IF ( wrtm.GT.SPIN(n) ) wrtm = SPIN(n)
            mstop = INT(wrtm-wrt+1.01)
            DO i = 1 , mstop
               CAT(is,1) = n
               CAT(is,2) = SPIN(n)
               CAT(is,3) = wrt + DBLE(i-1)
               IF ( n.EQ.1 .AND. ABS(CAT(is,3)-Polm).LT.1.E-6 ) Joj = is
               is = is + 1
            ENDDO
         ENDIF
         NSTART(n+1) = is
         NSTOP(n) = is - 1
      ENDDO
      ISMAX = is - 1
      IF ( ISMAX.LE.LP10 ) THEN
         IF ( Ient.EQ.3 ) RETURN
         nz = 0
         DO jj = 1 , 7
            DO jjj = 1 , MEMAX
               QAPR(jjj,1,jj) = 0.
               QAPR(jjj,2,jj) = 0.
            ENDDO
         ENDDO
         DO i = 1 , 8
            LZETA(i) = 0
         ENDDO
         DO i1 = 1 , LAMMAX
            lam = LAMDA(i1)
            IF ( Icg.NE.2 .OR. lam.LE.6 ) THEN
               la = lam
               IF ( lam.GT.6 ) lam = lam - 6
               rlam = DBLE(lam)
               ssqrt = SQRT(2.*rlam+1.)
               LZETA(la) = nz
               ir = 0
 10            ir = ir + 1
               IF ( ir.LE.ISMAX ) THEN
                  n = CAT(ir,1)
                  IF ( Icg.NE.1 ) THEN
                     IF ( MAGA(Iexp).EQ.0 .AND. ir.NE.IPATH(n) ) GOTO 10
                     IF ( ABS(ir-IPATH(n)).GT.1 ) GOTO 10
                  ENDIF
                  ld = LDNUM(la,n)
                  IF ( ld.EQ.0 ) THEN
                     ir = ir + NSTOP(n) - NSTART(n)
                  ELSE
                     CALL LSLOOP(ir,n,nz,ld,lam,la,ssqrt,Icg,Iexp)
                  ENDIF
                  GOTO 10
               ENDIF
            ENDIF
         ENDDO
         IF ( nz.GT.LP7 ) THEN
            WRITE (22,99001) LP7
99001       FORMAT (1x,
     &              'ERROR - NUMBER OF ELEMENTS IN ZETA ARRAY EXCEEDS',
     &              'ZEMAX',5X,'(ZEMAX =',I6,')')
         ELSE
            RETURN
         ENDIF
      ELSE
         WRITE (22,99002) LP10
99002    FORMAT (' ERROR-ISMAX EXCEEDS MAGMAX',5X,'(MAGMAX =',I4,')')
      ENDIF
      ERR = .TRUE.
      END