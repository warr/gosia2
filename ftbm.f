      SUBROUTINE FTBM(Icll,Chisq,Idr,Ncall,Chilo,Bten)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , AGELI , aval , Bten , CAT , CC , Chilo , 
     &       chis1 , CHIS11 , chish , Chisq , chisx , chx , CORF , 
     &       DIPOL , DYEX , EG , ELM , ELML
      REAL*8 ELMU , EMH , EN , EP , EPS , EROOT , fc , FIEX , fx , 
     &       polm , pr , prop , Q , SA , SPIN , TAU , TLBDG , UPL , 
     &       val , VINF
      REAL*8 wz , XA , XA1 , YEXP , YNRM , ZETA , ZPOL
      INTEGER*4 i1 , i11 , iapx , IAXS , Icll , idec , Idr , IDRN , 
     &          IEXP , iflg , IGRD , ii , ILE , ile1 , ile2 , ile3 , 
     &          ilin , indx , inko , INM
      INTEGER*4 inp , inpo , inpx , INTR , inzz , inzzz , IPATH , IPRM , 
     &          IPS1 , ISMAX , ISO , issp , ITAK2 , itemp , IVAR , ixx , 
     &          IY , IZ , IZ1 , izzz
      INTEGER*4 j , jj , jjgg , jjj , jk , jkl , jm , jmf , jmt , jmte , 
     &          jpp , jpz , JSKIP , jy , k , karm , kk , kk6 , kkx , kmt
      INTEGER*4 knm , KSEQ , kx , larm , lcc , lcou , LFL , LFL1 , 
     &          LFL2 , licz , lix , llx , lm , LMAX , LMAXE , lmh , 
     &          LNY , loc , loch , loct
      INTEGER*4 lp , LP1 , LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , 
     &          LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , lpit , lput , lpx , 
     &          lpxd , ls , lst
      INTEGER*4 luu , lx , LZETA , MAGA , MAGEXC , MEMAX , MEMX6 , 
     &          NANG , Ncall , NDIM , NEXPT , NICC , NLIFT , nlin , 
     &          NMAX , NMAX1 , nowr , npoz , nrest , NSTART
      INTEGER*4 NSTOP , NWR , nwyr , NYLDE
      COMPLEX*16 ARM
      DIMENSION jmte(6) , prop(6) , Bten(1200)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      COMMON /CCC   / EG(50) , CC(50,5) , AGELI(50,200,2) , Q(3,200,8) , 
     &                NICC , NANG(200)
      COMMON /ILEWY / NWR
      COMMON /CH1T  / CHIS11
      COMMON /IGRAD / IGRD
      COMMON /LCZP  / EMH , INM , LFL1 , LFL2 , LFL
      COMMON /UWAGA / ITAK2
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /KIN   / EPS(50) , EROOT(50) , FIEX(50,2) , IEXP , IAXS(50)
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /CLM   / LMAX
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /AZ    / ARM(600,7)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /PTH   / IPATH(75) , MAGA(75)
      COMMON /PRT   / IPRM(20)
      COMMON /CEXC  / MAGEXC , MEMAX , LMAXE , MEMX6 , IVAR(500)
      COMMON /SKP   / JSKIP(50)
      COMMON /LIFE  / NLIFT
      COMMON /LOGY  / LNY , INTR , IPS1
      issp = 0
      Chilo = 0.
      fx = 2.*SPIN(1) + 1.
      Chisq = 0.
      LFL = 0
      chis1 = 0.
      ixx = NDIM*MEMAX + LP11
      DO i1 = 1 , ixx
         ZETA(i1) = 0.
      ENDDO
      DO ii = 1 , LP6
         ILE(ii) = 1
      ENDDO
      itemp = 0
      NWR = 0
      iapx = 1
      DO jkl = 1 , NEXPT
         IEXP = jkl
         IGRD = 0
         LFL2 = 1
         IF ( ITAK2.EQ.-1 ) THEN
            DO larm = 1 , 4
               DO karm = 1 , LP10
                  ARM(karm,larm) = (0.,0.)
               ENDDO
            ENDDO
         ENDIF
         iflg = 0
         IF ( IEXP.NE.1 ) THEN
            kk = NANG(IEXP)
            DO jjj = 1 , LP6
               ILE(jjj) = ILE(jjj) + NYLDE(IEXP-1,jjj)
            ENDDO
         ENDIF
         lp = 3
         IF ( JSKIP(jkl).EQ.0 ) GOTO 200
         IF ( MAGA(IEXP).EQ.0 ) lp = 1
         IF ( Ncall.EQ.0 ) GOTO 150
         IF ( Icll.EQ.4 ) GOTO 100
 50      loch = LP3*(MEMAX-1) + NMAX + LP11
         DO k = 1 , loch
            ZETA(k) = 0.
         ENDDO
         CALL LOAD(IEXP,1,2,0.D0,jj)
         DO k = 1 , LMAX
            fc = 2.
            IF ( k.EQ.LMAX ) fc = 1.
            IF ( DBLE(INT(SPIN(1))).LT.SPIN(1) ) fc = 2.
            loc = 0
            polm = DBLE(k-1) - SPIN(1)
            CALL LOAD(IEXP,3,2,polm,jj)
            CALL PATH(jj)
            CALL LOAD(IEXP,2,2,polm,jj)
            CALL APRAM(IEXP,1,1,jj,ACCA)
            IF ( Ncall.NE.0 ) THEN
               IF ( Icll.NE.3 ) THEN
                  DO indx = 1 , MEMX6
                     CALL APRAM(IEXP,0,indx,jj,ACCA)
                     kx = 0
                     DO i11 = 1 , NMAX
                        IF ( NSTART(i11).NE.0 ) THEN
                           loc = LP3*(indx-1) + i11 + LP11
                           jpp = INT(2.*SPIN(i11)+1.)
                           lpx = MIN(lp,jpp)
                           IF ( ISO.NE.0 ) lpx = NSTOP(i11)
     &                          - NSTART(i11) + 1
                           DO lpxd = 1 , lpx
                              kx = kx + 1
                              ZETA(loc) = ZETA(loc) + fc*DBLE(ARM(kx,5))
     &                           *DBLE(ARM(kx,6))
     &                           /fx + fc*IMAG(ARM(kx,5))
     &                           *IMAG(ARM(kx,6))/fx
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
            CALL TENB(k,Bten,LMAX)
         ENDDO
         IF ( loc.NE.0 ) THEN
            REWIND 14
            WRITE (14,*) (ZETA(i11),i11=LP8,loch)
         ENDIF
         CALL TENS(Bten)
         IF ( Ncall.EQ.0 ) GOTO 200
         IF ( Icll.GE.2 ) GOTO 200
         llx = 28*NMAX
         DO lx = 1 , llx
            ZETA(LP9+lx) = ZETA(lx)
         ENDDO
         IF ( Icll.NE.1 ) GOTO 200
 100     iapx = 0
         issp = 1
         CALL LOAD(IEXP,1,1,0.D0,jj)
         CALL ALLOC(ACCUR)
         CALL SNAKE(IEXP,ZPOL)
         CALL SETIN
         DO k = 1 , LMAX
            polm = DBLE(k-1) - SPIN(1)
            CALL LOAD(IEXP,2,1,polm,kk)
            IF ( IPRM(7).EQ.-1 ) WRITE (22,99001) polm , IEXP
99001       FORMAT (1X//40X,'EXCITATION AMPLITUDES'//10X,'M=',1F4.1,5X,
     &              'EXPERIMENT',1X,1I2//5X,'LEVEL',2X,'SPIN',2X,'M',5X,
     &              'REAL AMPLITUDE',2X,'IMAGINARY AMPLITUDE'//)
            CALL STING(kk)
            CALL PATH(kk)
            CALL INTG(IEXP)
            CALL TENB(k,Bten,LMAX)
            IF ( IPRM(7).EQ.-1 ) THEN
               DO j = 1 , ISMAX
                  WRITE (22,99002) INT(CAT(j,1)) , CAT(j,2) , CAT(j,3) , 
     &                             DBLE(ARM(j,5)) , IMAG(ARM(j,5))
99002             FORMAT (7X,1I2,3X,1F4.1,2X,1F4.1,2X,1E14.6,2X,1E14.6)
               ENDDO
            ENDIF
         ENDDO
         CALL TENS(Bten)
         IF ( IPRM(7).EQ.-1 ) THEN
            DO jjgg = 2 , NMAX
               loct = (jjgg-1)*28 + 1
               WRITE (22,99003) jjgg , ZETA(loct)
99003          FORMAT (2X,'LEVEL',1X,1I2,10X,'POPULATION',1X,1E14.6)
            ENDDO
         ENDIF
         GOTO 200
 150     IF ( iflg.EQ.1 ) THEN
            itemp = 1
            iflg = 2
            GOTO 50
         ELSE
            IF ( iflg.EQ.2 ) GOTO 300
            itemp = 2
            iflg = 1
            GOTO 100
         ENDIF
 200     CALL CEGRY(Chisq,itemp,Chilo,Idr,nwyr,Icll,issp,0)
         issp = 0
         IF ( Ncall.EQ.0 .AND. JSKIP(jkl).NE.0 ) THEN
            IF ( Ncall.NE.0 ) GOTO 200
            GOTO 150
         ELSE
            NWR = NWR + nwyr
            IF ( Icll.LE.2 .AND. JSKIP(jkl).NE.0 ) THEN
               IF ( IEXP.EQ.1 ) chish = CHIS11
               IF ( Icll.EQ.1 ) chis1 = CHIS11
               IF ( Icll.EQ.0 ) chis1 = Chisq
               LFL2 = 0
               IGRD = 1
               IF ( ITAK2.EQ.-1 ) LFL = 1
               REWIND 14
               READ (14,*) (ZETA(i11),i11=LP8,loch)
               DO larm = 1 , 4
                  DO karm = 1 , LP10
                     ARM(karm,larm) = (0.,0.)
                  ENDDO
               ENDDO
               chisx = 0.
               llx = 28*NMAX
               DO lix = 1 , llx
                  ZETA(LP9+lix) = ZETA(lix)
               ENDDO
               CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,1)
               DO knm = 1 , MEMAX
                  INM = knm
                  chisx = 0.
                  EMH = ELM(INM)
                  ELM(INM) = 1.05*EMH
                  lcc = LP3*(INM-1) + LP11
                  DO lst = 2 , NMAX
                     wz = ZETA(lst+lcc)
                     inpx = (lst-1)*28
                     DO jy = 1 , 4
                        inp = inpx + (jy-1)*7
                        IF ( jy.EQ.1 ) pr = ZETA(LP13+inp) + 1.E-12
                        jmf = 2*jy - 1
                        IF ( IAXS(IEXP).EQ.0 ) jmf = 1
                        DO jm = 1 , jmf
                           inp = inp + 1
                           ZETA(inp) = ZETA(inp+LP9)*(1.+.1*EMH*wz/pr)
                        ENDDO
                     ENDDO
                  ENDDO
                  CALL CEGRY(chisx,itemp,Chilo,Idr,nwyr,0,0,0)
                  ELM(INM) = EMH
               ENDDO
               IF ( ITAK2.EQ.-1 .AND. LFL1.NE.0 ) THEN
                  IF ( IPRM(17).NE.0 ) THEN
                     kmt = ABS(IPRM(17))
                     WRITE (22,99004) IEXP
99004                FORMAT (1X///20X,'EXPERIMENT',11X,1I2,5X,
     &                       'D(LOG(P))/D(LOG(ME)) MAP'/20X,52('-')///)
                     nlin = (NMAX-2)/6 + 1
                     nrest = NMAX - 1 - 6*(nlin-1)
                     DO ilin = 1 , nlin
                        npoz = 6
                        IF ( ilin.EQ.nlin ) npoz = nrest
                        inpo = (ilin-1)*6 + 2
                        inko = inpo + npoz - 1
                        lpit = 0
                        DO lm = inpo , inko
                           lpit = lpit + 1
                           jmte(lpit) = lm
                        ENDDO
                        WRITE (22,99005) (jmte(lm),lm=1,lpit)
99005                   FORMAT (5X,'LEVEL',6(8X,1I2,9X))
                        WRITE (22,99006)
     &                         (ZETA(LP13+(jpz-1)*28),jpz=inpo,inko)
99006                   FORMAT (1X,'EXC.PROB.',6(5X,1E10.4,4X))
                        DO jmt = 1 , kmt
                           lput = 0
                           DO ls = inpo , inko
                              lput = lput + 1
                              prop(lput) = 0.
                              DO lm = 1 , MEMX6
                                 inzz = ls + LP3*(lm-1) + LP11
                                 inzzz = LP13 + (ls-1)*28
                                 IF ( ABS(ZETA(inzzz)).LT.1.E-20 )
     &                                ZETA(inzzz) = 1.E-20
                                 val = 2.*ELM(lm)*ZETA(inzz)/ZETA(inzzz)
                                 aval = ABS(val)
                                 IF ( aval.GT.ABS(prop(lput)) ) THEN
                                    prop(lput) = val
                                    lmh = lm
                                    jmte(lput) = lm
                                 ENDIF
                              ENDDO
                              izzz = (lmh-1)*LP3 + LP11 + ls
                              ZETA(izzz) = 0.
                           ENDDO
                           WRITE (22,99007)
     &                            (jmte(lcou),prop(lcou),lcou=1,npoz)
99007                      FORMAT (10X,6(2X,'(',1X,1I3,1X,1E8.2,')',2X))
                        ENDDO
                     ENDDO
                     REWIND 14
                     READ (14,*) (ZETA(i11),i11=LP8,loch)
                     IF ( IPRM(17).LT.0 ) GOTO 300
                  ENDIF
                  LFL = 0
                  WRITE (22,99008) IEXP
99008             FORMAT (10X,'EXPERIMENT',1X,1I2/10X,
     &                    'D(LOG(Y)/D(LOG(ME))',//)
                  ile1 = ILE(1) + NYLDE(IEXP,1) - 1
                  ile3 = ILE(1)
                  licz = 0
                  DO ile2 = ile3 , ile1
                     licz = licz + 1
                     idec = IY(ile2,1)
                     IF ( idec.GT.1000 ) idec = idec/1000
                     luu = 6*licz - 5
                     jk = (luu-1)/LP10 + 1
                     kk = luu - LP10*(jk-1)
                     kk6 = kk + 5
                     WRITE (22,99009) KSEQ(idec,3) , KSEQ(idec,4) , 
     &                                (INT(DBLE(ARM(kkx,jk))),
     &                                IMAG(ARM(kkx,jk)),kkx=kk,kk6)
99009                FORMAT (2X,1I2,'--',1I2,5X,
     &                       6('(',1I3,2X,1E8.2,')',3X))
                  ENDDO
               ENDIF
            ENDIF
         ENDIF
 300  ENDDO
      IF ( ITAK2.EQ.-1 .AND. Icll.LT.2 ) ITAK2 = 0
      IF ( Ncall.NE.0 ) THEN
         IF ( Icll.LE.2 ) THEN
            IF ( Icll.EQ.1 ) CALL CEGRY(Chisq,itemp,Chilo,Idr,nowr,7,
     &                                  issp,0)
         ENDIF
         CALL BRANR(Chisq,NWR,Chilo)
         CALL MIXR(NWR,0,Chisq,Chilo)
         CALL CHMEM(NWR,Chisq,Chilo)
         NWR = NWR + NLIFT
         Chisq = Chisq/NWR
         IF ( INTR.NE.0 ) THEN
            chx = Chisq
            Chisq = Chilo
            Chilo = chx
         ENDIF
      ENDIF
      END