      SUBROUTINE SEQ(Idr)
      IMPLICIT NONE
      REAL*8 ACCA , ACCUR , CONV , DELTA , DIPOL , ega , egs , emax , 
     &       EN , ENDEC , ENZ , F , FP , GF , GKP , SPIN , spinf , 
     &       spini , TAU , twoi
      REAL*8 ZPOL
      INTEGER*4 idecay , Idr , indx , inx , inx1 , ir , is , ISO , 
     &          istr1 , istr2 , ITMA , j , js , jsave , k , KLEC , kpa , 
     &          KSEQ , l , la
      INTEGER*4 la1 , LAMDA , LAMMAX , ld , LDNUM , LEAD , LEADF , LP1 , 
     &          LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , 
     &          LP6 , LP7 , LP8 , LP9
      INTEGER*4 m , m1 , m6 , MEM , mk , mule , mulm , MULTI , n , n1 , 
     &          NDIM , NMAX , NMAX1 , nob
      COMMON /COEX2 / NMAX , NDIM , NMAX1
      COMMON /TRA   / DELTA(500,3) , ENDEC(500) , ITMA(50,200) , 
     &                ENZ(200)
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      COMMON /COEX  / EN(75) , SPIN(75) , ACCUR , DIPOL , ZPOL , ACCA , 
     &                ISO
      COMMON /LEV   / TAU(75) , KSEQ(500,4)
      COMMON /CATLF / FP(4,500,3) , GKP(4,500,2) , KLEC(75)
      m6 = 0
      DO l = 1 , 6
         m6 = m6 + MULTI(l)
      ENDDO
      idecay = 0
      Idr = 0
      DO l = 1 , LP3
         KLEC(l) = 0
      ENDDO
      DO k = 1 , LP2
         DO j = 1 , 3
            DO l = 1 , 4
               FP(l,k,j) = 0.
               IF ( j.NE.3 ) GKP(l,k,j) = 0.
            ENDDO
            DELTA(k,j) = 0.
         ENDDO
      ENDDO
      DO n = 1 , NMAX
         TAU(n) = EN(n)
      ENDDO
      DO n = 1 , NMAX
         emax = 0.
         DO j = 1 , NMAX
            IF ( TAU(j).GE.emax ) THEN
               emax = TAU(j)
               jsave = j
            ENDIF
         ENDDO
         DO is = 1 , NMAX
            DO la = 1 , 8
               IF ( la.LE.3 .OR. la.EQ.7 .OR. la.EQ.8 ) THEN
                  ld = LDNUM(la,is)
                  IF ( ld.NE.0 ) THEN
                     DO ir = 1 , ld
                        m = LEADF(is,ir,la)
                        IF ( m.EQ.jsave .OR. is.EQ.jsave ) THEN
                           IF ( is.NE.jsave .OR. EN(m).LT.EN(is) ) THEN
                              IF ( m.NE.jsave .OR. EN(is).LT.EN(m) )
     &                             THEN
                                 indx = MEM(is,m,la)
                                 idecay = idecay + 1
                                 KSEQ(idecay,1) = m
                                 KSEQ(idecay,2) = is
                                 KSEQ(idecay,3) = indx
                                 KSEQ(idecay,4) = la + 10
                                 IF ( EN(m).LE.EN(is) ) THEN
                                    KSEQ(idecay,1) = is
                                    KSEQ(idecay,2) = m
                                 ENDIF
                              ENDIF
                           ENDIF
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
         TAU(jsave) = -1.
      ENDDO
      DO l = 1 , idecay
         istr1 = 0
         IF ( KSEQ(l,4).LT.10 ) GOTO 200
         istr2 = 0
         n = KSEQ(l,1)
         m = KSEQ(l,2)
         inx = KSEQ(l,3)
         la = KSEQ(l,4) - 10
         ega = EN(n) - EN(m)
         twoi = 1./SQRT(2.*SPIN(n)+1.)
         spini = SPIN(n) + .001
         spinf = SPIN(m) + .001
         egs = SQRT(ega)*twoi
         js = l + 1
         la1 = 0
         inx1 = 0
         DO j = js , idecay
            IF ( KSEQ(j,4).GE.10 ) THEN
               n1 = KSEQ(j,1)
               m1 = KSEQ(j,2)
               IF ( n1.EQ.n .AND. m1.EQ.m ) THEN
                  inx1 = KSEQ(j,3)
                  la1 = KSEQ(j,4) - 10
                  KSEQ(j,4) = KSEQ(j,4) - 10
               ENDIF
            ENDIF
         ENDDO
         KSEQ(l,4) = KSEQ(l,4) - 10
         Idr = Idr + 1
         mule = 0
         mulm = 0
         nob = 1
 50      IF ( la.LE.3 ) THEN
            IF ( la.EQ.1 ) THEN
               DELTA(Idr,1) = 399.05*ega*egs
               mule = 1
               istr1 = 1
            ELSEIF ( la.EQ.2 ) THEN
               DELTA(Idr,1) = 3.4928*egs*ega*ega
               mule = 2
               istr1 = 2
            ELSEIF ( la.EQ.3 ) THEN
               DELTA(Idr,1) = .02391*ega*ega*ega*egs
               mule = 3
               istr1 = 3
            ELSE
               GOTO 100
            ENDIF
            GOTO 150
         ENDIF
 100     la = la - 6
         IF ( la.EQ.2 ) THEN
            DELTA(Idr,2) = .0368*ega*ega*egs
            mulm = 2
            istr2 = 5
         ELSE
            DELTA(Idr,2) = 4.1952*ega*egs
            mulm = 1
            istr2 = 4
         ENDIF
 150     IF ( nob.NE.2 ) THEN
            IF ( mule.NE.1 ) THEN
               nob = nob + 1
               IF ( la.GT.3 ) inx1 = inx
               IF ( la1.NE.0 ) THEN
                  la = la1
                  GOTO 50
               ENDIF
            ENDIF
            inx1 = 0
         ENDIF
         DELTA(Idr,3) = DELTA(Idr,1)*DELTA(Idr,2)
         DELTA(Idr,1) = DELTA(Idr,1)*DELTA(Idr,1)
         DELTA(Idr,2) = DELTA(Idr,2)*DELTA(Idr,2)
         KSEQ(Idr,1) = inx
         KSEQ(Idr,2) = inx1
         KSEQ(Idr,3) = n
         KSEQ(Idr,4) = m
         IF ( inx.GT.m6 ) THEN
            KSEQ(Idr,2) = inx
            KSEQ(Idr,1) = 0
         ENDIF
         ENDEC(Idr) = EN(n) - EN(m)
         DO mk = 1 , 7 , 2
            kpa = mk/2 + 1
            k = mk - 1
            IF ( mule.GE.3 .OR. k.NE.6 ) THEN
               GKP(kpa,Idr,1) = GF(k,spini,spinf,mule)*DELTA(Idr,1)
     &                          *(1.+CONV(ega,istr1))
               GKP(kpa,Idr,2) = GF(k,spini,spinf,mulm)*DELTA(Idr,2)
     &                          *(1.+CONV(ega,istr2))
               FP(kpa,Idr,1) = F(k,spini,spinf,mule,mule)*DELTA(Idr,1)
               FP(kpa,Idr,3) = F(k,spini,spinf,mulm,mule)*DELTA(Idr,3)
               FP(kpa,Idr,2) = F(k,spini,spinf,mulm,mulm)*DELTA(Idr,2)
            ENDIF
         ENDDO
         DELTA(Idr,1) = DELTA(Idr,1)*(1.+CONV(ega,istr1))
         DELTA(Idr,2) = DELTA(Idr,2)*(CONV(ega,istr2)+1.)
         KLEC(n) = KLEC(n) + 1
 200  ENDDO
      NMAX1 = 0
      DO n = 1 , NMAX
         IF ( KLEC(n).NE.0 ) NMAX1 = NMAX1 + 1
      ENDDO
      END
