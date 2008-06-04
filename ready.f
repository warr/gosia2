 
C----------------------------------------------------------------------
C SUBROUTINE READY
C
C Called by: ADHOC
C Calls:     DECAY, SZEREG
C
C Purpose: To read experimental yields from a file.
C
C Uses global variables:
C      DYEX   - error on experimental yield
C      IY     - index of experimental yields
C      KSEQ   - index into ELM for pair of levels, and into EN or SPIN
C      LP6    - 32
C      NDST   - number of data sets
C      NEXPT  - number of experiments
C      NYLDE  - number of yields
C      YEXP   - experimental yield
C
C Formal parameters:
C      Idr    - number of decays
C      Ntap   - unit for yield file
C      Ipri   - printing flag (Ipri=1 gives additional output)
C
C NTAP is the unit number of the file from which we should read the
C experimental yields
 
      SUBROUTINE READY(Idr,Ntap,Ipri)
      IMPLICIT NONE
      REAL*8 ap , CORF , DYEX , EP , TAU , TLBDG , u , UPL , VINF , w , 
     &       waga , XA , XA1 , xep , YEXP , YNRM , zp
      INTEGER*4 idc , idc1 , idcx , Idr , IDRN , ii , ILE , Ipri , IY , 
     &          iytot , iytt , IZ , IZ1 , j , k , kk , kkl , KSEQ , 
     &          lbg , LP1
      INTEGER*4 LP10 , LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , 
     &          LP6 , LP7 , LP8 , LP9 , lxp , nanx , nde , nde1 , 
     &          ne , NEXPT , ns1
      INTEGER*4 ns2 , ns3 , ns4 , nsxh , nsyh , Ntap , nval , NYLDE
      DIMENSION iytot(32)
      COMMON /YEXPT / YEXP(32,1500) , IY(1500,32) , CORF(1500,32) , 
     &                DYEX(32,1500) , NYLDE(50,32) , UPL(32,50) , 
     &                YNRM(32,50) , IDRN , ILE(32)
      COMMON /CX    / NEXPT , IZ , XA , IZ1(50) , XA1(50) , EP(50) , 
     &                TLBDG(50) , VINF(50)
      COMMON /LEV   / TAU(75) , KSEQ(1500,4)
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      INCLUDE 'cccds.inc'
      
C     Rewind yield file
      REWIND Ntap

      DO k = 1 , LP6 ! LP6 = 32
         iytot(k) = 0
      ENDDO
      IF ( Ipri.EQ.1 ) WRITE (22,99001)
99001 FORMAT (5X/47X,'REPRINT OF EXPERIMENTAL DATA TO BE FITTED'//)
      DO lxp = 1 , NEXPT ! For each experiment
         DO kkl = 1 , LP6
            NYLDE(lxp,kkl) = 0
         ENDDO
         ii = NDST(lxp) ! Number of datasets
         DO kk = 1 , ii ! iexp, ng, zp, ag, ep, nd, wt
            READ (Ntap,*) ne , nanx , zp , ap , xep , nval , waga
            IF ( Ipri.EQ.1 ) WRITE (22,99002) ne , zp , ap , xep , 
     &                              NDST(ne) , waga
99002       FORMAT (1X,///5X,'EXPERIMENT',1X,1I2/2X,'PROJECTILE',1X,'(',
     &              1F4.0,',',1F4.0,')',1X,1F7.3,1X,'MEV',1X,'---',1I1,
     &              1X,'GE(LI) DETECTOR(S)',2X,'WEIGHT=',1E8.2/20X,
     &              '** EXPERIME','NTAL YIELDS **')
            IF ( Ipri.EQ.1 ) WRITE (22,99003)
99003       FORMAT (4X,'DECAY',1X,'IS',2X,'IF',1(9X,'YIELD+/-ERROR',9X)
     &              /)
            DO j = 1 , nval
               READ (Ntap,*) ns1 , ns2 , u , w ! ii, if, y, dy
               nsxh = ns1
               nsyh = ns2
               IF ( ns1.GE.100 ) THEN
                  ns1 = ns1/100
                  ns2 = ns2/100
               ENDIF
               DO nde = 1 , Idr
                  IF ( ns1.EQ.KSEQ(nde,3) .AND. ns2.EQ.KSEQ(nde,4) )
     &                 GOTO 10
               ENDDO
               IF ( Ipri.EQ.1 ) WRITE (22,99005) ns1 , ns2
               GOTO 40
 10            idc = nde
               iytot(kk) = iytot(kk) + 1
               idc1 = 0
               IF ( nsxh.GE.100 ) THEN
                  ns3 = nsxh - 100*ns1
                  ns4 = nsyh - 100*ns2
                  DO nde1 = 1 , Idr
                     IF ( ns3.EQ.KSEQ(nde1,3) .AND. ns4.EQ.KSEQ(nde1,4)
     &                    ) GOTO 20
                  ENDDO
                  IF ( Ipri.EQ.1 ) WRITE (22,99005) ns3 , ns4
               ENDIF
               GOTO 30
 20            idcx = idc*1000 + nde1
               IF ( idc.GT.nde1 ) idcx = nde1*1000 + idc
               idc = idcx
 30            idc1 = idc
               IF ( idc1.GT.1000 ) idc1 = idc/1000
               IF ( Ipri.EQ.1 ) WRITE (22,99004) idc , KSEQ(idc1,3) , 
     &                                 KSEQ(idc1,4) , u , w
99004          FORMAT (2X,1I6,2X,1I2,2X,1I2,1(1E14.6,3X,1E14.6))
               iytt = iytot(kk)
               YEXP(kk,iytt) = u
               DYEX(kk,iytt) = w/(SQRT(waga)+1.E-4)
               IY(iytt,kk) = idc
 40         ENDDO
            iytt = iytot(kk)
            lbg = iytt - nval + 1
            CALL SZEREG(lbg,iytt,kk)
            NYLDE(lxp,kk) = nval
         ENDDO ! For each dataset kk
      ENDDO ! Loop on experiments lxp
99005 FORMAT (1X///5X,'ERROR-NO MATRIX ELEMENT BETWEEN STATES',1X,1I2,
     &        ' AND ',1I2,/10X,'THIS TRANSITION IGNORED',//)
      END
