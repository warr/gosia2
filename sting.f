      SUBROUTINE STING(Irld)
      IMPLICIT NONE
      REAL*8 CAT , D2W , ELM , ELML , ELMU , rsg , SA , w0 , ZETA
      INTEGER*4 i , i57 , ibg , iend , IFLG , indx , IRA , Irld , is2 , 
     &          ISG , ISG1 , ISMAX , ISSTAR , ISSTO , j , j1 , jj , 
     &          KDIV , lam , LAMDA
      INTEGER*4 LAMMAX , LAMR , ld , LDNUM , LEAD , LZETA , maxh , 
     &          MAXLA , mm , MSTORE , MULTI , n , NDIV , NPT , NSW , nz
      COMPLEX*16 ARM , EXPO
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /AZ    / ARM(600,7)
      COMMON /COMME / ELM(500) , ELMU(500) , ELML(500) , SA(500)
      COMMON /ADBXI / EXPO(500)
      COMMON /FLA   / IFLG
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      COMMON /CCOUP / ZETA(50000) , LZETA(8)
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /CLCOM8/ CAT(600,3) , ISMAX
      COMMON /RNG   / IRA(8) , MAXLA
      maxh = MAXLA
 100  ISG = -1
      n = 1
      rsg = -1.
      IFLG = 1
      w0 = IRA(MAXLA)*.03 + .03
      DO j = 1 , ISMAX
         DO jj = 1 , 6
            ARM(j,jj) = (0.,0.)
         ENDDO
      ENDDO
      ARM(Irld,5) = (1.,0.)
      DO j = 1 , 8
         LAMR(j) = 0
      ENDDO
      LAMR(MAXLA) = 1
      NPT = IRA(MAXLA) + 1
      IF ( MAXLA.EQ.7 .AND. IRA(2).NE.0 ) THEN
         LAMR(2) = 1
         NPT = NPT - 1
         w0 = w0 - .03
      ENDIF
      NDIV = 0
      KDIV = 0
      DO j = 1 , 4
         NPT = NPT - 1
         DO j1 = 1 , LAMMAX
            lam = LAMDA(j1)
            IF ( LAMR(lam).NE.0 ) THEN
               CALL NEWLV(n,ld,lam)
               IF ( ld.NE.0 ) THEN
                  nz = LZETA(lam)
                  ld = LDNUM(lam,1)
                  i57 = 5
                  CALL LAISUM(Irld,n,rsg,lam,ld,nz,i57)
                  DO mm = 1 , ld
                     indx = MSTORE(2,mm)
                     ibg = ISSTAR(mm)
                     iend = ISSTO(mm)
                     DO is2 = ibg , iend
                        ARM(is2,4) = ARM(is2,4) + ARM(is2,6)*ELM(indx)
     &                               /EXPO(indx)
                        ARM(is2,6) = (0.,0.)
                     ENDDO
                  ENDDO
               ELSEIF ( j1.EQ.MAXLA ) THEN
                  IRA(MAXLA) = -IRA(MAXLA)
                  DO jj = 1 , LAMMAX
                     lam = LAMDA(jj)
                     IF ( IRA(lam).GT.0 ) GOTO 105
                  ENDDO
 105              MAXLA = LAMDA(jj)
                  GOTO 100
               ENDIF
            ENDIF
         ENDDO
         IF ( j.EQ.4 ) GOTO 200
         DO i = 1 , ISMAX
            ARM(i,j) = ARM(i,4)
            ARM(i,4) = (0.,0.)
         ENDDO
      ENDDO
 200  CALL LAIAMP(Irld,w0)
      MAXLA = maxh
      DO jj = 1 , 8
         IRA(jj) = ABS(IRA(jj))
      ENDDO
      END
