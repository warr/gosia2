      SUBROUTINE NEWLV(N,Ld,La)
      IMPLICIT NONE
      REAL*8 D2W
      INTEGER*4 i2 , IFLG , indx , ISG , ISG1 , ISSTAR , ISSTO , KDIV , 
     &          La , LAMDA , LAMMAX , LAMR , Ld , LDNUM , LEAD , LEADF , 
     &          m , MEM , MSTORE , MULTI
      INTEGER*4 N , NDIV , NPT , NSTART , NSTOP , NSW
      COMPLEX*16 EXPO , EXPON
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      COMMON /CAUX  / NPT , NDIV , KDIV , LAMR(8) , ISG , D2W , NSW , 
     &                ISG1
      COMMON /PINT  / ISSTAR(76) , ISSTO(75) , MSTORE(2,75)
      COMMON /ADBXI / EXPO(500)
      COMMON /FLA   / IFLG
      COMMON /CEXC0 / NSTART(76) , NSTOP(75)
      Ld = LDNUM(La,N)
      IF ( Ld.EQ.0 ) RETURN
      DO i2 = 1 , Ld
         m = LEADF(N,i2,La)
         ISSTAR(i2) = NSTOP(m)
         ISSTO(i2) = NSTART(m)
         MSTORE(1,i2) = m
         indx = MEM(N,m,La)
         MSTORE(2,i2) = indx
         IF ( IFLG.NE.0 ) THEN
            IF ( m.NE.N ) EXPO(indx) = EXPON(indx,NPT,ISG,ISG1,NDIV,KDIV
     &                                 )
         ENDIF
      ENDDO
      END
