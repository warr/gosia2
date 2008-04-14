 
C----------------------------------------------------------------------
C SUBROUTINE NEWLV
C
C Called by: AMPDER, STING
C Calls:     EXPON, LEADF, MEM
C
C Purpose: calculate and store the exponential:
C       exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C
C Uses global variables:
C      EXPO   - adiabatic exponential
C      IFAC   -
C      IFLG   - flag to determine whether to calculate exponential (so we don't calculate twice)
C      ISG    -
C      ISSTAR -
C      ISSTO  -
C      KDIV   -
C      LDNUM  - number of matrix elements with each multipolarity populating level
C      MSTORE -
C      NDIV   -
C      NPT    -
C      NSTART -
C      NSTOP  -
C
C Formal parameters:
C      N      - level number
C      Ld     -
C      La     - multipolarity
C
C Note that the exponential is calculated by EXPON. This file does the
C storage part.
      
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
