 
C----------------------------------------------------------------------
C SUBROUTINE NEWLV
C
C Called by: AMPDER, STING
C Calls:     EXPON, LEADF, MEM
C
C Purpose: Setup a new level which can be excited from ground state. We store
C       ISSTAR, ISSTO and MSTORE for the level and calculate and store the
C       exponential: exp(i \xi_{kn} (\epsilon \sinh(\omega) + \omega))
C
C Uses global variables:
C      EXPO   - adiabatic exponential
C      IFLG   - flag to determine whether to calculate exponential (so we don't calculate twice)
C      ISG    - index for sigma
C      ISG1   - 
C      ISSTAR - index of last substate for that level
C      ISSTO  - index of first substate for that level
C      KDIV   - index for division
C      LDNUM  - number of matrix elements with each multipolarity populating level
C      MSTORE - index of final level number and index of matrix element
C      NDIV   - number of divisions
C      NPT    - index in ADB array (this is omega / 0.03)
C      NSTART - index in CAT of first substate associated with a level
C      NSTOP  - index in CAT of last substate associated with a level
C
C Formal parameters:
C      N      - level number
C      Ld     - Number of matrix elements for level N multipolarity La
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

      Ld = LDNUM(La,N) ! Get number of levels connected to level N by multipolarity La
      IF ( Ld.EQ.0 ) RETURN ! Return if there aren't any

      DO i2 = 1 , Ld ! For each level
         m = LEADF(N,i2,La) ! Get the other level associated
         ISSTAR(i2) = NSTOP(m) ! Get the index of last substate for that level
         ISSTO(i2) = NSTART(m) ! Get the index of first substate for that level
         MSTORE(1,i2) = m ! Store the final level number
         indx = MEM(N,m,La) ! Index for matrix element from level N to level m with multipolarity La
         MSTORE(2,i2) = indx ! Store index of matrix element
         IF ( IFLG.NE.0 ) THEN
            IF ( m.NE.N ) EXPO(indx) = EXPON(indx,NPT,ISG,ISG1,NDIV,KDIV
     &                                 )
         ENDIF
      ENDDO
      END
