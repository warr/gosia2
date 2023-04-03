
C----------------------------------------------------------------------
C SUBROUTINE RANGEL
C
C Called by: ALLOC
C
C Purpose: to determine the range of the integration over omega.
C
C Uses global variables:
C      ACC50  - accuracy required for integration
C      IRA    - range for omega for each multipolarity
C      MULTI  - number of matrix elements with each multipolarity populating levels
C
C Formal parameters:
C      Acc1   - the desired accuracy
C
C \omega_max >= \alpha_\lambda - {1 \over \lambda} \ln(a_c)
C where a_c is Acc1 here.
C
C The gosia documentation gives a table for \alpha_\lambda: E1 = -0.693,
C E2 = 0.203, E3 = 0.536, E4 = 0.716, E5 = 0.829, E6 = 0.962, M1 = 0.203,
C M2 = 0.536.
C
C Note that first we work out omega, but then we work out the appropriate
C index, knowing that we are always using steps of 0.03.


      SUBROUTINE RANGEL(Acc1)
      IMPLICIT NONE
      REAL*8 Acc1 , ACC50 , acl , w
      INTEGER*4 i , IRA , LAMDA , LAMMAX , LDNUM , LEAD , MAXLA , MULTI
      COMMON /A50   / ACC50
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /CLCOM / LAMDA(8) , LEAD(2,1500) , LDNUM(8,75) , LAMMAX ,
     &                MULTI(8)

      acl = -LOG(Acc1)
      ACC50 = Acc1/50.
      DO i = 1 , 8 ! Loop over multipolarity 1..6 = E1..6, 7,8 = M1,M2
         IF ( MULTI(i).NE.0 ) THEN
            IF ( i.EQ.2 .OR. i.EQ.7 ) THEN ! E2 or M1
               w = acl/2. + .203
            ELSEIF ( i.EQ.3 .OR. i.EQ.8 ) THEN ! E3 or M2
               w = acl/3. + .536
            ELSEIF ( i.EQ.4 ) THEN ! E4
               w = acl/4. + .716
            ELSEIF ( i.EQ.5 ) THEN ! E5
               w = acl/5. + .829
            ELSEIF ( i.EQ.6 ) THEN ! E6
               w = acl/6. + .962
            ELSE
               w = acl - .693 ! E1
            ENDIF
            w = w/.03        ! We step in steps of \Delta\omega = 0.03
            IRA(i) = INT(w+1.5)
         ELSE
            IRA(i) = 0
         ENDIF
      ENDDO
      IF ( IRA(7).NE.0 ) IRA(7) = IRA(7) + 1
      IF ( IRA(8).NE.0 ) IRA(8) = IRA(8) + 1
      END
