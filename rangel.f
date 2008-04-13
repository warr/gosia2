      SUBROUTINE RANGEL(Acc1)
      IMPLICIT NONE
      REAL*8 Acc1 , ACC50 , acl , w
      INTEGER*4 i , IRA , LAMDA , LAMMAX , LDNUM , LEAD , MAXLA , MULTI
      COMMON /A50   / ACC50
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /CLCOM / LAMDA(8) , LEAD(2,500) , LDNUM(8,75) , LAMMAX , 
     &                MULTI(8)
      acl = -LOG(Acc1)
      ACC50 = Acc1/50.
      DO i = 1 , 8
         IF ( MULTI(i).NE.0 ) THEN
            IF ( i.EQ.2 .OR. i.EQ.7 ) THEN
               w = acl/2. + .203
            ELSEIF ( i.EQ.3 .OR. i.EQ.8 ) THEN
               w = acl/3. + .536
            ELSEIF ( i.EQ.4 ) THEN
               w = acl/4. + .716
            ELSEIF ( i.EQ.5 ) THEN
               w = acl/5. + .829
            ELSEIF ( i.EQ.6 ) THEN
               w = acl/6. + .962
            ELSE
               w = acl - .693
            ENDIF
            w = w/.03
            IRA(i) = INT(w+1.5)
         ELSE
            IRA(i) = 0
         ENDIF
      ENDDO
      IF ( IRA(7).NE.0 ) IRA(7) = IRA(7) + 1
      IF ( IRA(8).NE.0 ) IRA(8) = IRA(8) + 1
      END
