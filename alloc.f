      SUBROUTINE ALLOC(Accur)
      IMPLICIT NONE
      REAL*8 Accur , u , v
      INTEGER*4 iflag , IRA , j , k , k1 , load , LOCQ , LP1 , LP10 , 
     &          LP11 , LP12 , LP13 , LP14 , LP2 , LP3 , LP4 , LP6 , 
     &          LP7 , LP8 , LP9
      INTEGER*4 MAXLA
      COMMON /ALLC  / LOCQ(8,7)
      COMMON /RNG   / IRA(8) , MAXLA
      COMMON /MGN   / LP1 , LP2 , LP3 , LP4 , LP6 , LP7 , LP8 , LP9 , 
     &                LP10 , LP11 , LP12 , LP13 , LP14
      CALL RANGEL(Accur)
      load = 0
      iflag = 0
      DO j = 1 , 8
         DO k = 1 , 7
            LOCQ(j,k) = 0
         ENDDO
      ENDDO
      DO k = 1 , 6
         k1 = k + 1
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k)
         ENDDO
      ENDDO
      DO k = 7 , 8
         k1 = k - 6
         DO j = 1 , k1
            LOCQ(k,j) = load
            load = load + IRA(k)
         ENDDO
      ENDDO
      IF ( load.LE.LP14 ) RETURN
      WRITE (22,99001)
99001 FORMAT (5X,'NO SPACE FOR Q FUNCTIONS TABULATION'//5X,
     &        'SORRY,JOB WILL BE BRUTALLY TERMINATED!')
      v = -1.
      u = LOG10(v)
      u = SIN(u)
      END
