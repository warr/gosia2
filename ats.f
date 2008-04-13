      FUNCTION ATS(N)
      IMPLICIT NONE
      REAL*8 ATS , x , xm
      INTEGER*4 m , N
      IF ( N.LE.0 .OR. N.GT.96 ) THEN
         ATS = 0.
         RETURN
      ELSE
         x = N/2. + 1
         m = N/2 + 1
         xm = DBLE(m)
         IF ( ABS(x-xm).GE.1.E-9 ) THEN
            IF ( m.EQ.1 .OR. m.EQ.2 .OR. m.EQ.3 .OR. m.EQ.6 .OR. 
     &           m.EQ.7 .OR. m.EQ.10 .OR. m.EQ.15 .OR. m.EQ.16 .OR. 
     &           m.EQ.19 .OR. m.EQ.24 .OR. m.EQ.25 .OR. m.EQ.28 .OR. 
     &           m.EQ.31 .OR. m.EQ.35 .OR. m.EQ.37 .OR. m.EQ.40 .OR. 
     &           m.EQ.41 .OR. m.EQ.44 ) THEN
               ATS = .5
               RETURN
            ELSEIF ( m.EQ.4 .OR. m.EQ.5 .OR. m.EQ.8 .OR. m.EQ.9 .OR. 
     &               m.EQ.11 .OR. m.EQ.17 .OR. m.EQ.18 .OR. m.EQ.20 .OR. 
     &               m.EQ.26 .OR. m.EQ.27 .OR. m.EQ.36 .OR. m.EQ.42 .OR. 
     &               m.EQ.43 .OR. m.EQ.45 ) THEN
               ATS = 1.5
               RETURN
            ELSEIF ( m.EQ.12 .OR. m.EQ.14 .OR. m.EQ.21 .OR. m.EQ.23 .OR. 
     &               m.EQ.32 .OR. m.EQ.39 ) THEN
               ATS = 2.5
               RETURN
            ELSEIF ( m.EQ.13 .OR. m.EQ.22 .OR. m.EQ.38 ) THEN
               ATS = 4.5
               RETURN
            ELSEIF ( m.EQ.29 .OR. m.EQ.30 .OR. m.EQ.48 ) THEN
               ATS = 3.5
               RETURN
            ELSEIF ( m.EQ.33 ) THEN
               ATS = 7.5
               RETURN
            ELSEIF ( m.EQ.34 ) THEN
               ATS = 6.5
               GOTO 99999
            ELSEIF ( m.EQ.46 .OR. m.EQ.47 ) THEN
               ATS = 5.5
               RETURN
            ENDIF
         ENDIF
         m = m - 1
         IF ( m.EQ.4 .OR. m.EQ.8 .OR. m.EQ.17 .OR. m.EQ.26 .OR. 
     &        m.EQ.28 .OR. m.EQ.30 .OR. m.EQ.32 .OR. m.EQ.42 .OR. 
     &        m.EQ.45 .OR. m.EQ.48 ) THEN
            ATS = 2.
            RETURN
         ELSEIF ( m.EQ.10 .OR. m.EQ.36 ) THEN
         ELSEIF ( m.EQ.12 .OR. m.EQ.21 .OR. m.EQ.37 ) THEN
            ATS = 3.
            RETURN
         ELSEIF ( m.EQ.13 .OR. m.EQ.22 .OR. m.EQ.29 .OR. m.EQ.31 .OR. 
     &            m.EQ.34 .OR. m.EQ.38 .OR. m.EQ.47 ) THEN
            ATS = 4.
            RETURN
         ELSEIF ( m.EQ.33 ) THEN
            ATS = 8.
            RETURN
         ELSEIF ( m.EQ.46 ) THEN
            ATS = 6.
            RETURN
         ELSE
            ATS = 0.
            RETURN
         ENDIF
      ENDIF
      ATS = 1.
99999 END
