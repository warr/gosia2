C FUNCTION CCLKUP
C
C Called by: BRICC
C Calls:     SPLNER
C
C Purpose: looks up a single conversion coefficient for given Z, energy and
C          multipolarity in the BrIcc database, interpolating for the right
C          energy using a cubic spline on a log-log scale.
C
C Formal parameters:
C      Myz     - Z of nucleus to look up
C      Egamma  - Gamma ray energy for which we want to interpolate
C      Myimult - Multipolarity (1=E1,2=E2,3=E3,6=M1,7=M2)
      REAL*8 FUNCTION CCLKUP(Myz,Egamma,Myimult)
 
      IMPLICIT NONE
 
      INTEGER*4 iz , ia , flag(37) , nrec(37) , rec1(37) , Myz
      INTEGER*4 ishell , irec , ienergy , imult , Myimult
      REAL*8 x(500) , y(500) , result , Egamma
      REAL*4 binding_energy(500) , energy(37,500) , cc(37,500,10)
      CHARACTER*8 file(37) , name(37)
      CHARACTER*4 element
 
C     Initialise
      CCLKUP = 0.0D0

C     We can't calculate above 6000 keV
      IF ( Egamma.GE.6000 ) RETURN
  
C     Read the index record for that Z
      READ (30,REC=Myz,ERR=100) iz , element , ia , 
     &                          (flag(ishell),binding_energy(ishell),
     &                          file(ishell),name(ishell),nrec(ishell),
     &                          rec1(ishell),ishell=1,37)

C     Ignore internal pair conversion below 1100 keV
      IF ( Egamma.LE.1100 ) flag(37) = 0
      
C     Now read the internal conversion coefficient data records for each shell
      DO ishell = 1 , 37 ! Only first 37 subshells

C        Flag that shell doesn't contribute if we are below binding energy
         IF ( Egamma.LE.binding_energy(ishell) ) flag(ishell) = 0

         IF ( flag(ishell).EQ.-1 ) THEN ! If subshell is present in record

C           Read all the records for this subshell
            ienergy = 1
            DO irec = rec1(ishell) , rec1(ishell) + nrec(ishell) - 1
               READ (31,REC=irec,ERR=100) energy(ishell,ienergy) , 
     &               (cc(ishell,ienergy,imult),imult=1,10)
               ienergy = ienergy + 1
            ENDDO ! Loop on records
             
            IF ( cc(ishell,1,Myimult).EQ.0.0D0 ) flag(ishell) = 0
          ENDIF ! If subshell is present in record
      ENDDO ! Loop on subshells
 
C     Interpolate for each subshell
      DO ishell = 1 , 37 ! Only first 37 subshells                       

        IF ( flag(ishell).EQ.-1 ) THEN ! If subshell is present in record
 
C           If the energy is less than the first data point use its ICC, or
C           if it is more than the last data point, otherwise we interpolate
          IF ( Egamma.LE.energy(ishell,1) ) THEN
               result = DBLE(CC(ishell,1,Myimult))
               WRITE (22,'(A,F7.4,3A)') 'Warning Egamma=',Egamma/1.D3,
     &           ' is in regime where solid state effects dominate',
     &           ' conversion coefficients for shell ',name(ishell)
            ELSEIF ( Egamma .GE. energy(ishell,nrec(ishell)) ) THEN
               result = DBLE(cc(ishell,nrec(ishell),Myimult))
               WRITE (22,'(A,F7.4,3A)') 'Warning Egamma=',Egamma/1.D3,
     &           ' exceeds range of conversion coefficients table',
     &           ' for shell ',name(ishell)
            ELSE
C             Set up for interpolation
              DO ienergy = 1 , nrec(ishell)
                 x(ienergy) = LOG(DBLE(energy(ishell,ienergy)))
                 y(ienergy) = DBLE(cc(ishell,ienergy,Myimult)) ! Log for this is done in FUNC
              ENDDO
 
C             Perform spline over data
              CALL SPLNER(x,y,nrec(ishell),LOG(Egamma),result,2)
            ENDIF

C           Add the conversion coefficients of each subshell
            CCLKUP = CCLKUP + result
            
         ENDIF ! If subshell is present in record
      ENDDO ! Loop on subshells
 
      RETURN
 
 100  WRITE (*,*) 'ERROR - No data found for this Z ', Myz
      CCLKUP = -1.0D0
      RETURN
      END
 
