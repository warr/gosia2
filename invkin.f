C----------------------------------------------------------------------
C SUBROUTINE INVKIN
C
C Called by: GOSIA
C
C Purpose: calculate the angle of the scattered projectile in the lab frame
C          when the user gave the angle of the recoiling target nucleus in
C          the lab frame. There are two solutions to this problem, so Iflag
C          = 1 selects the larger angle one and Iflag = 2 the smaller one.
C          Note that the smaller angle (Iflag = 2) corresponds to very low
C          energies of the recoiling target nucleus, which probably either
C          don't get out of the target or don't get detected. So Iflag = 2
C          is probably not very useful! Also, this routine calculates the
C          correct value of the kinematic flag IKIN.
C
C Formal parameters:
C      E_p     - Beam energy in MeV (readonly)
C      E_x     - energy of excited state to use for kinematic in MeV (readonly)
C      I_Z     - Projectile/target flag. -ve if projectile excitation
C      M_inv   - mass of investigated nuclei in AMU (readonly)
C      M_non   - mass of non-investigated nuclei in AMU (readonly)
C      Theta_t - theta of recoiling target nucleus in lab frame (readonly)
C      Theta_p - theta of scattered projectile in lab frame (writeonly)
C      Iflag   - flag to select one of two possible solutions (readonly)
C      Ikin    - kinematic flag (writeonly)
      
      SUBROUTINE INVKIN(E_p, E_x , I_z , M_inv , M_non , Theta_t ,
     &                  Theta_p , Iflag , Ikin)

      REAL*8 E_p , M_inv , M_non , Theta_t , Theta_p , E_x , M_p , M_t
      REAL*8 ared , epmin , t , x(2), y , thres , tau , taup
      INTEGER*4 Iflag , Ikin , I_z

C     Sort out which is the projectile and which is the target
      
      IF ( I_z.LT.0 ) THEN
         M_p = M_inv ! Projectile is investigated
         M_t = M_non ! Target is non investigated
      ELSE
         M_p = M_non ! Projectile is non investigated
         M_t = M_inv ! Target is investigated
      ENDIF

C     Reduced mass
      
      ared = 1 + M_p / M_t
      
C     Excitation energy of inelastically scattered particle when state at
C     energy E_x is excited
      
      epmin = E_p - E_x * ared
      
C     Tau
      
      taup = sqrt(E_p / epmin)
      tau = taup * M_p / M_t
      
C     Calculate the two solutions
      
      y = tan(theta_t/57.2957795)
      t = taup * taup * y * y * y * y -
     &      (1 + y * y) * (taup * taup * y * y - 1)
      t = sqrt(t)
      x(1) = (taup * y * y + t) / (1 + y * y)
      x(1) = atan2(sqrt(1 - x(1) * x(1)), tau + x(1))
      x(2) = (taup * y * y - t) / (1 + y * y)
      x(2) = atan2(sqrt(1 - x(2) * x(2)), tau + x(2))
      
C     Select the solution we want according to the flag. Note that the
C     solution with the lower angle corresponds to target recoils which
C     are probably undetectable.
      
      IF ( Iflag.EQ.1 ) THEN
         theta_p = MAX(x(1),x(2))*57.2957795
      ELSE
         theta_p = MIN(x(1),x(2))*57.2957795
      ENDIF
      
C     Calculate angle of scattered projectile in centre of mass frame, for
C     which the maximum laboratory scattering angle is reached.
      
      t = acos(-1./tau)
      
C     Now calculate the arctangent of the corresponding angle for the
C     recoiling target nuclei in the laboratory frame
      
      thres = sin(t)/(taup-cos(t))
      
C     So now, if y = tan(theta_t_lab) > thres, we are above the maximum and
C     need the larger value of theta_p_cm, so we set Ikin to 1. Otherwise we
C     are below the maximum and need the smaller value so we choose Ikin = 0.
      
      IF ( y.GT.thres ) THEN
         Ikin = 1
      ELSE
         Ikin = 0
      ENDIF
      
      END
