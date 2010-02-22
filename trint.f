 
C----------------------------------------------------------------------
C SUBROUTINE TRINT
C
C Called by: STAMP
C Calls:     POL4
C
C Purpose: calculate sine and cosine integrals (Si and Ci). Note, that we
C actually calculate pi/2-Si and -Ci. Note also that a constant added to these
C values doesn't make any difference, because STAMP always subtracts them
C pairwise from each other.
C
C Formal parameters:
C      Arg    - value of x for which to evaluate sine and cosine integrals
C      Si     - returned sine integral at that value (actually pi/2 - Si)
C      Ci     - returned cosine integral at that value (actually -Ci)
C
C For small x we use the series expansion. See Abramowitz and Stegun Handbook
C of Mathematical Functions with Formulas, Graphs and Mathematical Tables,
C National Bureau of Standards, 8th Ed. P232 Eqs. 5.2.14 and 5.2.16, except we
C calculate pi/2 - Si and -Ci:
C
C pi/2 - Si = pi/2 - x + x^3 / (3! * 3) - x^5 / (5! * 5) + x^7 / (7! * 7) + ...
C pi/2                       = 1.57079632679
C 1 / (3! * 3) = 1 / 18      = 0.05555555
C 1 / (5! * 5) = 1 / 600     = 1.666667E-3
C 1 / (7! * 7) = 1 / 35280   = 2.83446E-5
C
C -Ci = -Gamma - ln(x) + x^2 / (2! * 2) - x^4 / (4! * 4) + x^6 / (6! * 6) - ...
C Gamma        = Euler Gamma = 0.577215664902
C 1 / (2! * 2) = 1 / 4       = 0.25
C 1 / (4! * 4) = 1 / 96      = 0.0104166
C 1 / (6! * 6) = 1 / 4320    = 2.31481E-4
C 1 / (8! * 8) = 1 / 322560  = 3.10019E-6
C
C For large x we use the rational approximations. See Abramowitz and Stegun
C Handbook of Mathematical Functions with Formulas, Graphs and Mathematical
C Tables, National Bureau of Standards, 8th Ed. P233 Eqs. 5.2.38 and 5.2.39 to
C calculate the auxillary functions f and g and then use 5.2.8 and 5.2.9 to
C obtain the values of pi/2 - Si and -Ci from f and g.
      
      SUBROUTINE TRINT(Arg,Si,Ci)
      IMPLICIT NONE
      REAL*8 a , Arg , c , Ci , f , g , POL4 , s , Si

      a = Arg*Arg

C     If Arg is small, use the polynomial expansion. The coefficients are
C     evaluated from Abramowitz and Stegun 5.2.14 and 5.2.16 as shown above:
      IF ( Arg.LT.1. ) THEN
         Si = POL4(0.D0,2.83446712D-5,-1.66666667D-3,.055555555D0,-1.D0,
     &        a)
         Si = Si*Arg
         Si = Si + 1.57079632679D0 ! This is actually pi/2 - Si
         Ci = POL4(-3.100198413D-6,2.314814815D-4,-.0104166667D0,.25D0,
     &        0.D0,a)
         Ci = Ci - LOG(Arg) - 0.577215664902D0 ! This is actually -Ci
         GOTO 99999
      ENDIF

C     Otherwise use the expansion in terms of sine and cosine
      s = SIN(Arg)
      c = COS(Arg)

C     Here we use an approximation. If Arg is quite large, a is very large 
C     and the four polynomials are all huge. Moreover, the four polynomials 
C     are almost identical, so the ratios are unity. So in this case, 
C     f = 1./Arg and g=1./a is a good approximation.
      
      f = 1.
      g = 1.

C     From Abramowitz and Stegun 5.2.38 and 5.2.39 we have the following
C     relations for the auxillary functions f and g, using the coefficients
C     from that reference:
      IF ( a.LE.1.D+8 ) THEN
         f = POL4(1.D0,38.027246D0,265.187033D0,335.67732D0,38.102495D0,
     &       a)
         f = f/POL4(1.D0,40.021433D0,322.624911D0,570.23628D0,
     &       157.105423D0,a)
         g = POL4(1.D0,42.242855D0,302.757865D0,352.018498D0,
     &       21.821899D0,a)
         g = g/POL4(1.D0,48.196927D0,482.485984D0,1114.978885D0,
     &       449.690326D0,a)
      ENDIF
      
      f = f/Arg
      g = g/a

C     From Abramowitz and Stegun 5.2.8 and 5.2.9:      
      Si = f*c + g*s ! This is actually pi/2 - Si compared to Abramowitz and Stegun
      Ci = g*c - f*s ! This is actually -Ci compared to Abramowitz and Stegun
99999 END
