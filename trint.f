 
C----------------------------------------------------------------------
C SUBROUTINE TRINT
C
C Called by: STAMP
C Calls:     POL4
C
C Purpose: calculate sine and cosine integrals (Si and Ci)
C
C Formal parameters:
C      Arg    - value of x for which to evaluate sine and cosine integrals
C      Si     - returned sine integral at that value
C      Ci     - returned cosine integral at that value
C
C For small x:
C Si = x - x^3 / (3! * 3) + x^5 / (5! * 5) - x^7 / (7! * 7) + ...
C 1 / (3! * 3) = 1 / 18    = 0.05555555
C 1 / (5! * 5) = 1 / 600   = 1.666667E-3
C 1 / (7! * 7) = 1 / 35280 = 2.83446E-5
C
C Ci = -gamma - ln(x) + x^2 / (2! * 2) - x^4 / (4! * 4) + x^6 / (6! * 6) - ...
C where gamma is the Euler gamma = 0.5772156649
C 1 / (2! * 2) = 1 / 4      = 0.25
C 1 / (4! * 4) = 1 / 96     = 0.0104166
C 1 / (6! * 6) = 1 / 4320   = 2.31481E-4
C 1 / (8! * 8) = 1 / 322560 = 3.10019E-6
C
C For large x:
C Si = pi / 2 - cos(x) / x * (1 - 2! / x^2 + 4! / x^4 - ...)
C             - sin(x) / x * (1 / x - 3! / x^3 + 5! / x^5 - ...)
C
C Ci = cos(x) / x * (1 / x - 3! / x^3 + 5! / x^5 - ...)
C      - sin(x) / x * (1 - 2! / x^2 + 4! / x^4 - ...)
C
C Note that this function seems to differ a little from the standard version
C described above.
      
      SUBROUTINE TRINT(Arg,Si,Ci)
      IMPLICIT NONE
      REAL*8 a , Arg , c , Ci , f , g , POL4 , s , Si

      a = Arg*Arg

C     If Arg is small, use the polynomial expansion
      IF ( Arg.LT.1. ) THEN
         Si = POL4(0.D0,2.83446712D-5,-1.66666667D-3,.055555555D0,-1.D0,
     &        a)
         Si = Si*Arg
         Ci = POL4(-3.100198413D-6,2.314814815D-4,-.0104166667D0,.25D0,
     &        0.D0,a)
         Ci = Ci - LOG(Arg)
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

C See Abramowitz and Segun - Handbook of Mathematical Functions with Formulas, Graphs
C and Mathematical Tables, National Bureau of Standards, 8th Ed. P 233 for the following
C coefficients.
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
      
      Si = f*c + g*s
      Ci = g*c - f*s
99999 END
