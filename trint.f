      SUBROUTINE TRINT(Arg,Si,Ci)
      IMPLICIT NONE
      REAL*8 a , Arg , c , Ci , f , g , POL4 , s , Si
      a = Arg*Arg
      IF ( Arg.LT.1. ) THEN
         Si = POL4(0.D0,2.83446712D-5,-1.66666667D-3,.055555555D0,-1.D0,
     &        a)
         Si = Si*Arg
         Ci = POL4(-3.100198413D-6,2.314814815D-4,-.0104166667D0,.25D0,
     &        0.D0,a)
         Ci = Ci - LOG(Arg)
         GOTO 99999
      ENDIF
      s = SIN(Arg)
      c = COS(Arg)
      f = POL4(1.D0,38.027246D0,265.187033D0,335.67732D0,38.102495D0,a)
      f = f/POL4(1.D0,40.021433D0,322.624911D0,570.23628D0,157.105423D0,
     &    a)/Arg
      g = POL4(1.D0,42.242855D0,302.757865D0,352.018498D0,21.821899D0,a)
      g = g/POL4(1.D0,48.196927D0,482.485984D0,1114.978885D0,
     &    449.690326D0,a)/a
      Si = f*c + g*s
      Ci = g*c - f*s
99999 END
