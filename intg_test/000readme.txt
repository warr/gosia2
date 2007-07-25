This test is to test the OP,INTG option.

The experiment is a 92Kr radioactive beam on a 109Ag target measured with
Miniball at REX-ISOLDE.

The files are:
ag109_intg.inp      - input file for 109Ag
kr92_intg.inp       - input file for 92Kr
ag109.yie           - experimental yields for 109Ag
kr92.yie            - experimental yields for 92Kr
miniball_gdet.f8    - miniball detector configuration (unit 8)
miniball_gdet.f9    - miniball detector configuration (unit 9)
ag109_intg.out.good - expected output for 109Ag (verified against gosia1)
kr92_intg.out.good  - expected output for 92Kr (verified against gosia1)

To run this example: make
   This should generate the files ag109_intg.out and kr92_intg.out, which
can be compared to those with the .good suffix.

If you want to try it with gosia1, you need to remove the first line of the
two .inp files and run gosia once on each file.
