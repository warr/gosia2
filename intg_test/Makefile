INP += ag109_intg.inp
INP += kr92_intg.inp
INP += ag109.yie
INP += kr92.yie
INP += miniball_gdet.f8
INP += miniball_gdet.f9
INP += ag109_intg.out.good
INP += kr92_intg.out.good

OUT += ag109.cor
OUT += kr92.cor
OUT += ag109.me
OUT += kr92.me
OUT += ag109_intg.out
OUT += kr92_intg.out

all: compare

ag109_intg.inp: unpack

kr92_intg.inp: unpack

ag109.yie: unpack

kr92.yie: unpack

miniball_gdet.f8: unpack

miniball_gdet.f9: unpack

unpack: intg_test.tgz
	tar zxvf intg_test.tgz

../gosia2:
	make -C ..

ag109_intg.out: ../gosia2 ag109_intg.inp
	../gosia2 < ag109_intg.inp

kr92_intg.out: ../gosia2 kr92_intg.inp
	../gosia2 < kr92_intg.inp

compare: ag109_intg.out kr92_intg.out
	diff ag109_intg.out ag109_intg.out.good
	diff kr92_intg.out kr92_intg.out.good

clean:
	rm -f *~ fort.* $(INP) $(OUT)

tidy:
	rm -f *~ fort.*

