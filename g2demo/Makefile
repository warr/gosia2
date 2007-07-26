GOSIA2DIR=/home/warr/devel/gosia2

INP += g2demo-beam.inp
INP += g2demo-targ.inp
INP += g2demo-beam.yld
INP += g2demo-targ.yld
INP += det.gdt
INP += g2demo-beam.cor
INP += g2demo-targ.cor
INP += g2demo-beam.map
INP += g2demo-targ.map
INP += g2demo-beam.out.good
INP += g2demo-targ.out.good

OUT += cnor.dat
OUT += crf.dat
OUT += dum.14
OUT += dum.15
OUT += dum.17
OUT += dum.18
OUT += g2demo-beam.out
OUT += g2demo-targ.out
OUT += g2demo-beam.bst
OUT += g2demo-beam.smr
OUT += g2demo-targ.bst
OUT += g2demo-targ.smr

all: compare

g2demo-beam.inp: unpack

g2demo-targ.inp: unpack

unpack: g2demo.tgz
	tar zxvf g2demo.tgz

g2demo-beam.out: g2demo-beam.inp
	$(GOSIA2DIR)/gosia2 < g2demo-beam.inp

g2demo-targ.out: g2demo-targ.inp
	$(GOSIA2DIR)/gosia2 < g2demo-targ.inp

compare: g2demo-beam.out g2demo-targ.out
	diff g2demo-beam.out*
	diff g2demo-targ.out*

clean:
	rm -f *~ $(INP) $(OUT)

