GOSIA2DIR=/home/warr/devel/gosia2

INP += ag109_mini0.inp
INP += kr92_mini0.inp
INP += ag109.yie
INP += kr92.yie
INP += ag109.cor
INP += kr92.cor
INP += ag109.map
INP += kr92.map
INP += ag109_mini0.out.good
INP += kr92_mini0.out.good
INP += miniball_gdet.f8
INP += miniball_gdet.f9

OUT += ag109_mini0.out
OUT += kr92_mini0.out
OUT += ag109.me
OUT += kr92.me
OUT += cnor.dat

all: ag109_mini0.out kr92_mini0.out

$(GOSIA2DIR)/gosia2:
	make -C $(GOSIA2DIR)

ag109_mini0.inp: unpack

kr92_mini0.inp: unpack

unpack:
	tar zxvf mini_test.tgz

ag109_mini0.out: $(GOSIA2DIR)/gosia2 ag109_mini0.inp
	$(GOSIA2DIR)/gosia2 < ag109_mini0.inp

kr92_mini0.out: $(GOSIA2DIR)/gosia2 kr92_mini0.inp
	$(GOSIA2DIR)/gosia2 < kr92_mini0.inp


clean:
	rm -f *~ fort.* $(INP) $(OUT)

tidy:
	rm -f *~ fort.*

