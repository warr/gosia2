BINDIR=$(ROOT)/usr/bin
MANDIR=$(ROOT)/usr/share/man/man1

UID=0
GID=0

EXE=gosia2
MAN=gosia2.1

FC=g77

# Turn on debugging - note that -Wall and -O2 together gives warnings
# about variables being possibly used without being initialised. These
# cases seem to be harmless. i.e. if (condition) then a=1 else a=2 endif
# and then using a gives this warning, but in fact a is always set to
# something
#FFLAGS = -g -Wall

# Turn on optimisation
FFLAGS += -O -funroll-loops

DEPS=Makefile

ALL: $(EXE)

OBJS += adhoc.o
OBJS += alloc.o
OBJS += ampder.o
OBJS += angula.o
OBJS += apram.o
OBJS += arccos.o
OBJS += arctg.o
OBJS += ats.o
OBJS += branr.o
OBJS += cclkup.o
OBJS += cegry.o
OBJS += chmem.o
OBJS += cmlab.o
OBJS += code7.o
OBJS += conv.o
OBJS += coord.o
OBJS += decay.o
OBJS += djmm.o
OBJS += double.o
OBJS += effix.o
OBJS += elmt.o
OBJS += expon.o
OBJS += f.o
OBJS += fakp.o
OBJS += faza.o
OBJS += faza1.o
OBJS += fhip.o
OBJS += fiint.o
OBJS += fiint1.o
OBJS += ftbm.o
OBJS += func.o
OBJS += func1.o
OBJS += fxis1.o
OBJS += fxis2.o
OBJS += gamatt.o
OBJS += gcf.o
OBJS += gf.o
OBJS += gkk.o
OBJS += gkvac.o
OBJS += gosia2.o
OBJS += half.o
OBJS += intg.o
OBJS += klopot.o
OBJS += kontur.o
OBJS += lagran.o
OBJS += laiamp.o
OBJS += laisum.o
OBJS += leadf.o
OBJS += limits.o
OBJS += load.o
OBJS += lsloop.o
OBJS += mem.o
OBJS += mini.o
OBJS += mixr.o
OBJS += mixup.o
OBJS += newcat.o
OBJS += newcnv.o
OBJS += newlv.o
OBJS += openf.o
OBJS += path.o
OBJS += podziel.o
OBJS += pol4.o
OBJS += pomnoz.o
OBJS += prelm.o
OBJS += prim.o
OBJS += pticc.o
OBJS += qe.o
OBJS += qfit.o
OBJS += qm.o
OBJS += qrange.o
OBJS += rangel.o
OBJS += ready.o
OBJS += recoil.o
OBJS += reset.o
OBJS += rk4.o
OBJS += rndm.o
OBJS += rotate.o
OBJS += seq.o
OBJS += setin.o
OBJS += simin.o
OBJS += sixel.o
OBJS += snake.o
OBJS += splner.o
OBJS += spline.o
OBJS += splint.o
OBJS += stamp.o
OBJS += sting.o
OBJS += szereg.o
OBJS += tacos.o
OBJS += tapma.o
OBJS += tasin.o
OBJS += tcabs.o
OBJS += tcexp.o
OBJS += tenb.o
OBJS += tens.o
OBJS += trint.o
OBJS += wsixj.o
OBJS += wthrej.o
OBJS += xstatic.o
OBJS += ylm.o
OBJS += ylm1.o

SINGLE_FILE = gosia2_single_file.f

gosia2: $(OBJS) $(DEPS)
	$(FC) $(LDFLAGS) -o gosia2 $(OBJS)

clean:
	rm -f *~ *.o $(EXE) $(SINGLE_FILE)

install: $(EXE) $(MAN)
	install -m 755 -o $(UID) -g $(GID) -d $(BINDIR)
	install -m 755 -o $(UID) -g $(GID) -s $(EXE) $(BINDIR)
	install -m 755 -o $(UID) -g $(GID) -d $(MANDIR)
	install -m 644 -o $(UID) -g $(GID) $(MAN) $(MANDIR)
	gzip -f $(MANDIR)/$(MAN)

single_file:
	cat gosia2.f arccos.f arctg.f load.f lsloop.f leadf.f \
	mem.f cmlab.f qe.f qm.f snake.f fhip.f alloc.f \
	rangel.f qrange.f ampder.f laisum.f expon.f faza.f \
	setin.f sting.f laiamp.f faza1.f trint.f pol4.f \
	stamp.f reset.f half.f double.f path.f intg.f newlv.f \
	code7.f apram.f newcat.f pomnoz.f tenb.f tens.f djmm.f \
	ftbm.f mini.f cegry.f fakp.f prim.f seq.f gf.f f.f \
	conv.f wthrej.f wsixj.f lagran.f func.f func1.f \
	gkvac.f gkk.f xstatic.f ats.f ylm.f decay.f angula.f \
	ready.f branr.f limits.f szereg.f sixel.f prelm.f \
	recoil.f rotate.f ylm1.f fiint.f fiint1.f tapma.f \
	simin.f mixup.f fxis1.f fxis2.f podziel.f klopot.f \
	mixr.f coord.f chmem.f pticc.f rndm.f kontur.f rk4.f \
	qfit.f gamatt.f gcf.f tcexp.f tcabs.f tasin.f tacos.f \
	openf.f effix.f adhoc.f elmt.f newcnv.f \
	splner.f spline.f splint.f cclkup.f> $(SINGLE_FILE)

