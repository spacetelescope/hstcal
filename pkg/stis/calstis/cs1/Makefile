CC := gcc
CFLAGS = -I ../../

OBJECTS := blevdrift.o blevfit.o calstis1.o checkvolt.o do2d.o \
doatod.o dobias.o doblev.o dodark.o dodqi.o doflat.o dolores.o \
donoise.o dononlin.o dophot.o doppconv.o doshad.o findbin.o \
findblev.o findover.o flagfilter.o getccdtab.o getepctab1.o \
getflags1.o getgrpinfo1.o getkeyinfo1.o getlintab.o \
gettdccorr.o gettds1.o history1.o makedopp.o photmode.o sanity1.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
