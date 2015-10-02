CC := gcc
CFLAGS = -I ../../

OBJECTS := calstis4.o convslit.o echshift.o findbars.o findedge.o \
findends.o findwl.o flagcr.o getangle4.o getapdes4.o getdisp4.o \
getflags4.o getgrpinfo4.o getinang4.o getkeyinfo4.o getlamp.o gettrace4.o \
getwcp.o history4.o invert.o maketemplate.o newtrace4.o peakquad3.o \
scaleref.o spatialshift.o updateshift.o wavecal.o waveshift.o whichslit.o \
xc4.o xcechelle.o xcpeak.o xcwave.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
