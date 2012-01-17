CC := gcc
CFLAGS = -I ../../

OBJECTS := absflux.o adjustdisp.o calstis7.o datamasked.o do2dx.o \
geocorr7.o getabsphot.o getapdes7.o getapoffset.o getapthr.o \
getdisp.o getflags7.o getgrpinfo7.o getidc.o getinang.o getkeyinfo7.o \
getpct.o getsdc.o gettrace.o heliofactor.o history7.o newcoord.o \
newdisp.o newtrace.o putgrpinfo7.o scalewcs.o x2dcorr7.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
