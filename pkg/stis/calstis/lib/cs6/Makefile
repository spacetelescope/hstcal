CC := gcc
CFLAGS = -I ../../

OBJECTS := absflux6.o addoffsets6.o adjustdisp6.o background.o \
buildoptprof.o buildprof.o calstis6.o calstis6std.o checkoptimal.o \
commline.o crosscor.o defineback.o do1dx.o findlya.o \
gcrosscor.o getabsphot6.o getapdes6.o getapoffset6.o getapthr6.o \
getccdtab6.o getcheckref.o getdisp6.o getextract.o getflags6.o \
getgac6.o getgrpinfo6.o getinang6.o getintens.o getkeyinfo6.o \
getpct6.o getprofile.o getrefcline.o getsdc6.o gettrace6.o \
heliofactor6.o history6.o init.o lee.o memory.o message6.o newcoord6.o \
newdisp6.o newextract.o newprofile.o newtrace6.o output.o poly.o \
selectalg.o setranges.o timestamp6.o wave.o x1dspec.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
