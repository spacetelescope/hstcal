# Use 'make. to compile.
# Use 'make cs0.e' to link.

CFITSIO := /user/hodge/tables/cfitsio

CC := gcc
CFLAGS = -I ../../
LIBS   = ../../libhstio.a ../../libcvos.a ../../libapplib.a \
../../libxtables.a $(CFITSIO)/libcfitsio.a

OBJECTS := calstis0.o checkwav.o compfiles.o getflags.o getnewref.o \
getsciinfo.o getwavinfo.o refexist.o savename.o sciflags.o \
stisinit.o wavflags.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

cs0.e:  cs0.o ../calstis.a
	$(CC) -o cs0.e cs0.o ../calstis.a $(LIBS) -lm

cs0.o: cs0.c
	$(CC) -c $(CFLAGS) cs0.c
