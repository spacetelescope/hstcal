CC := gcc
CFLAGS = -I ../../

OBJECTS := add2d.o addk2d.o bin2d.o bincoords.o binupdate.o \
blazecorr.o checkimsetok.o defswitch.o div2d.o dostat.o extrap1d.o \
fileexists.o getcorner.o getdettemp.o getlt0.o getlt.o getmoc.o \
getrefname.o getswitch.o gettds.o gotfilename.o imghistory.o \
imgpedigree.o initreftab.o interp1d.o interp2d.o key.o median.o \
mkname.o mkoutname.o mult2d.o multk2d.o omitstep.o orbitaldopp.o \
prinfo.o pseudoap.o radialvel.o reffiles.o rotatetrace.o rowpedigree.o \
sameint.o spline.o strcmptail.o streqic.o strtor.o sub2d.o \
tabhistory.o tabpedigree.o tdscorrect.o timestamp.o ucalver.o \
ufilename.o unbin2d.o whicherror.o xcfft.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
