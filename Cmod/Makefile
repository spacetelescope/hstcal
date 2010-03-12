CFITSIO := /user/hodge/tables/cfitsio

CC := gcc
CFLAGS = -Wall -I ./ -I $(CFITSIO)
LIBS   = $(CFITSIO)/libcfitsio.a

OBJECTS := vfn2osfn.o c_imt.o \
addCol.o cerror.o ctables.o c_tbagt.o c_tbapt.o \
c_tbcdef1.o c_tbcfnd1.o c_tbciga.o c_tbcigi.o c_tbcigt.o c_tbcinf.o \
c_tbcnum.o c_tbegt.o c_tbept.o c_tbfpri.o c_tbfres.o \
c_tbhad.o c_tbhcal.o c_tbhgcm.o c_tbhgnp.o c_tbhgt.o c_tbhpcm.o c_tbhptt.o \
c_tbparse.o c_tbpsta.o c_tbrcsc.o c_tbrgtr.o c_tbrudf.o \
c_tbtacc.o c_tbtclo.o c_tbtcre.o c_tbtnam.o c_tbtopn.o \
initCol.o tbCopyTmpl.o tbHeader.o tbl_util.o tbSaveInfo.o c_phopar.o

# .INTERMEDIATE: Test.o $(OBJECTS)

test:  Test.o $(OBJECTS)
	$(CC) -o test Test.o $(OBJECTS) $(LIBS) -lm

Test.o: Test.c
	$(CC) -c $(CFLAGS) Test.c

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

clean:
	-rm test
