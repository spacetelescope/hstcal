CC := gcc
CFLAGS = -I ../../

OBJECTS := addshifts.o calstis12.o getkeyinfo12.o getscigrp.o \
getwavgrp.o history12.o matchwav.o targpos.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
