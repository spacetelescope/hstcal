CFITSIO := /user/hodge/tables/cfitsio

CC := gcc
CFLAGS = -Wall -I ./ -I ../ -I ../cfitsio/

OBJECTS := hstio.o keyword.o numeric.o

# .INTERMEDIATE: $(OBJECTS)

../libhstio.a:  $(OBJECTS)
	ar -r -s ../libhstio.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
