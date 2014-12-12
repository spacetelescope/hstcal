CC := gcc
CFLAGS = -I ../../../

OBJECTS := calstis6idt.o echscatread.o idtmemory.o rebindata.o \
splice.o x1dtable.o xc6.o

# .INTERMEDIATE: $(OBJECTS)

../../calstis.a:  $(OBJECTS)
	ar -r -s ../../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
