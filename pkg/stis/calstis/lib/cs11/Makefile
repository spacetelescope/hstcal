CC := gcc
CFLAGS = -I ../../

OBJECTS := binsubtract.o calstis11.o getkeyinfo11.o gettimes11.o \
matchsci.o subsci.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
