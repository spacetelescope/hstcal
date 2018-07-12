CC := gcc
CFLAGS = -I ../../

OBJECTS := calstis2.o cr_history.o cr_mode.o crrej_check.o crrej_do.o \
crrej_init.o crrej_loop.o crrej_sky.o cr_scaling.o cs2_command.o \
o_cal2_in.o piksrt.o

# .INTERMEDIATE: $(OBJECTS)

../calstis.a:  $(OBJECTS)
	ar -r -s ../calstis.a $(OBJECTS)

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@
