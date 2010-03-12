CC := gcc

../libapplib.a:  ncarfft.o
	ar -r -s ../libapplib.a ncarfft.o

ncarfft.o: ncarfft.f
	f77 -c ncarfft.f
