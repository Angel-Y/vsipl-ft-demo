RDIR=..
CC=cc
INCLUDEDIR=-I$(RDIR)/include
LIBDIR=-L$(RDIR)/lib
LIBS=-lvsip -lfftw3f -lm
OPTIONS=-O3

example: test_pulse_compression_4096.c
	$(CC) -o test_pulse_compression_4096 test_pulse_compression_4096.c $(OPTIONS) $(INCLUDEDIR) $(LIBDIR) $(LIBS) 

clean:
	rm -f test_pulse_compression_4096 test_pulse_compression_4096.exe
