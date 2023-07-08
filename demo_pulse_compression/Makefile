RDIR=..
CC=cc
INCLUDEDIR=-I$(RDIR)/include
LIBDIR=-L$(RDIR)/lib
LIBS=-lvsip -lfftw3f -lm
OPTIONS=-O3

example: test_pulse_compression.c
	$(CC) -o test_pulse_compression test_pulse_compression.c $(OPTIONS) $(INCLUDEDIR) $(LIBDIR) $(LIBS) 

clean:
	rm -f test_pulse_compression test_pulse_compression.exe
