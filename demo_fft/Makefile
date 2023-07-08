RDIR=..
CC=cc
INCLUDEDIR=-I$(RDIR)/include
LIBDIR=-L$(RDIR)/lib
LIBS=-lvsip -lfftw3f -lm
OPTIONS=-O3

example: test_vsipl_fft.c
	$(CC) -o test_vsipl_fft test_vsipl_fft.c $(OPTIONS) $(INCLUDEDIR) $(LIBDIR) $(LIBS) 

clean:
	rm -f test_vsipl_fft test_vsipl_fft.exe
