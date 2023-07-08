
RDIR=..
## C compiler
CC=cc
INCLUDEDIR=-I$(RDIR)/include
LIBDIR=-L$(RDIR)/lib
LIBS=-lvsip -lm
OPTIONS=-O3

test_vsin: test_vsin.c
	$(CC) -o test_vsin test_vsin.c $(OPTIONS) $(INCLUDEDIR) $(LIBDIR) $(LIBS) 

clean:
	rm -f test_vsin test_vsin.exe
