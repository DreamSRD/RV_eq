HEADER = -I/usr/local/include
LIBB = -L/usr/local/lib
LIBRA = -lm -lfftw3
SOURCES = dzeq.c
CFLAGS = -O3

all:
	gcc -std=c99  $(CFLAGS) $(SOURCES) $(HEADER) $(LIBB) $(LIBRA) -o dzeq
clean:
	rm -rf *.o dzeq
