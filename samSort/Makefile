CC=gcc
CFLAGS=-c -Wall -O2
CFLAGS2=-c -Wall -O2 -lgutil -lm

all: samSort

samSort: redblackSN.o samSort.o
	$(CC) redblackSN.o samSort.o -o samSort

redblackSN.o: redblackSN.c
	$(CC) $(CFLAGS) redblackSN.c

samSort.o: samSort.c
	$(CC) $(CFLAGS2) samSort.c

clean:
	rm -rf *o
