CC=g++
CFLAGS=-g -c -O2 -std=c++11

all: main
 
main: main.o peptide.o mhstate.o wl.o scorer.o se.o metropolis.o
	g++ -g main.o peptide.o mhstate.o wl.o scorer.o se.o metropolis.o -o main 

.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *.o main
	
