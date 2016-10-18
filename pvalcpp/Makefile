CC=g++
CFLAGS=-g -c -O2 -std=c++11


all: main
	./main

main: main.o peptide.o wl.o psm.o main_functions.o se.o
	g++ -g main.o peptide.o wl.o psm.o main_functions.o se.o -o main 


.SUFFIXES: .cpp .o

.cpp.o:
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm -rf *.o main_tmp
	
