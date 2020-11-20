F = -g 
H = init.h solve.h
all: a.out

a.out: main.o solve.o init.o Makefile
	g++ main.o solve.o init.o -o a.out
main.o: main.cpp $H
	g++ $F -c main.cpp
solve.o: solve.cpp solve.h
	g++ -O3 -Ofast -march=native $F -c solve.cpp
init.o: init.cpp init.h
	g++ $F -c init.cpp

clear:
	rm *.o a.out

