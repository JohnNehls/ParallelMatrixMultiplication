CC=g++
CFLAGS=-Wall -I. -fopenmp -O3  -Wno-write-strings #-pedantic-errors
LIBS=-lm

mpi:
	mpic++ -o mpi.exe $(CFLAGS) mpiCode.cpp

omp:
	$(CC) -o omp.exe $(CFLAGS) ompCode.cpp

serial:	
	$(CC) -o serial.exe $(CFLAGS) serialCode.cpp

clean:
	rm *.exe
