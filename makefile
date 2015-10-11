CC=g++
# CFLAGS=-Wall -I. -fopenmp -O3 -pedantic-errors -Wno-write-strings
CFLAGS=-Wall -I. -fopenmp -O3  -Wno-write-strings
LIBS=-lm

# SRCS=gridtmz.c materialtmz.c pmltmz.c printData.c sources.c updatetmz.c inputtmz.c main.c setup.c
# OBJS= $(SRCS:.c=.o)
#ALLDEFS  -DPOINT_SOURCE -DOPENMP -DDRUDE-DPLANE_SOURCE -DCENTER_RICKER_SOURCE -DLEFT_RICKER_SOURCE
# STDDEFS= -DPOINT_SOURCE -DOPENMP
# DRU = -DDRUDE
hybrid:
	mpicxx -o hybrid.exe $(CFLAGS) hybridCode.cpp

mpi:
	mpic++ -o mpi.exe $(CFLAGS) mpiCode.cpp

omp:
	$(CC) -o omp.exe $(CFLAGS) ompCode.cpp

serial:	
	$(CC) -o serial.exe $(CFLAGS) serialCode.cpp
clean:
	rm *.o	
	rm *.exe

tidy:
	rm *~
