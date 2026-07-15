CC = mpicxx
INCLUDES = -I../../PARALLEL_GA/SOURCE/
CFLAGS = -O3 --std=c++17 -Wall $(INCLUDES)

all: GA_par.exe

GA_par.exe : esercizio_10.1.o population.o random.o
	$(CC) esercizio_10.1.o population.o random.o -o GA_par.exe $(CFLAGS)
esercizio_10.1.o : esercizio_10.1.cpp ../../PARALLEL_GA/SOURCE/population.h ../../PARALLEL_GA/SOURCE/random.h
	$(CC) -c esercizio_10.1.cpp -o esercizio_10.1.o $(CFLAGS)
population.o : ../../PARALLEL_GA/SOURCE/population.cpp ../../PARALLEL_GA/SOURCE/population.h ../../PARALLEL_GA/SOURCE/random.h
	$(CC) -c ../../PARALLEL_GA/SOURCE/population.cpp -o population.o $(CFLAGS)
random.o : ../../PARALLEL_GA/SOURCE/random.cpp ../../PARALLEL_GA/SOURCE/random.h
	$(CC) -c ../../PARALLEL_GA/SOURCE/random.cpp -o random.o $(CFLAGS)
clean : 
	rm -f *.o GA_par.exe