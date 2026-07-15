CC = mpicxx
CFLAGS = -O3 --std=c++17 -Wall

all: GA_par.exe

GA_par.exe : esercizio_10.1.o population.o random.o
	$(CC) esercizio_10.1.o population.o random.o -o GA_par.exe $(CFLAGS)
esercizio_10.1.o : esercizio_10.1.cpp population.h
	$(CC) -c esercizio_10.1.cpp -o esercizio_10.1.o $(CFLAGS)
population.o : population.cpp population.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean : 
	rm -f *.o GA_par.exe