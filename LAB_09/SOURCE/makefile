CC = g++
CFLAGS = -O3 --std=c++17 -Wall

all: GA.exe

GA.exe : esercizio_09.1.o population.o random.o
	$(CC) esercizio_09.1.o population.o random.o -o GA.exe $(CFLAGS)
esercizio_09.1.o : esercizio_09.1.cpp population.h
	$(CC) -c esercizio_09.1.cpp -o esercizio_09.1.o $(CFLAGS)
population.o : population.cpp population.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean : 
	rm -f *.o GA.exe