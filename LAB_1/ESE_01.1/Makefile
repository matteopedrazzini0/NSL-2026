CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio_01.1.exe : main.o random.o function.o#dice che per creare l'eseguibile servono main.o e random .o, e sotto dice come creare main.exe a partitre dai quei due
	$(CC) main.o random.o function.o -o esercizio_01.1.exe $(CFLAGS)
main.o : esercizio_01.1.cpp random.h function.h #qui dice cosa serve per creare main.o e sotto come costruitlo. il tag -c serve a dire di non linkare
	$(CC) -c esercizio_01.1.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
function.o : function.cpp function.h random.h
	$(CC) -c function.cpp -o function.o $(CFLAGS)
clean :
	rm *.o esercizio_01.1.exe seed.out
