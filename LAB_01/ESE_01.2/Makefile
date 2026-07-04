CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio_01.2.exe : main.o random.o #dice che per creare l'eseguibile servono main.o e random .o, e sotto dice come creare main.exe a partitre dai quei due
	$(CC) random.o main.o -o esercizio_01.2.exe $(CFLAGS)
main.o : esercizio_01.2.cpp #qui dice cosa serve per creare main.o e sotto come costruitlo. il tag -c serve a dire di non linkare
	$(CC) -c esercizio_01.2.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o esercizio_01.2.exe seed.out
