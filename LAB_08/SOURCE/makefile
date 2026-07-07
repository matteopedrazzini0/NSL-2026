CC = g++
ARMA_HOME = /opt/homebrew/opt/armadillo
AFLAGS =-L$(ARMA_HOME)/lib -larmadillo
CFLAGS = -Wall -O3 --std=c++17 -I$(ARMA_HOME)/include

esercizio_08.1.exe : main.o random.o VMC_SA.o #dice che per creare l'eseguibile servono main.o e random .o, e sotto dice come creare main.exe a partitre dai quei due
	$(CC) random.o main.o VMC_SA.o -o esercizio_08.1.exe $(AFLAGS)
main.o : esercizio_08.1.cpp VMC_SA.h #qui dice cosa serve per creare main.o e sotto come costruitlo. il tag -c serve a dire di non linkare
	$(CC) -c esercizio_08.1.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
VMC_SA.o : VMC_SA.cpp VMC_SA.h
	$(CC) -c VMC_SA.cpp -o VMC_SA.o $(CFLAGS)
clean :
	rm *.o esercizio_08.1.exe seed.out
