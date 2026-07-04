CC = g++
CFLAGS = -Wall -O3 --std=c++17 -I$(ARMA_HOME)/include
ARMA_HOME = /opt/homebrew/opt/armadillo
AFLAGS = -L$(ARMA_HOME)/lib -larmadillo

esercizio_02.2.exe : main.o random.o rw_functions.o #dice che per creare l'eseguibile servono main.o e random .o, e sotto dice come creare main.exe a partitre dai quei due
	$(CC) random.o main.o rw_functions.o -o esercizio_02.2.exe $(CFLAGS) $(AFLAGS)
main.o : esercizio_02.2.cpp random.h rw_functions.h #qui dice cosa serve per creare main.o e sotto come costruitlo. il tag -c serve a dire di non linkare
	$(CC) -c esercizio_02.2.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
rw_functions.o : rw_functions.cpp rw_functions.h random.h
	$(CC) -c rw_functions.cpp -o rw_functions.o $(CFLAGS)
clean :
	rm *.o esercizio_02.2.exe seed.out