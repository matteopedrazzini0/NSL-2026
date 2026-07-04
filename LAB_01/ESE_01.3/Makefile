CC = g++
INCLUDES = -I../UTILS
CFLAGS = -Wall -O3 --std=c++11 $(INCLUDES)

esercizio_01.3.exe : main.o random.o function.o#dice che per creare l'eseguibile servono main.o e random .o, e sotto dice come creare main.exe a partitre dai quei due
	$(CC) main.o random.o function.o -o esercizio_01.3.exe $(CFLAGS)
main.o : esercizio_01.3.cpp ../UTILS/random.h ../UTILS/function.h #qui dice cosa serve per creare main.o e sotto come costruitlo. il tag -c serve a dire di non linkare
	$(CC) -c esercizio_01.3.cpp -o main.o $(CFLAGS)
random.o : ../UTILS/random.cpp ../UTILS/random.h
	$(CC) -c ../UTILS/random.cpp -o random.o $(CFLAGS)
function.o : ../UTILS/function.cpp ../UTILS/function.h ../UTILS/random.h
	$(CC) -c ../UTILS/function.cpp -o function.o $(CFLAGS)
clean :
	rm *.o esercizio_01.3.exe seed.out
