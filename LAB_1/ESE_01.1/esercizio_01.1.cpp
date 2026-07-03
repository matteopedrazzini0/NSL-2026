#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include"../UTILS/random.h"
#include"../UTILS/function.h"


using namespace std;
 
int main (int argc, char *argv[]){
   Random rnd;
   rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");

   //punto 1) calcolo media ed errore in funzione del numero di blocchi
   int M{100000};
   int N{100}; //numero di blocchi usati

   vector<double> data(M);
   vector<double> sum_prog(N), err_(N);
   for(int i=0; i<M; i++){
      data[i] = rnd.Rannyu();
   }
   compute_prog_mean(data, sum_prog, err_);
   save("../OUTPUT/meanvalue.out", sum_prog, err_);

   //punto 2) calcolo l'errore quadratico medio e il suo errore
   for(int i=0; i<M; i++){
      data[i] = pow( (rnd.Rannyu()-0.5), 2);
   }
   compute_prog_mean(data, sum_prog, err_);
   save("../OUTPUT/variance_value.out", sum_prog, err_);
   
   //PUNTO 3)
   int Q{100}, n{10000}; //Q = numero di sottointerballi per dividere [0, 1], n = numero lanci
   double xmin{}, xmax{1};
   vector<double> X2(N);
   //compute X2 for all blocks
   for(int i=0; i<N; i++){
      X2[i] = compute_X2(Q, n, xmin, xmax, rnd);
   }
   ofstream wr;
   wr.open("../OUTPUT/X2.out");
   double chi_mean{};
   if (wr.is_open()){
      for(int i{}; i<Q; i++)
      {
         wr << X2[i] << "\n";
         chi_mean += X2[i];
      }
      chi_mean /= Q;
   } else cerr << "Problem to open X2.out" << "\n";
   wr.close();
   
   cout << "<X^2> = " << chi_mean << "\n";

   rnd.SaveSeed();
   return 0;
}