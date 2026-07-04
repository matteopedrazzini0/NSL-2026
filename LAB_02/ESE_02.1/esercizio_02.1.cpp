#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include"random.h"
#include"function.h"

using namespace std;
 
int main (int argc, char *argv[]){
    Random rnd;
    rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");
    int M{100000}; // number of throws
    int N{100}; //numer of blocks
    vector<double> data(M), I(N), err_I(N);

    // 1) using normal sampling
    for(int i{}; i<M; i++){
        double x = rnd.Rannyu();
        data[i] = (M_PI/2) * cos(x * ((M_PI/2)));
    }
    compute_prog_mean(data, I, err_I);
    save("../OUTPUT/uniform_sampling.out", I, err_I);

    // 2) using importance sampling with first order Taylor appr in x = 1
    for(int i{}; i<M; i++){
        double x = retta_prob(rnd);
        data[i] = (M_PI/4) * cos(x * ((M_PI/2))) / (1 - x);
    }
    compute_prog_mean(data, I, err_I);
    save("../OUTPUT/importance_sampling.out", I, err_I);

    // show how the results change with an other probability distribution like the exponential one
    for(int i{}; i<M; i++){
        double x = rnd.exp(1);
        data[i] = (M_PI/2) * cos(x * ((M_PI/2))) / (exp(-x));
    }
    compute_prog_mean(data, I, err_I);
    save("../OUTPUT/expo_sampling.out", I, err_I);

    rnd.SaveSeed();
    return 0;
}