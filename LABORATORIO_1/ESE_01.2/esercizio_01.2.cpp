/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include<iostream>
#include<fstream>
#include<string>
#include"../UTILS/random.h"
#include<cmath>
#include<vector>
#include<array>

using namespace std;

int main(int argc, char* argv[]){
    Random rnd;
    rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");

    int n_samples{10000};
    vector<int> N_Blocks{1, 2, 10, 100};
    double S_N_uni{}, S_N_exp{}, S_N_Lorentz{};
    string S_uni = "../OUTPUT/S_uniform.out", S_exp = "../OUTPUT/S_exp.out", S_Lor = "../OUTPUT/S_Lorentz.out";

    ofstream out_uni, out_exp, out_Lor;
    out_uni.open(S_uni);
    out_exp.open(S_exp);
    out_Lor.open(S_Lor);

    for(int i=0; i<(int) N_Blocks.size(); i++){
        for(int j=0; j<n_samples; j++){
            S_N_uni = 0;
            S_N_exp = 0;
            S_N_Lorentz = 0;
            for(int k=0; k<N_Blocks[i]; k++){
                S_N_uni += rnd.Rannyu();
                S_N_exp += rnd.exp(1);
                S_N_Lorentz += rnd.Lorentz(0, 1);
            }
            out_uni << N_Blocks[i] << " " << S_N_uni / N_Blocks[i] << endl;
            out_exp << N_Blocks[i] << " " << S_N_exp / N_Blocks[i] << endl;
            out_Lor << N_Blocks[i] << " " << S_N_Lorentz / N_Blocks[i] << endl;
        }
    }
    out_uni.close();
    out_exp.close();
    out_Lor.close();

    rnd.SaveSeed();
    return 0;
}