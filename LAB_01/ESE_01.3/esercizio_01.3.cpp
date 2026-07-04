#include<iostream>
#include<fstream>
#include<string>
#include<cmath>
#include<vector>
#include<array>
#include"../UTILS/random.h"
#include"../UTILS/function.h"

using namespace std;

double buffon_exp(double d, double l, int needle_throws, Random& r);

int main(int argc, char* argv[]){
    Random rnd;
    rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");
    double d{10};
    double l{8};
    int N_blocks{100};
    int n_thr_per_blocks{10000};
    vector<double> PI(N_blocks);
    vector<double> ave(N_blocks);
    vector<double> PI2(N_blocks);
    vector<double> ave2(N_blocks);
    vector<double> err(N_blocks);

    //loop over all bklocks
    double glob_acc = 0.0, glob_acc2 = 0.0;
    for(int i=0; i<N_blocks; i++){
        // loop inside single block
        PI[i] = buffon_exp(d, l, n_thr_per_blocks, rnd);
        PI2[i] = PI[i] * PI[i];
        // cumulative sum
        glob_acc += PI[i];
        glob_acc2 += PI2[i];
        ave[i] = glob_acc / (i+1);
        ave2[i] = glob_acc2 / (i+1);
        err[i] = error(ave[i], ave2[i], i);
    }
    save("../OUTPUT/buffon.out", ave, err);

    rnd.SaveSeed();
    return 0;
}

double buffon_exp(double d, double l, int needle_throws, Random& r){
    int succ = 0;
    for(int i=0; i<needle_throws; i++){
        double x = r.Rannyu(0, d/2);
        // generate sin(theta) without using π
        double x_c, y_c, rad;
        do{
            x_c = r.Rannyu(-1, 1);
            y_c = r.Rannyu(-1, 1);
            rad = sqrt(x_c*x_c + y_c*y_c);
        }while(rad==0 || rad>1);
        double sin_theta = y_c / rad;
        if(x<=(l/2)*fabs(sin_theta)) succ++;
    }
    return (succ>0) ? (2*l*needle_throws) / (d*succ) : 0.0;
}