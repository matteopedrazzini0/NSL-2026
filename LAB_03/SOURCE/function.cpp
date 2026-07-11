#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<cmath>
#include"random.h"
#include"function.h"

using namespace std;

double error(double av, double av2, int n){  // mean dev std
    return (n==0) ? 0 : sqrt((av2 - pow(av, 2)) / n);
}

void compute_prog_mean(const vector<double>& data, vector<double>& sum_prog, vector<double>& err_prog){
    int n_Blocks = int(sum_prog.size());
    int n_measure = int(data.size());
    int L = n_measure / n_Blocks; // Block size
    vector<double> ave(n_Blocks, 0.0), ave2(n_Blocks, 0.0), sum_prog2(n_Blocks, 0.0);
    double cum_sum = 0.0, cum_sum2 = 0.0;
    
    for(int i=0; i<n_Blocks; i++){
        double sum = 0;
        for(int j=0; j<L; j++){
            sum += data[i*L+j];
        }
        ave[i] = sum/L;
        ave2[i] = ave[i]*ave[i];

        cum_sum += ave[i];
        cum_sum2 += ave2[i];

        sum_prog[i] =cum_sum / (i+1);
        sum_prog2[i] = cum_sum2 / (i+1);
        err_prog[i] = error(sum_prog[i], sum_prog2[i], i);
    }
}

void save(const string filename, const vector<double>& sum_prog, vector<double>& err_prog)
{
    ofstream wr;
    wr.open(filename);
    if(wr.is_open()){
        wr << "#blocchi" << " " << "mean_prog" << " " << "err_prog" << "\n";
        for(size_t i=0; i<sum_prog.size(); i++){
            wr << i+1 << " " << sum_prog[i] << " " << err_prog[i] << "\n";
        }
    }else cerr << "Problem to open " << filename << "\n";
   wr.close();
}

double call_price(double S, double K){
    return (S - K > 0) ? S - K : 0;
}

double put_price(double S, double K){
    return (K - S > 0) ? K - S : 0;
}