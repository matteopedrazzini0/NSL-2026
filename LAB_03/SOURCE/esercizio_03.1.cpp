#include<iostream>
#include<fstream>
#include"random.h"
#include"function.h"

using namespace std;

int main(int argc, char* argv[]){
    Random rnd;
    rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");

    double S_0 = 100.0;
    double T = 1.0;
    double k = 100.0;
    double r = 0.1;
    double sigma = 0.25;
    double disc_fact = exp(-r * T);

    int N = 100000;
    int N_blocks = 100;
    int N_step = 100;
    double dt = T / N_step;
    vector<double> C_i(N), P_i(N);
    vector<double> C(N_blocks), err_C(N_blocks), P(N_blocks), err_P(N_blocks);

    // 1) GMB
    for(int i=0; i<N; i++){
        double z_i = rnd.Gauss(0, 1);
        double S = S_0 * exp((r-0.5*sigma*sigma)*T + sigma*z_i*sqrt(T));
        C_i[i] = disc_fact * call_price(S, k); // sistema st disc factr che forse si può fare una volta in meno
        P_i[i] = disc_fact * put_price(S, k);
    }
    compute_prog_mean(C_i, C, err_C);
    compute_prog_mean(P_i, P, err_P);
    save("../OUTPUT/GBM_call.dat", C, err_C);
    save("../OUTPUT/GBM_put.dat", P, err_P);

    //2) discrete case
    ofstream asset_ev;
    asset_ev.open("../OUTPUT/asset_price_evolution.dat");
    asset_ev << "# TIME: " << setw(10) << "PRICE:" << endl;
    asset_ev << 0 << setw(10) << S_0 << endl;
    int wprint = int(rnd.Rannyu(0, N-1));

    for(int i=0; i<N; i++){
        double S = S_0;
        for(int j=1; j<=N_step; j++){
            double z_i = rnd.Gauss(0, 1);
            S = S * exp((r-0.5*sigma*sigma)*dt + sigma*z_i*sqrt(dt));
            if(i==wprint) asset_ev << j << setw(10) << S << endl;
        }
        C_i[i] = disc_fact * call_price(S, k); // sistema st disc factr che forse si può fare una volta in meno
        P_i[i] = disc_fact * put_price(S, k);
    }
    asset_ev.close();
    compute_prog_mean(C_i, C, err_C);
    compute_prog_mean(P_i, P, err_P);
    save("../OUTPUT/discrete_call.dat", C, err_C);
    save("../OUTPUT/discrete_put.dat", P, err_P);

    rnd.SaveSeed();
    return 0;
}