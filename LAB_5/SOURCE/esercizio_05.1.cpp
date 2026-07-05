#include<iostream>
#include<fstream>
#include<vector>
#include"random.h"
#include"function.h"

using namespace std;

int main (int argc, char *argv[]){
    Random rnd;
    rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");
    
    int L{10000}; // block size
    int N{100}; //numer of blocks
    int eq_step = 500;
    vector<double> data(L*N), r(N), err_r(N);
    
    double rate = 2.7;
    bool Gauss = false;
    array<int, 3> nlm = {2, 1, 0};
    bool is100 = (nlm[0]==1 && nlm[1]==0 && nlm[2]==0);
    
    Files filenames = get_Files(is100, Gauss, rate);
    ofstream cout_eq_acc(filenames.eq_acceptance), cout_points(filenames.points), cout_acc_run(filenames.run_acceptance);
    is_open(cout_eq_acc, filenames.eq_acceptance);
    is_open(cout_points, filenames.points);
    is_open(cout_acc_run, filenames.run_acceptance);
    cout_eq_acc << "#step" << setw(12) << "acceptance" << endl;
    cout_points << "#" << setw(14) << "x" << setw(14) << "y" << setw(14) << "z" << endl;
    cout_acc_run << "# Block" << setw(12) << "acceptance" << endl;

    int accettanza = 0;
    array<double, 3> x = {10.0, 0, 0};
    //equilibrazione
    for(int i{}; i<eq_step; i++){
        x = sample(rnd, x, rate, accettanza, nlm, Gauss);
        cout_eq_acc << i << setw(12) << double(accettanza) / (i+1) << endl;
    }
    cout_eq_acc.close();
    cout << "accettanza finale equilibrazione: " << double(accettanza)/eq_step << endl;

    //sample points
    for(int i{}; i<N; i++){
        accettanza = 0;
        for(int j=0; j<L; j++){
            x = sample(rnd, x, rate, accettanza, nlm, Gauss);
            data[i*L+j] = sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
            if(j%100 == 0) cout_points << setw(14) << x[0] << setw(14) << x[1] << setw(14) << x[2] << endl;
        }
        cout_acc_run << i << setw(12) << double(accettanza) / double(L) << endl;
    }
    cout_acc_run.close();
    cout_points.close();

    compute_prog_mean(data, r, err_r);
    save(filenames.radius, r, err_r);
    
    rnd.SaveSeed();
    return 0;
}