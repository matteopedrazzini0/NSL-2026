#include<iostream>
#include<fstream>
#include<cmath>
#include"rw_functions.h"

using namespace std;
using namespace arma;

double error(double acc, double acc2, int n){
    return (n==0) ? 0 : sqrt(fabs(acc2/double(n) - pow( acc/double(n) ,2)) / n);
}

array<double, 3> operator+(const array<double, 3>& a, const array<double, 3>& b) {
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

void rw_lattice(Random& r, vector<array<double, 3>>& rw, double a){
    rw[0] = {0, 0, 0};
    int time = int(rw.size());
    for(int i{1}; i<time; i++){
        int sign = (r.Rannyu() <= 0.5) ? -1 : 1;
        double ax = r.Rannyu(0, 3); //selziono una direzione a random
        array<double, 3> incr;

        if(ax<=1) incr = {a*sign, 0, 0};
        else if (1<ax && ax<=2) incr = {0, a*sign, 0};
        else incr = {0, 0, a*sign};

        rw[i] = rw[i-1] + incr;
    }
}

void rw_continuum(Random& r, vector<array<double, 3>>& rw, double a){
    rw[0] = {0, 0, 0};
    int time = int(rw.size());
    for(int i{1}; i<time; i++){
        double theta = r.Rannyu(0, M_PI);
        double phi = r.Rannyu(0, 2*M_PI);
        array<double, 3> incr{a * cos(phi) * sin(theta), a * sin(phi) * sin(theta), a * cos(theta)};
        rw[i] = rw[i-1] + incr;
    }
}

double compute_square_displ(array<double, 3>& r_i){
    return r_i[0]*r_i[0]+ r_i[1]*r_i[1] + r_i[2]*r_i[2];
}

void print_position(vector<array<double, 3>>& rw, string filename){
    int time = int(rw.size());
    ofstream coutf;
    coutf.open(filename);
    if(!coutf.is_open()){
        cerr << "Problem: Unable to open file " << filename << endl;
        exit(1);
    }
    coutf << "# x"<< setw(14) << "y" << setw(14) << "z" << endl;
    for(int i{}; i<time; i++){
        coutf << "  " << rw[i][0] << setw(14) << rw[i][1] << setw(14) << rw[i][2] << endl;
    }
    coutf.close();
}

void print_distances(vec acc, vec acc2, int blk, string filename){
    int N_Step = int(acc.n_elem);
    ofstream coutf;
    coutf.open(filename);
    if(!coutf.is_open()){
        cerr << "Problem: Unable to open file " << filename << endl;
        exit(1);
    }
    coutf << "#step" << setw(14) << "√r^2" << setw(14) << "err" << endl;
    for(int k=0; k<N_Step; k++){
        coutf << k << setw(14) << acc.at(k) / blk << setw(14) << error(acc.at(k), acc2.at(k), blk) << endl;
    }
}