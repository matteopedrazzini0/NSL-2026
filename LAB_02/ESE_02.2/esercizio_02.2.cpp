#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cmath>
#include"random.h"
#include"rw_functions.h"
#include<armadillo>

using namespace std;
using namespace arma;

int main(int argc, char* argv[]){
    Random rnd;
    rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");

    double a = 1.0;
    int N_blocks = 100;
    int N_Step = 100;
    int block_size = 100;

    vector<array<double, 3>> rw(N_Step);
    rw_lattice(rnd, rw, a);
    print_position(rw, "../OUTPUT/coord_lattice.xyz");

    rw_continuum(rnd, rw, a);
    print_position(rw, "../OUTPUT/coord_ontinuum.xyz");

    vec block_acc(N_Step), global_ave(N_Step), global_ave2(N_Step);

    // discrete
    for(int i=0; i<N_blocks; i++){
        block_acc.zeros();
        for(int j=0; j<block_size; j++){
            rw_lattice(rnd, rw, a);
            for(int k=0; k<N_Step; k++){
                block_acc.at(k) += compute_square_displ(rw[k]) / block_size;
            }
        }
        global_ave += arma::sqrt(block_acc);
        global_ave2 += block_acc;
    }
    print_distances(global_ave, global_ave2, N_blocks, "../OUTPUT/rw_lattice_dist.out");

    // continuum
    global_ave.zeros(); // resetting acc
    global_ave2.zeros();
    for(int i=0; i<N_blocks; i++){
        block_acc.zeros();
        for(int j=0; j<block_size; j++){
            rw_continuum(rnd, rw, a);
            for(int k=0; k<N_Step; k++){
                block_acc.at(k) += compute_square_displ(rw[k]) / block_size;
            }
        }
        global_ave += arma::sqrt(block_acc);
        global_ave2 += block_acc;
    }
    print_distances(global_ave, global_ave2, N_blocks, "../OUTPUT/rw_continuum_dist.out");

    rnd.SaveSeed();
    return 0;
}