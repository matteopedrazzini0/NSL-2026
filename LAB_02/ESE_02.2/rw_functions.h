#pragma once

#include<vector>
#include<array>
#include<string>
#include"random.h"
#include<armadillo>

using namespace std;
using namespace arma;

/*
class R_Walk{
    public:
    R_Walk(int N_steps, double a) : N_steps{N_steps}, step_length{a}, rw(N_steps, {0.0, 0.0, 0.0}){
        _rnd.Initialize_RNG("../INPUT/Primes.txt", "../INPUT/seed.in");
    }

    private:
    Random _rnd;
    int N_steps;
    double step_length;
    vector<array<double, 3>> rw;
};
*/
void rw_lattice(Random& r, vector<array<double, 3>>& rw, double a);

void rw_continuum(Random& r, vector<array<double, 3>>& rw, double a);

double compute_square_displ(array<double, 3>& r_i);

void print_position(vector<array<double, 3>>& rw, string filename);

void print_distances(vec acc, vec acc2, int blk, string filename);

double error(double acc, double acc2, int n);