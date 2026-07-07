#pragma once

#include<iostream>
#include<array>
#include<string>
#include"random.h"
#include<algorithm>
#include<armadillo>

using namespace std;
using namespace arma;

class VMC_SA{

    public:

    void prepare();
    int get_type() { return _type; }
    int get_Nbl(); // return the number of blocks used
    int get_Nsteps(); // return the number of steps per block
    void block_reset(int blk); // reset the measurment container
    void reset_global();
    void average(int blk);
    double psi_T(double x);
    double prob(double x);
    double error(double av, double av2, int block);
    //void save(const string filename, const vector<double>& sum_prog, vector<double>& err_prog);
    double T();
    double V();
    void H();
    array<double, 2> Energy_SA();
    bool metro(double ratio);
    void move(); // Used to move in the coordinate space via a specific probability distribution
    void move_param(array<double, 2>& E_old); // Used to optimize parameter moving in the parameters space with the MRT alg
    void optimize();
    void sample_pdf(int n_sample);
    void Load_param(string filename); // if a previous simulation have found the proper optimized parameters
    void finalize();
    void printmusigma();

    private:

    int _type;
    double _hbar, _m;
    Random _rnd;
    double _sigma, _mu;
    int _n_step, _n_blocks, _n_blocks_SA; // step at each block of calcuating energy, number of blocks
    double _block_ave; // contain the single measurement inside the single block
    double _average, _global_av, _global_av2; // contain the progressive average made on the blocks
    int _n_accepted, _n_attempts; // For the acceptance in the computation of energy
    int _n_accepted_p, _n_attempts_p; // for the acceptance in the optimization of parameters
    double _cum_acc_temp;
    double _rate, _rate_param; // len of the metropolis step for computing energy / optimizing parameters respectively
    double _x, _x_min, _x_max; //position of the system
    double _T_in, _beta_in; // Initial temperature and inverse temperature to start SA
    double _beta; // Actual inverse temperature of the system
    double _T_f, _beta_f; // Final temperature and inverse temperature to stop SA
    double _cooling; // Coefficient to cool the system during optimization
    int _N_0; // Number of step to performe parameters espoloration at T_start
    int _SA_STEP; // How many temperatures we pass given start T, final T and cooling coefficient
    int _nbins;
    double _bin_size;
    string _filename;
    vec pdf;
};