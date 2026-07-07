#include"VMC_SA.h"
#include<cmath>
#include<fstream>

using namespace std;

void VMC_SA :: block_reset(int blk){ // Reset block accumulators to zero and restart from random point
  ofstream coutf;
  if(blk>0){
    coutf.open("../OUTPUT/output.dat",ios::app);
    coutf << "Block completed: " << blk << endl;
    coutf.close();
  }
  _block_ave = 0;
  _n_accepted = 0;
  _n_attempts = 0;
}

double VMC_SA :: error(double sum_av, double sum_av2, int n){
    return (n==0) ? 0 : sqrt(fabs(sum_av2 / double(n) - pow(sum_av / double(n), 2)) / n);
}

void VMC_SA :: prepare(){
    int p1, p2;
    ifstream Primes("../INPUT/Primes.txt");
    if(Primes.is_open()){
        Primes >> p1 >> p2;
    } else cerr << "PROBLEM: Unable to open Primes" << endl;
    Primes.close();
    int seed[4];
    ifstream input("../INPUT/seed.in");
    string property;
    if(input.is_open()){
        while(!input.eof()){
            input >> property;
            if(property == "RANDOMSEED"){
                input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
                _rnd.SetRandom(seed, p1, p2);
            }
        }
        input.close();
    } else cerr << "PROBLEM: Unable to open seed.in" << endl;

    ofstream couta("../OUTPUT/acceptance.dat"); // Set the heading line in file ../OUTPUT/acceptance.dat
    if(!couta.is_open()){
        cerr << "PROBLEM: Unable to open acceptance.dat" << endl;
        exit(1);
    }
    couta << "# BLOCK:" << setw(14) << "ACCEPTANCE:" << endl;
    couta.close();

    ifstream inp("../INPUT/input.txt");
    if(!inp.is_open()){
        cerr << "PROBLEM: Unable to open input.txt" << endl;
        exit(1);
    }
    ofstream coutf("../OUTPUT/output.dat");
    if(!coutf.is_open()){
        cerr << "PROBLEM: Unable to open output.dat" << endl;
        exit(1);
    }
    while(inp >> property){
        if(property == "TYPE") inp >> _type;
        else if(property == "STEP_NUMBER") inp >> _n_step; // step per il calcolo dell'energia, non dipende dalla temperatura questo metropolis
        else if(property == "BLOCK_NUMBER") inp >> _n_blocks;
        else if(property == "BLOCK_NUMBER_SA") inp >> _n_blocks_SA;
        else if(property == "RATE") inp >> _rate;
        else if(property == "RATE_PARAM") inp >> _rate_param;
        else if(property == "N_STEP_T_START") inp >> _N_0;
        else if(property == "START_TEMPERATURE"){
            inp >> _T_in;
            _beta_in = 1.0 / _T_in;
            _beta = _beta_in;
        }
        else if(property == "COOLING") inp >> _cooling;
        else if(property == "END_TEMPERATURE"){
            inp >> _T_f;
            _beta_f = 1.0 / _T_f;
        }
        else if(property =="NUM_BINS"){
            inp >> _nbins;
            _x_min = -2.5, _x_max = 2.5;
            _bin_size = (_x_max - _x_min) / _nbins;
            pdf.zeros(_nbins);
        }
        else if(property == "ENDINPUT") break;
        else{
            cerr << "PROBLEM: unknown input" << endl;
            exit(1);
        }
    }
    inp.close();
    coutf << "System initialized" << endl;
    coutf.close();

    if(_type==-1){
        _filename = "../OUTPUT/Mean_energy_no_opt_p.dat";
        ofstream coutf(_filename);
        if(!coutf.is_open()){
            cerr << "Problem to open " << _filename << endl;
            exit(1);
        }
        coutf << setw(14) << "#BLOCK" << setw(14) << "ACTUAL_E" << setw(14) << "E_MEAN" << setw(14) << "ERR_E" << endl;
        coutf.close();
    }
    if(_type==0){
        // schema: T_(n+1) = T_n / c (beta_(n+1)=beta_n*c) con c>1, allora T_f = T_0 * c^(-N)
        _SA_STEP = int(log(_T_in / _T_f) / log(_cooling)) + 1;
        ofstream coutf, outf, out_acc_SA, acc_space;

        coutf.open("../OUTPUT/output.dat", ios::app);
        coutf << "Simulated annealing optimization running " << endl;

        outf.open("../OUTPUT/parameters_evolution.dat");
        outf << "#T" << setw(14) << "mu" << setw(14) << "sigma" <<setw(14) << "E" << setw(14) << "ERR_E" << endl;

        out_acc_SA.open("../OUTPUT/acceptance_param.dat");
        out_acc_SA << "#T" << setw(14) << "A" << endl;

        acc_space.open("../OUTPUT/acceptance_space_temp.dat");
        acc_space << "#T" << setw(20) << "MEAN_ACC_SPACE" << endl;

        coutf.close();
        outf.close();
        out_acc_SA.close();
        acc_space.close();
    }
    if(_type==1 || _type==2){
        _filename = "../OUTPUT/Mean_energy_opt_p.dat";
        ofstream coutf(_filename);
        if(!coutf.is_open()){
            cerr << "Problem to open " << _filename << endl;
            exit(1);
        }
        coutf << setw(14) << "#BLOCK" << setw(14) << "ACTUAL_E" << setw(14) << "E_MEAN" << setw(14) << "ERR_E" << endl;
        coutf.close();
    }
    //hbar and mass
    _hbar = 1.0;
    _m = 1.0;
    //initialize acceptance
    _n_accepted = 0;
    _n_attempts = 0;
    _n_accepted_p = 0;
    _n_attempts_p = 0;
    // initialize mu and sigma, and x centered in mu with randomness
    _mu = 1.2, _sigma = 0.5;
    _x = _mu + _sigma * _rnd.Rannyu(-1.0, 1.0);
}

double VMC_SA :: psi_T(double x){
    return exp(-pow((x - _mu), 2) / (2 * _sigma * _sigma)) +
           exp(-pow((x + _mu), 2) / (2 * _sigma * _sigma));
}

double VMC_SA :: prob(double x){
    return psi_T(x)*psi_T(x);
}

bool VMC_SA :: metro(double ratio){
    double A = min(1.0, ratio);
    return _rnd.Rannyu() <= A;
}

void VMC_SA :: move(){ // allow to sample points from the pdf defined by psi^2 via M(RT)^2
    double delta_x = _rate * _sigma; // cambio dinamicamente il passo su x in base alla larghezza della funzione d'onda da campionare
    double new_x = _x + delta_x*_rnd.Rannyu(-1.0, 1.0);
    _n_attempts++;
    if(this->metro(prob(new_x) / prob(_x))){
        _x = new_x;
        _n_accepted++;
    }
}

void VMC_SA :: move_param(array<double, 2>& E_old){ // passo E_old come parametro per averla più efficente e non ricalcolare due energie due volte
    double x_old = _x; // lo salvo perchè se rifiuto la mossa riparto con un x distribuito secondo la vecchia distribuzione
    double old_mu = _mu, old_sigma = _sigma;
    // New parameters proposal
    _mu += _rate_param * _rnd.Rannyu(-1.0,1.0); //trial random parameter
    _sigma += _rate_param * _rnd.Rannyu(-1.0,1.0);
    //sigma must be positive
    if(_sigma <= 0.1){
        _mu = old_mu, _sigma = old_sigma;
        _x = x_old;
        _n_attempts++;
        return;
    }
    array<double, 2> E_new = Energy_SA();
    double ratio = exp(-_beta * (E_new[0]-E_old[0]));
    if(!this->metro(ratio)){
        _mu = old_mu, _sigma = old_sigma;
        _x = x_old;
    }
    else{
        _n_accepted_p++;
        E_old = E_new;
    }
    _n_attempts_p++;
}

double VMC_SA :: T(){
    return -0.5*(exp(-pow((_x - _mu), 2) / (2.0 * pow(_sigma, 2))) *
        (pow((_x - _mu), 2) / pow(_sigma, 4) - 1.0 / pow(_sigma, 2)) +
        exp(-pow((_x + _mu), 2) / (2.0 * pow(_sigma, 2))) *
        (pow((_x + _mu), 2) / pow(_sigma, 4) - 1.0 / pow(_sigma, 2)))/psi_T(_x);
    /*
    double sigma2{_sigma * _sigma};
    auto psi_1_T = [&](int sign){ return (-1/sigma2) * ((_x-_mu)/sigma2 * exp(-pow(_x-_mu, 2) / sigma2) + // first derivatives of psi
                                                        sign * (_x+_mu) * exp(-pow(_x+_mu, 2) / sigma2));};
    return (-_hbar*_hbar/(2*_m)) * (-1/sigma2 - (_x/sigma2 * psi_T(_x)) * psi_1_T(1) + (_mu/(sigma2 * psi_T(_x))) * psi_1_T(-1)); // second derivatives
    */
}

double VMC_SA :: V(){
    return (pow(_x, 4) - 2.5 * pow(_x, 2));
}

void VMC_SA :: H(){ // perform one measurement of energy
    double E = 0.0;
    for(int i=0; i<_n_step; i++){
        this->move();
        E += T() + V();
    }
    _block_ave += E / double(_n_step);
}

array<double, 2> VMC_SA :: Energy_SA(){ // uso semrpe nstep per calcolare l'energia perchè il suo calcolo non dipende dalla temperatura, fissati dei parametri si calcola con la psi
    // Breve transiente di equilibrazione con i nuovi parametri
    for(int i=0; i<500; i++) this->move();

    _n_attempts = 0, _n_accepted = 0;
    double block_E = 0.0, block_E2 = 0.0;
    for(int i=0; i<_n_blocks_SA; i++){
        double E = 0.0;
        for(int j=0; j<_n_step; j++){
            this->move();
            E += T() + V();
        }
        E /= double(_n_step);
        block_E += E;
        block_E2 += E*E;
    }
    _cum_acc_temp += double(_n_accepted) / double(_n_attempts);
    return {block_E / _n_blocks_SA, this->error(block_E / _n_blocks_SA, block_E2 / _n_blocks_SA, _n_blocks_SA)};
}

void VMC_SA :: average(int blk){ // stampa energie per blocco e somma cumulata, inoltre stampa anche l'accettanza per blocco
    ofstream coutf(_filename, ios::app);

    _average = _block_ave;
    _global_av += _average;
    _global_av2 += _average * _average;

    coutf << setw(14) << blk
          << setw(14) << _average
          << setw(14) << _global_av / double(blk)
          << setw(14) << this->error(_global_av, _global_av2, blk)
          << endl;
    coutf.close();

    // ACCEPTANCE
    coutf.open("../OUTPUT/acceptance.dat", ios::app);
    double acceptation = (_n_attempts > 0) ? double(_n_accepted) / double(_n_attempts) : 0;
    coutf << blk << setw(14) << acceptation << endl;
    coutf.close();
}

void VMC_SA :: optimize(){
    ofstream outf, out_acc_SA, coutf, cout_acc_space;
    cout_acc_space.open("../OUTPUT/acceptance_space_temp.dat", ios::app);
    coutf.open("../OUTPUT/output.dat", ios::app);
    outf.open("../OUTPUT/parameters_evolution.dat", ios::app);
    out_acc_SA.open("../OUTPUT/acceptance_param.dat", ios::app);

    array<double,2> E = Energy_SA();
    for(int i=0; i<_SA_STEP; i++){
        _n_attempts_p = 0, _n_accepted_p = 0;
        _cum_acc_temp = 0;
        int N_moves = int(_N_0 * pow(_beta / _beta_in, 0.3));
        for(int j=0; j<N_moves; j++){ // A bassa teperatura devo chiamare più volte move param perchè l'accettazione di move param dipende dalla temperatura, non quella del calcolo dell'energia
            move_param(E);
        }
        cout_acc_space << 1.0 / _beta << setw(14) << _cum_acc_temp / double(N_moves) << endl;
        outf << 1.0 / _beta << setw(14) << _mu << setw(14) << _sigma << setw(14) << E[0] << setw(14) << E[1] << endl;
        out_acc_SA << 1.0 / _beta << setw(14) << double(_n_accepted_p) / double(_n_attempts_p) << endl;
        _beta *= _cooling;
    }

    coutf << "Optimization completed!\n" << "Optimal parameters: mu = " << _mu << " sigma = " << _sigma << endl;
    coutf.close();
    outf.close();
    cout_acc_space.close();
    out_acc_SA.close();
}

void VMC_SA :: sample_pdf(int n_sample){ // after optimize parameters, we move the system and sample the pdf counting how many times we have a certain position
    ofstream outf("../OUTPUT/pdf.dat");
    for(int i=0; i<n_sample; i++){
        this->move();
        int bin = int( (_x-_x_min)/_bin_size );
        if(bin>=0 && bin<_nbins) pdf(bin)++;
    }
    // normalisation
    double norm = arma::sum(pdf) * _bin_size;
    pdf /= norm;
    for(int i=0; i<_nbins; i++){
        outf << pdf(i) << endl;
    }
    outf.close();
}

void VMC_SA :: Load_param(string filename){
    ifstream load(filename);
    if(!load.is_open()){
        cerr << "Impossible to open file " << filename << endl;
        exit(1);
    }
    string dummy;
    getline(load, dummy);
    double T, mu, sigma, e, err_e;
    while(load >> T >> mu >> sigma >> e >> err_e){
        _mu = mu;
        _sigma = sigma;
    } 
    cout << "E_min = " << e << "  mu = " << _mu << "  sigma = " << _sigma << endl;
    cout << "Ho finito di leggere" << endl;
    load.close();
}

void VMC_SA :: finalize(){
  _rnd.SaveSeed();
  ofstream coutf;
  coutf.open("../OUTPUT/output.dat",ios::app);
  coutf << "Simulation completed!" << endl;
  coutf.close();
}

void VMC_SA :: reset_global(){
    _average = 0.0;
    _global_av = 0.0;
    _global_av2 = 0.0;
    _block_ave = 0.0;
}

int VMC_SA :: get_Nbl(){ return _n_blocks; }

int VMC_SA :: get_Nsteps(){ return _n_step; }

void VMC_SA :: printmusigma(){
    cout << _mu << " " << _sigma << endl;
}