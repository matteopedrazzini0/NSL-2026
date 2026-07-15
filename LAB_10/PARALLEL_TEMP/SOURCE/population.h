#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<functional>
#include"random.h"

using namespace std;

class Individual{
    public:
    Individual(int N_genes) : genes(N_genes), cost(0.0) {}
    vector<int> genes;
    double cost;
    double beta;
};

class population{
    public:
    population(int N_indi, int N_gen) : N_ind(N_indi), N_genes(N_gen),
    Pop(N_indi, Individual(N_gen)), elitario(N_gen),
    map(N_gen, vector<double>(2,0.0))
    {
        if(N_indi<=0){
            cerr << "Number of individuals must be > 0" << endl;
            exit(1);
        }
        if(N_gen<=1){
            cerr << "Number of genes must be > 1" << endl;
            exit(1);
        }
    }
    ~population(){};

    void initialize(int rank, int size); // R i meant to be the radius of the circle or the side of the square
    void read_map(string filename);
    double Cost_function(const vector<int>& chromo);
    //function<double(int, int)> distance;
    void sort_pop();
    Individual& selector();
    void mutation(int mut_type, Individual& selected);
    bool check_function(const vector<int>& chr);
    void cloning();
    void evolution(int ev_step, ofstream& coutf);
    double get_cost(int index);
    double get_half_average();
    void Print_final_Path(string filename, int rank);
    double distance(int city1, int city2);
    Individual get_path(int idx);
    void substitute(vector<int>& path, int idx);
    void take_path(vector<int>& path, string filename);
    double rannyu(double a, double b);
    bool metro(double ratio);
    double get_beta();
    int get_step();
    int get_n_migr();
    void reset_acc();
    double get_acc();
    string get_filename_loss();

    private:
    double _beta, _beta_min, _beta_max; // Inverse temperature of the system
    int _nattempts, _naccepted;
    int _N_0, _N_moves, _N_mut, _N_migr;
    int _rank;
    int _size;
    string _type;
    double _r;
    Random _rnd;
    int N_ind, N_genes;
    vector<Individual> Pop;
    Individual elitario;
    vector<vector<double>> map;
    string _filename_loss, _filename_path;
};

void partnership(int rank, int size, int& next_part, int tag);
void exchange(int rank, int partner, population& pop, int Ncities, ofstream& coutll, int genblock);
void select_unique(int rank, int size, int have_min, vector<int>& who_have_mini, int& select);