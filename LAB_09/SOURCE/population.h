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

    void initialize(); // R i meant to be the radius of the circle or the side of the square
    double Cost_function(const vector<int>& chromo);
    //function<double(int, int)> distance;
    void sort_pop();
    Individual selector();
    void crossover(Individual& selected1, Individual& selected2); // in this scheme crossover and mutation completely substitute the previous generators
    void mutation(int mut_type, Individual& selected);
    bool check_function(const vector<int>& chr);
    void cloning();
    void evolution(int ev_step);
    double get_cost(int index);
    double get_half_average();
    void Print_final_Path();

    double distance(int city1, int city2);
    //double distance_circ(int city1, int city2);

    private:
    string _topology;
    double _r;
    double P_swap, P_shift, P_perm, P_inv, P_cross;
    Random _rnd;
    int N_ind, N_genes;
    vector<Individual> Pop;
    Individual elitario;
    vector<vector<double>> map;
    string _filename_loss, _filename_path;
};