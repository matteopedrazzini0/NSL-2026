#pragma once

#include<vector>
#include<string>
#include"random.h"
#include<algorithm>

using namespace std;

struct Files{
    string radius;
    string points;
    string eq_acceptance;
    string run_acceptance;
};

Files get_Files(bool is100, bool gauss, double rate);
void is_open(ofstream& coutf, string filename);
double error(double av, double av2, int n);
void compute_prog_mean(const vector<double>& data, vector<double>& sum_prog, vector<double>& err_prog);
void save(const string filename, const vector<double>& sum_prog, vector<double>& err_prog);
double prob(array<double, 3>& x, int n, int l, int m);
array<double, 3> step_rw(Random& r, array<double, 3>& x_i, double a, bool Gauss);
bool metro(Random& rn, array<double, 3>& old_point, array<double, 3>& point, int n, int l, int m);
array<double, 3> sample(Random& rn, array<double, 3>& actual_point, double rate, int& accettanza, array<int, 3>& nlm, bool Gauss);


array<double, 3> operator+(const array<double, 3>& a, const array<double, 3>& b);