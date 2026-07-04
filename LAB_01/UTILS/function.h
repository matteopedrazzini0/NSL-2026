#pragma once

#include<vector>
#include<string>
#include"random.h"

using namespace std;

double error(double av, double av2, int n);
void compute_prog_mean(const vector<double>& data, vector<double>& sum_prog, vector<double>& err_prog);
void save(const string filename, const vector<double>& sum_prog, vector<double>& err_prog);
double compute_X2(const int M, const int n, const double x_min, const double x_max, Random& r);