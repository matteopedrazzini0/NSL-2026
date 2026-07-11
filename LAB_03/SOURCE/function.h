#pragma once

#include<vector>
#include<string>

using namespace std;

double error(double av, double av2, int n);
void compute_prog_mean(const vector<double>& data, vector<double>& sum_prog, vector<double>& err_prog);
void save(const string filename, const vector<double>& sum_prog, vector<double>& err_prog);
double call_price(double S, double K);
double put_price(double S, double K);