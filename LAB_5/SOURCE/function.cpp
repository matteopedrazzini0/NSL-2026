#include<iostream>
#include<fstream>
#include"random.h"
#include<cmath>
#include<algorithm>
#include"function.h"

double error(double av, double av2, int n){
    return (n==0) ? 0 : sqrt((av2 - pow(av, 2)) / n);
}

array<double, 3> operator+(const array<double, 3>& a, const array<double, 3>& b){
    return {a[0]+b[0], a[1]+b[1], a[2]+b[2]};
}

void compute_prog_mean(const vector<double>& data, vector<double>& sum_prog, vector<double>& err_prog){
    int n_Blocks = int(sum_prog.size());
    int n_measure = int(data.size());
    int L = n_measure / n_Blocks; // Block size
    vector<double> ave(n_Blocks, 0.0), ave2(n_Blocks, 0.0), sum_prog2(n_Blocks, 0.0);
    double cum_sum = 0.0, cum_sum2 = 0.0;
    
    for(int i=0; i<n_Blocks; i++){
        double sum = 0;
        for(int j=0; j<L; j++){
            sum += data[i*L+j];
        }
        ave[i] = sum/L;
        ave2[i] = ave[i]*ave[i];

        cum_sum += ave[i];
        cum_sum2 += ave2[i];

        sum_prog[i] =cum_sum / (i+1);
        sum_prog2[i] = cum_sum2 / (i+1);
        err_prog[i] = error(sum_prog[i], sum_prog2[i], i);
    }
}

void save(const string filename, const vector<double>& sum_prog, vector<double>& err_prog)
{
    ofstream wr;
    wr.open(filename);
    if(wr.is_open()){
        wr << "#blocchi" << " " << "mean_prog" << " " << "err_prog" << "\n";
        for(size_t i=0; i<sum_prog.size(); i++){
            wr << i+1 << " " << sum_prog[i] << " " << err_prog[i] << "\n";
        }
    }else cerr << "Problem to open " << filename << endl;
   wr.close();
}

// Genero un rw con una certa probabilità di fare un determinato passo, con probabilità simmetrica
array<double, 3> step_rw(Random& r, array<double, 3>& x_i, double a, bool Gauss){
    array<double, 3> incr;
    if(!Gauss){
        double theta = r.Rannyu(0, M_PI);
        double phi = r.Rannyu(0, 2*M_PI);
        incr = {a * cos(phi) * sin(theta), a * sin(phi) * sin(theta), a * cos(theta)};
    }
    else{
        double x{r.Gauss(0, 1)}, y{r.Gauss(0, 1)}, z{r.Gauss(0, 1)};
        incr = {a*x, a*y, a*z};
    }
    return x_i + incr;
}

array<double, 3> sample(Random& rn, array<double, 3>& actual_point, double rate, int& accettanza, array<int, 3>& nlm, bool Gauss){
    int n{nlm[0]}, l{nlm[1]}, m{nlm[2]};
    array<double, 3> new_point = step_rw(rn, actual_point, rate, Gauss);
    if(metro(rn, actual_point, new_point, n, l, m)){
        accettanza++;
        return new_point;
    }
    return actual_point;
}

bool metro(Random& rn, array<double, 3>& old_point, array<double, 3>& point, int n, int l, int m){
    double A = min(1.0, prob(point, n, l, m)/prob(old_point, n, l, m));
    return rn.Rannyu() <= A; // perchè se A==1 è automaticamnete accettato, dato che Rannyu() genera fino a 1 escluso
}

double prob(array<double, 3>& x, int n, int l, int m){
    double r2{x[0]*x[0] + x[1]*x[1] + x[2]*x[2]}, r{sqrt(r2)};
    if(n == 1 and l == 0 and m == 0){
        return exp(-2 * r); //the factor r2 comes from the integration measure of the radial probability, that is 4 pi r^2
    }
    else if(n == 2 and l == 1 and m == 0){
        return exp(-r) * x[2] * x[2];
    }
    else{
        return exp(-r/n); //* pow(r, l) * pow(assoc_laguerre(n-l-1, 2*l+1, 2*r/n), 2) * pow(assoc_legendre(l, m, x[2] / r), 2);
    }
}

Files get_Files(bool is100, bool gauss, double rate){
    string orb = is100 ? "100" : "210";
    string tag = gauss ? "_Gauss" : "";
    int scaled = int(round(rate * 10));
    int int_part = scaled / 10;
    int dec_part = scaled % 10;
    string r = to_string(int_part) + "." + to_string(dec_part);
    return{
        "../OUTPUT/mean_r_"+orb+tag+".dat",
        "../OUTPUT/points_"+orb+tag+".xyz",
        "../OUTPUT/equil_acc_"+orb+"_rate_"+r+tag+".dat",
        "../OUTPUT/run_acc_"+orb+"_rate_"+r+tag+".dat"
    };
}

void is_open(ofstream& coutf, string filename){
    if(!coutf.is_open())
    {
        cerr << "Problem to open " << filename << endl;
        exit(1);   
    }
}