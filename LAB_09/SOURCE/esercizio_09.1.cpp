#include<iostream>
#include<chrono>
#include"population.h"

using namespace std;

int main(int argc, char *argv[]){
    auto start=chrono::high_resolution_clock::now();
    //double epsilon = 1e-4;
    //double best_dist = mypop.get_best_path();
    //int count = 0;
    //int max_count = 2000;
    //int gen = 0;
    int N_cities = 34;
    int N_individuals = 300;
    int n_gen_max = 1000;

    // ibrido, mi fermo se miglioro realmente o se raggiungo un max iter
    /*while(count < max_count && gen < n_gen_max){
        mypop.evolution();
        double act_dist = mypop.get_best_path();
        
        if(act_dist < best_dist - epsilon){
            best_dist = act_dist;
            count = 0;
        }
        else count++;
    }*/
    /* iterativo puro */
    population mypop(N_individuals, N_cities);
    mypop.initialize();
    for(int g=0; g<n_gen_max; g++){
        mypop.evolution(g);
    }/**/
    mypop.Print_final_Path();
    auto end=chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed=end-start;
    cout<<"Tempo di esecuzione: "<<elapsed.count() / 60 <<" minuti\n";
    return 0;
}