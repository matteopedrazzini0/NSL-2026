#include"VMC_SA.h"
#include<chrono>

using namespace std;
 
int main (int argc, char *argv[]){
    auto start=chrono::high_resolution_clock::now();
    VMC_SA sy;
    sy.prepare();
    if(sy.get_type()==-1){
        for(int i=0; i<sy.get_Nbl(); i++){
            sy.H();
            sy.average(i+1);
            sy.block_reset(i+1);
        }
        sy.sample_pdf(100000);
    }
    if(sy.get_type() == 0){
        sy.optimize();
    }
    else if(sy.get_type() == 1){
        sy.reset_global();
        sy.Load_param("../OUTPUT/parameters_evolution.dat");
        sy.printmusigma();
        for(int i=0; i<sy.get_Nbl(); i++){
            sy.H();
            sy.average(i+1);
            sy.block_reset(i+1);
        }
        sy.sample_pdf(100000);
    }
    else if(sy.get_type() == 2){
        sy.optimize();
        for(int i=0; i<sy.get_Nbl(); i++){
            sy.H();
            sy.average(i+1);
            sy.block_reset(i+1);
        }
        sy.sample_pdf(100000);
    }
    sy.finalize();
    
    auto end=chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed=end-start;
    cout<<"Tempo di esecuzione: "<<elapsed.count() / 60<<" minuti\n";
    return 0;
}