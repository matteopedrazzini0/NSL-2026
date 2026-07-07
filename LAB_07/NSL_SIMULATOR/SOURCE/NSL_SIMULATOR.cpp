/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include<iostream>
#include"system.h"
#include<chrono>

using namespace std;

int main (int argc, char *argv[]){
  auto start=chrono::high_resolution_clock::now();

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  //equilibrazione
  int equi_step = 5000;
  for(int i=0; i<equi_step; i++){
    SYS.step();
  }
  for(int i=0; i<SYS.get_nbl(); i++){
    for(int j=0; j<SYS.get_nsteps(); j++){
      SYS.step();
      SYS.measure();
      if(j%50==0){
        //SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf++; 
      }
    }
    SYS.averages(i+1);
    SYS.block_reset(i+1);
  }
  SYS.finalize();

  auto end=chrono::high_resolution_clock::now();
  chrono::duration<double> elapsed=end-start;
  cout<<"Tempo di esecuzione: "<<elapsed.count() / 60<<" minuti\n";

  return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
