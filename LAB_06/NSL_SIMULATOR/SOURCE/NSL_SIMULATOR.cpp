/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include"system.h"

using namespace std;

int main (int argc, char *argv[]){

  int nconf = 1;
  System SYS;
  SYS.initialize();
  SYS.initialize_properties();
  SYS.block_reset(0);

  //equilibrazione
  int initial_step = SYS.get_nsteps();
  int equi_step = initial_step / 2.0;
  double T_min = 0.5, T_max = 3.0;
  // cambio temperature
  for(double T=T_min; T<=T_max; T+=0.02){
    SYS.reset_global();
    SYS.setTemp(T);
    for(int i{}; i<equi_step; i++){
      SYS.step();
    }
    SYS.Reset_acceptance();
    // evolvo
    SYS.set_step(int(initial_step * pow(T_max/T, 1.2)));
    for(int i=0; i < SYS.get_nbl(); i++){ //loop over blocks
      for(int j=0; j < SYS.get_nsteps(); j++){ //loop over steps in a block
        SYS.step();
        SYS.measure();
        if(j%50 == 0){
        // SYS.write_XYZ(nconf); //Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
          nconf++;
        }
      }
      SYS.averages(i+1);
      SYS.block_reset(i+1);
    }
  }
  SYS.finalize();
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
