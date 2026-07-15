#include<iostream>
#include<fstream>
#include<mpi.h>
#include"population.h"

using namespace std;

int main(int argc, char *argv[]){
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    double tstart = MPI_Wtime();

    int N_cities = 110;
    int N_individuals = 1;

    ofstream coutf("../OUTPUT/acceptance_"+to_string(rank)+"_evolution.dat");
    coutf << "#GEN" << setw(10) << "T" << setw(10) << "A" << endl;
    // Ogni rank svuota il proprio file all'avvio del programma
    ofstream inizialize_file("../OUTPUT/acceptance_rank_" + to_string(rank) + ".dat", ios::trunc);
    if(inizialize_file.is_open()){
        inizialize_file << "#BLOCK" << setw(10) << "PARTNER" << setw(10) << "MY_BETA" << setw(10) << "PART_BETA" << setw(10) << "DELTA_B" << "DELTA_E" << setw(10) << "ACCEPTED" << endl;
        inizialize_file.close();
    }
    
    population mypop(N_individuals, N_cities);
    mypop.initialize(rank, size);
    ofstream coutL(mypop.get_filename_loss());
    int partner;
    for(int gen_block=0; gen_block<mypop.get_n_migr(); gen_block++){
        // local evolution
        mypop.reset_acc();
        for(int gen=0; gen<mypop.get_step(); gen++){
            mypop.evolution(gen_block*mypop.get_step()+gen, coutL);
        }
        coutf << gen_block << setw(14) << 1.0 / mypop.get_beta() << setw(14) << mypop.get_acc() << endl;
        mypop.sort_pop(); // inutile perchè un solo elemento, lascio per successive
        // scambio informazioni
        partnership(rank, size, partner, gen_block%2);
        exchange(rank, partner, mypop, N_cities, coutL, gen_block);
    }
    coutL.close();
    // Now any rank sort is population based on loss and fill an array of the best path's loss sended to the firs rank
    mypop.sort_pop();
    double Loss = mypop.get_cost(0);
    double global_min;
    MPI_Allreduce(&Loss, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    // il rank che ha la funzione di costo più bassa di tutti stampa la soluzione finale, dovrebbe essere quello a T più bassa
    int ho_min = (Loss==global_min) ? 1 : 0;
    vector<int> who_have_min(size);
    int selected = -1;
    select_unique(rank, size, ho_min, who_have_min, selected);
    if(rank==selected) mypop.Print_final_Path("../OUTPUT/Miglior_percorso_prov_ita_"+to_string(rank)+".xy", rank);

    double tend = MPI_Wtime();
    cout << "Tempo di esecuzione del rank " << rank << ": " << (tend - tstart) / 60 << "min" << endl;
    MPI_Finalize();
    return 0;
}