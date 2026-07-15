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
    int N_cities, N_individuals, N_ev_per_migr, N_migr, n_migrants;

    ifstream inpu("../INPUT/settings.txt");
    string prop;
    while(inpu>>prop){
        if(prop=="N_CITIES") inpu >> N_cities;
        else if(prop=="N_INDIVIDUALS") inpu >> N_individuals;
        else if(prop=="N_EV_PER_MIGR") inpu >> N_ev_per_migr;
        else if(prop=="N_MIGR") inpu >> N_migr;
        else if(prop=="N_MIGRANTS") inpu >> n_migrants;
        else if(prop=="ENDINPUT") break;
        else{
            cerr << "Unknow input" << endl;
            exit(1);
        }
    }

    vector<int> best_path(N_cities);
    vector<int> recv_path(N_cities);
    vector<int> migration(size);
    
    population mypop(N_individuals, N_cities);
    mypop.initialize(rank);
    ofstream coutf(mypop.get_filename_loss());
    for(int gen_block=0; gen_block<N_migr; gen_block++){
        // local evolution
        for(int gen=0; gen<N_ev_per_migr; gen++){
            mypop.evolution(gen_block*N_ev_per_migr+gen, coutf);
        }
        if(rank==0){
            if(gen_block%10==0) cout << "[rank 0] Avanzamento: " << gen_block << " blocchi" << endl;
        }
        mypop.sort_pop();
        // migrations
        int next_partner, prev_partner;
        partnership(rank, size, mypop, migration, next_partner, prev_partner);
        for(int l=0; l<n_migrants; l++){
            best_path = mypop.get_path(l).genes;
            // Sending messages
            MPI_Sendrecv(&best_path[0], N_cities, MPI_INT, next_partner, 0,
                         &recv_path[0], N_cities, MPI_INT, prev_partner, 0,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            mypop.substitute(recv_path, N_individuals-1-l);
        }
    }
    coutf.close();
    // Now any rank sort is population based on loss and fill an array of the best path's loss sended to the firs rank
    mypop.sort_pop();
    double Loss = mypop.get_cost(0);
    double global_min;
    MPI_Allreduce(&Loss, &global_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    // il rank che ha la funzione di costo più bassa di tutti stampa la soluzione finale
    // Unicità di stampa
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