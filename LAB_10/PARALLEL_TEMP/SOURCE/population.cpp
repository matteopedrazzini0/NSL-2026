#include"population.h"
#include<fstream>
#include<algorithm>
#include<iterator>
#include<cmath>
#include<mpi.h>

using namespace std;

void population :: initialize(int rank, int size){
    // initialize the random generator
    _rank = rank;
    _size = size;
    int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("../INPUT/Primes");
    for(int i=0; i<=rank; i++){
        Primes >> p1 >> p2;
    }
    Primes.close();

    int seed[4]; // Read the seed of the RNG
    ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed,p1,p2);
    Seed.close();

    ifstream in("../INPUT/input.txt");
    string prop;
    while(in>>prop){
        if(prop=="MAX_TEMPERATURE"){
            double T;
            in >> T;
            _beta_max = 1.0 / T;
        }
        else if(prop=="MIN_TEMPERATURE"){
            double T;
            in >> T;
            _beta_min = 1.0 / T;
        }
        else if(prop=="TYPE") in >> _type;
        else if(prop=="RADIUS") in >> _r;
        else if(prop=="N_STEP_T_MAX") in >> _N_0;
        else if(prop=="N_MIGR") in >> _N_migr;
        else if(prop=="N_MUT_PER_EV") in >> _N_mut;
        else if(prop=="ENDINPUT") break;
        else{
            cerr << "Unknow type" << endl;
            exit(1);
        }

    }
    in.close();
    // define the map
    if(_type == "circo"){
        for(int i=0; i<N_genes; i++){
            double theta = _rnd.Rannyu(0, 2 * M_PI);
            map[i] = {_r * cos(theta), _r * sin(theta)};
        }
    }
    else if(_type == "square"){
        for(int i=0; i<N_genes; i++){
            double x = _rnd.Rannyu(-_r/2, _r/2);
            double y = _rnd.Rannyu(-_r/2, _r/2);
            map[i] = {x, y};
        }
    }
    else{
        this->read_map("../INPUT/"+_type); // in questo caso type sarà un filename
        _type = "cap_prov_ita_";
    }

    // initialize the population
    for(int i=0; i<N_genes; i++){
        Pop[0].genes[i] = i;
    }
    Pop[0].cost = Cost_function(Pop[0].genes);
    for(int i=1; i<N_ind; i++){
        Pop[i].genes = Pop[0].genes;
        for(int j=1; j<N_genes-1; j++){ // shuffle genero un intero tra le posizioni non ancora toccate, e poi scambio i due,
            int k = int(_rnd.Rannyu(j, N_genes)); // così l'elemento già modificato non viene chiamato, metto (0, Ngenes) perchè rannyu genera escludendo l'ultimo elemento
            int temp = Pop[i].genes[j];
            Pop[i].genes[j] = Pop[i].genes[k];
            Pop[i].genes[k] = temp;
        }
        Pop[i].cost = Cost_function(Pop[i].genes);
    }
    // prepare the output file
    _filename_loss = "../OUTPUT/Cost_values_" + _type + to_string(rank) + ".dat";
    ofstream coutf;
    coutf.open(_filename_loss);
    coutf << "# EV_STEP" << setw(16) << "L_best" << setw(16) << "<L>" << endl;
    coutf.close();

    _nattempts = 0;
    _naccepted = 0;

    double beta[11] = {10, 7.84974, 5.63345, 3.95345, 2.73345, 1.88345, 1.33345, 1.13349, 0.96394, 0.82365, 0.714286};
    _beta = beta[rank];
    //_beta = _beta_max + (_beta_min - _beta_max) * pow(1 - rank / double(size-1), 2.5); // Geometrica
    //_beta = _beta_min +(_beta_max-_beta_min)*double(rank) / double(size-1); // Lineare
    _N_moves = _N_0; //int(_N_0 * pow(_beta / _beta_max, 0.5));

    cout << "Sono il rank " << rank << " e ho un beta di: " << _beta << " E un numero di step pari a: " << _N_moves << endl;
}

void population :: read_map(string filename){
    ifstream in;
    in.open(filename);
    if(in.is_open()){
        double x, y;
        for(int i=0; i<N_genes; i++){
            in >> x >> y;
            map[i][0] = x;
            map[i][1] = y;
        }
    }
    else cerr << "PROBLEM: Unable to open INPUT file: " << filename << endl;
    in.close();
}

void population :: sort_pop(){
    if(N_ind!=1){
        sort(Pop.begin(), Pop.end(), [](const Individual& a, const Individual& b){
                                    return a.cost < b.cost;});
    }
    else return;
}

Individual& population :: selector(){ // Only if the population have been already sorted
    if(N_ind==1) return Pop[0];
    else{
        int j = int(N_ind * pow(_rnd.Rannyu(), 2.5)); // per penalizzare elementi in fondo all'array
        return Pop[j];
    }
}

bool population :: check_function(const vector<int>& chr){
    vector<bool> seen(N_genes, false);
    for(int j=0; j<N_genes; j++){
        int city = chr[j];
        if(city<0 || city>=N_genes) return false;
        if(seen[city]) return false;
        seen[city] = true;
    }
    return true;
}

void population :: mutation(int mut_type, Individual& selected){ // every mutation implemented is to considering always form the second element,
                                           // introduce a starting index
    if(mut_type == 0){
        int i{}, j{};
        do{
            i = int(_rnd.Rannyu(1, N_genes));
            j = int(_rnd.Rannyu(1, N_genes));
        }while(i == j);
        int temp{selected.genes[i]};
        selected.genes[i] = selected.genes[j];
        selected.genes[j] = temp;
    }
    else if(mut_type == 1){ // Shift blocks
        int group_size = int(_rnd.Rannyu(1, N_genes-1));
        int shift = int(_rnd.Rannyu(1, N_genes-1-group_size));
        // si può fare con rotate, cioè prendo una porzione da uno a groupsize e la scambio con quella che parte da groupsize+1 incluso e termina in groupsize+shift +1 incluso
        rotate(selected.genes.begin()+1, selected.genes.begin()+group_size+1, selected.genes.begin()+group_size+shift+1);
    }
    else if(mut_type == 2){ // Permutazion of two blocks of size m
        int m = int(_rnd.Rannyu(1, N_genes / 2.0));
        //int start = int(_rnd.Rannyu(1, N_genes - 2*m + 1)); // Per randomizzare anche il punto di partenza, e poi metto un + start dopo begin
        rotate(selected.genes.begin()+1, selected.genes.begin()+m+1, selected.genes.begin()+2*m+1);
    }
    else if(mut_type == 3){ // inverting
        int m = int(_rnd.Rannyu(1, N_genes-1));
        reverse(selected.genes.begin()+1, selected.genes.begin()+m+1);
    }
}

void population :: evolution(int ev_step, ofstream& coutf){

    if(ev_step%_N_moves==0) coutf << ev_step << setw(16) << Pop[0].cost << endl;
    Individual& parent = selector();
    Individual copy = parent;
    double old_loss = parent.cost;

    for(int j=0; j<_N_mut; j++) mutation(int(4 * pow(_rnd.Rannyu(), 2.5)), parent);
    
    double new_loss = Cost_function(parent.genes);
    double ratio = exp(-_beta * (new_loss - old_loss));
    if(check_function(parent.genes) && this->metro(ratio)){
        _naccepted++;
        parent.cost = new_loss;
    }
    else parent = copy;
    _nattempts++;
}

bool population :: metro(double ratio){
    double A = min(1.0, ratio);
    return _rnd.Rannyu()<=A;
}

double population :: Cost_function(const vector<int>& chromo){
    double L= 0.0;
    for(int i=0; i<N_genes; i++){
        int city1 = chromo[i];
        int city2 = (i != N_genes-1) ? chromo[i+1] : chromo[0];
        double dist = distance(city1, city2);
        L += dist;
    }
    return L;
}

double population :: get_cost(int index){
    return Pop[index].cost;
}

double population :: get_half_average(){
    double sum = 0.0;
    for(int i=0; i<N_ind/2; i++){
        sum += Pop[i].cost;
    }
    return sum / (N_ind * 0.5);
}

double population :: distance(int city1, int city2){
    vector<double> x1 = map[city1];
    vector<double> x2 = map[city2];
    double dx = x1[0] - x2[0];
    double dy = x1[1] - x2[1];
    return sqrt(dx*dx + dy*dy);
}

void population :: cloning(){ // Only if the population have been already sorted
    elitario = Pop[0];
}

void population :: Print_final_Path(string filename, int rank){
    sort_pop();
    ofstream coutf;
    coutf.open(filename);
    if(coutf.is_open()){
        coutf << "# City" << setw(16) << "x" << setw(16) << "y" << endl;
        for(int i=0; i<N_genes; i++){
            int city = Pop[0].genes[i];
            coutf << city << setw(16) << map[city][0] << setw(16) << map[city][1] << endl;
        }
    }
    else cerr << "Problem: Unable to open file " << filename << endl;
    coutf.close();
    cout << "Best cost function from rank " << rank << ": " << Pop[0].cost << endl;
    _rnd.SaveSeed();
}

Individual population :: get_path(int idx){
    return Pop[idx];
}

void population :: substitute(vector<int>& path, int idx){
    Pop[idx].genes = path;
    Pop[idx].cost = Cost_function(path);
}

double population :: rannyu(double a, double b){
    return _rnd.Rannyu(a, b);
}

void population :: take_path(vector<int>& path, string filename){
    ofstream coutf(filename);
    coutf << "City" << setw(16) << "x" << setw(16) << "y" << endl;
    if(coutf.is_open()){
        for(int i = 0; i<N_genes; i++){
            int city = path[i];
            coutf << city << setw(16) << map[city][0] << setw(16) << map[city][1] << endl;
        }
    }
    else cerr << "Problem: Unable to open file " << filename << endl;
    coutf.close();
    cout << "Best path cost function's: " << Cost_function(path);
}

void population :: reset_acc(){
    _nattempts = 0;
    _naccepted = 0;
}

double population :: get_acc(){
    return double(_naccepted) / double(_nattempts);
}

double population :: get_beta(){
    return _beta;
}

int population :: get_step(){
    return _N_moves;
}

int population :: get_n_migr(){
    return _N_migr;
}

string population :: get_filename_loss(){
    return _filename_loss;
}

void partnership(int rank, int size, int& part, int tag){
    if(tag==0){ // 0-1 / 2-3 / 4-5 / ...
        if(rank%2==0){ // Se sei pari il tu partner è il successivo
            part = (rank+1 < size) ? rank+1 : rank; // l'ultimo rimane accoppiato a se stesso e nel main salterà lo scambio
        }
        else part = rank-1;
    }
    else{ // 1-2 / 3-4 / ...
        if(rank==0) part = rank;
        else if(rank%2==1){ // Se sei dispari il tu partner è il successivo
            part = (rank+1 < size) ? rank+1 : rank;
        }
        else part = rank-1;
    }
}

void exchange(int rank, int partner, population& pop, int Ncities, ofstream& coutll, int genblock){ // da aggiungere come parametri int& nattemptss, int& nacceptedd,
    // Controllo che i rank disaccoppiati non facciano nulla
    if(rank==partner) return;
    double mybeta = pop.get_beta();
    double mycost = pop.get_cost(0);

    double beta_partner;
    double cost_partner;
    MPI_Sendrecv(&mybeta, 1, MPI_DOUBLE, partner, 0, &beta_partner, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Sendrecv(&mycost, 1, MPI_DOUBLE, partner, 1, &cost_partner, 1, MPI_DOUBLE, partner, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    bool accepted = false;
    if(rank<partner){ // il rank minore decide se lo scambio si fa
        double delta = (mybeta-beta_partner) * (mycost-cost_partner);
        accepted = pop.metro(exp(delta));
        MPI_Send(&accepted, 1, MPI_CXX_BOOL, partner, 2, MPI_COMM_WORLD);
    }
    else{
        MPI_Recv(&accepted, 1, MPI_CXX_BOOL, partner, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if(accepted){
        vector<int> best_path = pop.get_path(0).genes;
        vector<int> recv_path(Ncities);

        MPI_Sendrecv(&best_path[0], Ncities, MPI_INT, partner, 3, &recv_path[0], Ncities, MPI_INT, partner, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        pop.substitute(recv_path, 0);
    }
    pop.sort_pop();
    coutll << genblock*pop.get_step()+pop.get_step() - 1 << setw(16) << pop.get_cost(0) << endl;
    // Alla fine di replica_exchange, subito prima di uscire:
    ofstream cout_partnership("../OUTPUT/acceptance_rank_" + to_string(rank) + ".dat", ios::app);
    if(cout_partnership.is_open()){
        // Salviamo: BloccoCorrente, Partner, BetaMio, BetaPartner, TentativoAccettato(0 o 1)
        int current_block = genblock; // o come tieni traccia del blocco corrente
        cout_partnership << current_block << setw(10)
                         << partner << setw(10)
                         << mybeta << setw(10)
                         << beta_partner << setw(10)
                         << mybeta - beta_partner << setw(10)
                         << mycost-cost_partner << setw(10)
                         << (accepted ? 1 : 0) << endl;
    cout_partnership.close();
    }
}

void select_unique(int rank, int size, int have_min, vector<int>& who_have_mini, int& select){
    MPI_Gather(&have_min, 1, MPI_INT, &who_have_mini[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(rank==0){
        for(int i=0; i<size; i++){
            if(who_have_mini[i]==1){
                select = i;
                break;
            }
        }
    }
    MPI_Bcast(&select, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
