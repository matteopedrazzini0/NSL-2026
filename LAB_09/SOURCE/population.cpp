#include"population.h"
#include<fstream>
#include<algorithm>
#include<iterator>
#include<cmath>

using namespace std;

// SONO DEL LAB 9

void population :: initialize(){
    // set mutation probability
    // initialize the random generator
    int p1, p2; // Read from ../INPUT/Primes a pair of numbers to be used to initialize the RNG
    ifstream Primes("../INPUT/Primes");
    Primes >> p1 >> p2;
    Primes.close();

    int seed[4]; // Read the seed of the RNG
    ifstream Seed("../INPUT/seed.in");
    Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
    _rnd.SetRandom(seed,p1,p2);
    Seed.close();

    ifstream in("../INPUT/input.txt");
    string prop;
    while(in >> prop){
        if(prop=="TOPOLOGY") in >> _topology;
        else if(prop=="SIZE") in >> _r;
        // set mutation probability
        else if(prop=="P_SWAP") in >> P_swap;
        else if(prop=="P_SHIFT") in >> P_shift;
        else if(prop=="P_PERM") in >> P_perm;
        else if(prop=="P_INV") in >> P_inv;
        else if(prop=="P_CORSS") in >> P_cross;
        else if(prop=="ENDINPUT") break;
        else{
            cerr << "Unknow type" << endl;
            exit(1);
        }
    }
    in.close();

    // define the map
    if(_topology == "circo"){
        // set the correct distance
        //distance = [this](int a, int b){return distance_circ(a, b);};
        _filename_path = "../OUTPUT/Path_circo.xy";
        _filename_loss = "../OUTPUT/Loss_values_circo.dat";
        for(int i=0; i<N_genes; i++){
            double theta = _rnd.Rannyu(0, 2 * M_PI);
            map[i] = {_r * cos(theta), _r * sin(theta)};
        }
    }
    else if(_topology == "square"){
        //distance = [this](int a, int b){return distance_plane(a, b);};
        _filename_path = "../OUTPUT/Path_square.xy";
        _filename_loss = "../OUTPUT/Loss_values_square.dat";
        for(int i=0; i<N_genes; i++){
            double x = _rnd.Rannyu(-_r/2, _r/2);
            double y = _rnd.Rannyu(-_r/2, _r/2);
            map[i] = {x, y};
        }
    }
    else{
        cerr << "Invalid argument, unknown topology" << endl;
        exit(1);
    }

    // initialize the population
    for(int i=0; i<N_genes; i++){
        Pop[0].genes[i] = i;
    }
    Pop[0].cost = Cost_function(Pop[0].genes);
    for(int i=1; i<N_ind; i++){
        Pop[i].genes = Pop[0].genes;
        for(int j=1; j<N_genes-1; j++){ // shuffle genero un intero tra le posizioni non ancora toccate, e poi scambio i due,
            int k = int(_rnd.Rannyu(j, N_genes-1)); // così l'elemento già modificato non viene chiamato, metto (0, Ngenes) perchè rannyu genera escludendo l'ultimo elemento
            int temp = Pop[i].genes[j];
            Pop[i].genes[j] = Pop[i].genes[k];
            Pop[i].genes[k] = temp;
        }
        Pop[i].cost = Cost_function(Pop[i].genes);
    }
    // prepare the output file
    ofstream coutf;
    coutf.open(_filename_loss);
    coutf << "# EV_STEP" << setw(16) << "L_best" << setw(16) << "<L>" << endl;
    coutf.close();
}

void population :: sort_pop(){
    sort(Pop.begin(), Pop.end(), [](const Individual& a, const Individual& b){
                                    return a.cost < b.cost;});
}

Individual population :: selector(){ // Only if the population have been already sorted
    int j = int(N_ind * pow(_rnd.Rannyu(), 2.0)); // per penalizzare elementi in fondo all'array
    return Pop[j];
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
/*
bool population :: check_function(vector<int> chr){
    sort(chr.begin(), chr.end());
    auto y = unique(chr.begin(), chr.end()); // unique returns the pointer to the last elemets of the list, after compressed it with the unique elements. If the two list have different lenghts the pointer to the final element
    return y == chr.end();
}*/

void population :: crossover(Individual& selected1, Individual& selected2){
    vector<int> s1 = selected1.genes, s2 = selected2.genes; // progenie
    int cut = int(_rnd.Rannyu(1, N_genes));
    vector<bool> used_1(N_genes, false), used_2(N_genes, false); // vettori che mi dicono se una città è contenuta entro un centro indice
    // Potresti sistemarlo togliendo la prima città, facendo partire i cicli da 1 e mettendo gli used con una dimensione in meno (N_genes-1)
    for(int i=0; i<cut; i++){ // metto true dove trovo la città già vista
        used_1[selected1.genes[i]] = true;
        used_2[selected2.genes[i]] = true;
    }

    int index = cut;
    for(size_t i=0; i<selected2.genes.size(); i++){ // for(int gene : selected2)
        int gene = selected2.genes[i];
        if(!used_1[gene]) s1[index++] = gene;
    }

    index = cut;
    for(size_t i=0; i<selected1.genes.size(); i++){
        int gene = selected1.genes[i];
        if(!used_2[gene]) s2[index++] = gene;
    }

    selected1.genes = std::move(s1);
    selected2.genes = std::move(s2);
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

void population :: evolution(int ev_step){
    vector<Individual> new_Pop;
    new_Pop.reserve(N_ind);
    sort_pop();

    ofstream coutf;
    coutf.open(_filename_loss, ios::app);
    coutf << ev_step << setw(16) << Pop[0].cost << setw(16) << get_half_average() << endl;
    coutf.close();

    cloning();
    if(_rnd.Rannyu()<0.6) new_Pop.push_back(elitario); // elitismo con probabilità 60%
    while(new_Pop.size() < N_ind){ // parent è ambiguo come nome, ma in realtà viene modificato al suo interno
        Individual parent1 = selector();
        Individual parent2 = selector();
        if(_rnd.Rannyu()<P_cross){
            crossover(parent1, parent2);
        }
        if(_rnd.Rannyu()<P_swap) mutation(0, parent1);
        if(_rnd.Rannyu()<P_swap) mutation(0, parent2);
        
        if(_rnd.Rannyu()<P_shift) mutation(1, parent1);
        if(_rnd.Rannyu()<P_shift) mutation(1, parent2);
        
        if(_rnd.Rannyu()<P_perm) mutation(2, parent1);
        if(_rnd.Rannyu()<P_perm) mutation(2, parent2);
        
        if(_rnd.Rannyu()<P_inv) mutation(3, parent1);
        if(_rnd.Rannyu()<P_inv) mutation(3, parent2);
        
        parent1.cost = Cost_function(parent1.genes);
        parent2.cost = Cost_function(parent2.genes);
        if(check_function(parent1.genes) && new_Pop.size()<N_ind) new_Pop.push_back(parent1);
        if(check_function(parent2.genes) && new_Pop.size()<N_ind) new_Pop.push_back(parent2);        
    }
    Pop = std::move(new_Pop);
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
    return sqrt(pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2));
}

void population :: cloning(){ // Only if the population have been already sorted
    elitario = Pop[0];
}

void population :: Print_final_Path(){
    sort_pop();
    ofstream coutf;
    coutf.open(_filename_path);
    if(coutf.is_open()){
        coutf << "# City" << setw(16) << "x" << setw(16) << "y" << endl;
        for(int i=0; i<N_genes; i++){
            int city = Pop[0].genes[i];
            coutf << city << setw(16) << map[city][0] << setw(16) << map[city][1] << endl;
        }
    }
    else cerr << "Problem: Unable to open file " << _filename_path << endl;
    coutf.close();
    cout << "Best cost function: " << Pop[0].cost << endl;
    _rnd.SaveSeed();
}


/*vecchia implementaziopne di mutation 1
        vector<int> copy1(group_size);
        vector<int> copy2(Ngenes-group_size);
        for(int i=0; i<group_size; i++){
            copy1[i] = selected[i+1];
        }
        for(int j=0; j<Ngenes-group_size-1; j++){
            copy2[j] = selected[j+group_size+1];
        }
        for(int i=0; i<Ngenes-1; i++){
            selected[i+1] = (i<Ngenes-group_size-1) ? copy2[i] : copy1[i-(Ngenes-group_size-1)];
        }
        */

/* vecchia implementazione di mutation 3
        int fermo = (m%2 == 0) ? m / 2 : (m-1) / 2;
        for(int i=1; i<=fermo; i++){
            int temp = selected[m-i+1];
            selected[m-i+1] = selected[i];
            selected[i] = temp;
        }
        */