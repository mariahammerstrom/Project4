/* ISING MODEL IN TWO DIMENSIONS (NO MAGNTIC FIELD)

   Calculates expectation values <E>, <|M|>, Cv, chi as functions of T for a square lattice with 2 possible
   spin values (+1, -1) using the Metropolis algorithm with periodic boundary conditions.

   Coupling constat J = 1
   Boltzmann constant k = 1

   The code is parallelized using Open MPI.

   All constants are set inside the program.
*/

#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <random>
#include <mpi.h>

using namespace std;

inline int periodic(int i, int limit, int add){return (i+limit+add)%(limit);}
void Metropolis(int,int **,double &,double &,double *,int&);
void initialize(int,int**,double&,double&);
void initialize_random(int n_spins,int **spin_matrix,double& E, double& M);
void analytical_cf(double T,double& Z,double& E_exp,double& M_exp,double& Cv,double& chi,int N);
void print(double E_exp,double M_abs,double Cv,double chi,int N);
void write_to_file(ofstream &file,int mc, double T,double *average,int accepted_configs, int n_spins);
void ExpectationValues_MC(double T,ofstream &file,ofstream &fileE,int n_spins,int mc,int**spin_matrix,bool &rand,double *w,int my_rank);
void ExpectationValues_T(double T,ofstream &file,int n_spins,int mc,int**spin_matrix,bool &rand,int myloop_begin,int myloop_end,int my_rank);
void print_MC(int mc, double T,double *average,int n_spins);


int main(int argc, char *argv[]) // Leave second argument blank
{
    // INITIAL CONDITIONS
    double T = 1.0;                                     // Temperature [kT/J]
    int n_spins = 20;                                   // Number of spins
    int N = n_spins*n_spins;                            // Lattice dimensions (square)
    int MC_cycles = 1000000;                            // Number of Monte Carlo cycles

    bool random = false; // true = random spin matrix, false = all spins pointing upwards

    double temp_step = 0.05;                            // Steps in temperature
    double initial_temp = 2.0;                          // Initial temperature
    double final_temp = 2.4;                            // Final temperature


    // MPI INITIALIZATIONS
    int my_rank,numprocs;
    MPI_Init (&argc, &argv);
    MPI_Comm_size (MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);

    int no_intervalls = MC_cycles/numprocs;
    int myloop_begin = my_rank*no_intervalls + 1;
    int myloop_end = (my_rank+1)*no_intervalls;
    if ((my_rank == numprocs-1) &&( myloop_end < MC_cycles)) myloop_end = MC_cycles;

    // Broadcast to all nodes common variables
    MPI_Bcast (&n_spins, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast (&initial_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&final_temp, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast (&temp_step, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if(n_spins == 2 && my_rank == 0){
        // ANALYTICAL SOLUTION: CLOSED-FORM
        cout << "ANALYTICAL SOLUTION: CLOSED-FORM" << endl;

        double Z_cf,E_exp_cf,M_exp_cf,C_V_cf,chi_cf; // Partition function, expectation values <E>, <M>, Cv, X
        Z_cf = E_exp_cf = M_exp_cf = C_V_cf = chi_cf = 0.0;

        analytical_cf(T,Z_cf,E_exp_cf,M_exp_cf,C_V_cf,chi_cf,N); // Calculate analytical solution
        print(E_exp_cf,M_exp_cf,C_V_cf,chi_cf,N); // Print results
    }


    // SPIN MATRIX
    int **spin_matrix;
    spin_matrix = new int*[n_spins];

    for(int i=0;i<n_spins;i++)
        spin_matrix[i] = new int[n_spins];

    for(int i=0;i<n_spins;i++){
        for(int j=0;j<n_spins;j++){
            spin_matrix[i][j] = 1;
        }
    }

    // METROPOLIS ALGORITHM
    if(my_rank == 0){
        cout << "METROPOLIS ALGORITHM" << endl;
        cout << "MC cycles" << "\t" << MC_cycles << endl;
        cout << "T" << "\t" << "\t" << T << endl;}

    // Set up array for possible energy changes
    double w[17];
    for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
    for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T);

    // Create files
    char filename_MC[1000];
    sprintf(filename_MC, "ExpectationValues_MC_%d_%.1f_%d.txt", n_spins, T, random);
    char filename_E[1000];
    sprintf(filename_E,"Energy_MC_%d_%.1f_%d.txt",n_spins,T,random);

    ofstream file_MC(filename_MC); // File for expectation values
    ofstream file_E(filename_E); // File for energies

    // Expectation values as a function # of MC cycles
    int numberOfMonteCarloCycles = myloop_end - myloop_begin;
    ExpectationValues_MC(T,file_MC,file_E,n_spins,numberOfMonteCarloCycles,spin_matrix,random,w,my_rank);

    file_MC.close();
    file_E.close();

    // Create file for temperature variations
    char filename_T[1000];
    sprintf(filename_T, "ExpectationValues_temp_%d_%d.txt", n_spins, random);
    ofstream file_T(filename_T);

    // Expectation values as a function of temperature variations
    for(T=initial_temp;T<=final_temp;T+=temp_step){
        cout << "T = " << T << endl;
        ExpectationValues_T(T,file_T,n_spins,MC_cycles,spin_matrix,random,myloop_begin,myloop_end,my_rank);}
    file_T.close();


    // END MPI
    MPI_Finalize ();


    // CLEAR MEMORY
    for(int i=0;i<n_spins;i++) {
        delete[] spin_matrix[i];
    }
    delete[] spin_matrix;

    return 0;
}

void analytical_cf(double T,double& Z,double& E_exp,double& M_exp,double& Cv,double& chi,int N){
    // Function to calculate the analytical solutions in the L = 2 case

    Z = (4*cosh(8/T) + 12); // Partition function
    double Z_inv_cf = 1./Z;
    E_exp = -8*tanh(8/T)/N;
    Cv = (8/T/cosh(8/T))*(8/T/cosh(8/T))/N;
    M_exp = 8*(exp(8/T) + 2)*Z_inv_cf/N;
    chi = (32*Z_inv_cf*(exp(8/T) + 1) - Z_inv_cf*Z_inv_cf*(8*exp(8/T) + 4)*(8*exp(8/T) + 4))/T/N;
}

void initialize(int n_spins,int **spin_matrix,double& E, double& M){
    // Function that initializes the spin matrix (ordered), magnetization and energy

    // Setup spin matrix and initial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            spin_matrix[y][x] = 1; // spin orientation for the ground state
            M +=  (double) spin_matrix[y][x];
        }
    }

    // Setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*(spin_matrix[periodic(y,n_spins,-1)][x] + spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}

void initialize_random(int n_spins,int **spin_matrix,double& E, double& M){
    // Function that initializes the spin matrix with randomly oriented spins, magnetization and energy

    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);

    // Setup spin matrix
    double random;
    for(int x=0;x<n_spins;x++){
        for(int y=0;y<n_spins;y++){
            random = distribution(generator);
            if(random<0.5) spin_matrix[x][y] = -1;
            else spin_matrix[x][y] = 1;
        }
    }
    // Setup intial magnetization
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            M +=  (double) spin_matrix[y][x];
        }
    }
    // Setup initial energy
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            E -=  (double) spin_matrix[y][x]*(spin_matrix[periodic(y,n_spins,-1)][x] + spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}


void Metropolis(int n_spins, int **spin_matrix, double& E, double&M, double *w, int& accepted_configs){
    // Function that implements the Metropolis algorithm.

    // Generate random number generator
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();    // generate seed
    default_random_engine generator(seed);                                          // start random number generator
    uniform_real_distribution<double> distribution(0.0,1.0);                        // pick type of random distribution

    // Loop over all spins
    for(int y=0; y<n_spins*n_spins; y++) {

        // Find random position
        int ix = (int) (distribution(generator)*(double)n_spins);
        int iy = (int) (distribution(generator)*(double)n_spins);
        int deltaE = 2*spin_matrix[iy][ix]*(spin_matrix[iy][periodic(ix,n_spins,-1)]+ spin_matrix[periodic(iy,n_spins,-1)][ix] + spin_matrix[iy][periodic(ix,n_spins,1)] + spin_matrix[periodic(iy,n_spins,1)][ix]);

        // Perform the Metropolis test
        if (distribution(generator) <= w[deltaE+8]) {
            spin_matrix[iy][ix] *= -1;                      // Flip one spin and accept new spin config

            // Update energy and magnetization
            M += (double) 2*spin_matrix[iy][ix];
            E += (double) deltaE;
            accepted_configs += 1;
        }
    }
}


void ExpectationValues_T(double T,ofstream &file,int n_spins,int mc,int**spin_matrix,bool &rand,int myloop_begin,int myloop_end,int my_rank){
    // A function that calculates expectation values as a function of MC cycles when looping over temperatures
    // and writes calculations to data.

    // Set up array for possible energy changes
    double w[17];
    for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
    for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T);

    // Initialize array for expectation values
    double average[5];
    double total_average[5];

    // Initialize sums
    double M = 0;
    double E = 0;
    int accepted_configs = 0; // Initialize count of accepted configurations

    // Initialize arrays
    for(int i=0;i<5;i++) average[i] = 0.0;
    for(int i=0;i<5;i++) total_average[i] = 0.0;

    // Setup spin matrix
    if(rand) initialize_random(n_spins,spin_matrix,E,M);
    else initialize(n_spins,spin_matrix,E,M);

    //for(int cycles=1; cycles <= mc;cycles++)
    for (int cycles = myloop_begin; cycles <= myloop_end; cycles++){
        Metropolis(n_spins,spin_matrix,E,M,w,accepted_configs);

        // Update expectation values
        average[0] += E; average[1] += E*E;
        average[2] += M; average[3] += M*M; average[4] += fabs(M);
    }

    // Find total average
    for(int i=0;i<5; i++){
        MPI_Reduce(&average[i], &total_average[i], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Write to file
    if (my_rank == 0) {
        write_to_file(file,mc,T,total_average,accepted_configs,n_spins);
    }

}

void ExpectationValues_MC(double T,ofstream &file,ofstream &fileE,int n_spins,int mc,int**spin_matrix,bool &rand,double *w,int my_rank){
    // A function that calculates expectation values as a function of MC cycles and writes calculations to data.

    // Initialize array for expectation values
    double average[5];

    // Initialize sums
    double M = 0;
    double E = 0;
    int accepted_configs = 0; // Initialize count of accepted configurations
    int countstart = 0;

    for(int i=0;i<5;i++) average[i] = 0.0;

    // Setup spin matrix
    if(rand) initialize_random(n_spins,spin_matrix,E,M);
    else initialize(n_spins,spin_matrix,E,M);

    bool count = false;
    bool first = true;
    double test, Eprev;
    int minimum_mcs;

    for(int cycles=1; cycles <= mc;cycles++){
        Eprev = E; // save previous energy value

        Metropolis(n_spins,spin_matrix,E,M,w,accepted_configs);
        if(count && my_rank == 0){
            fileE << E/n_spins/n_spins << endl; // Write energies to file
        }

        // Update expectation values
        average[0] += E; average[1] += E*E;
        average[2] += M; average[3] += M*M; average[4] += fabs(M);

        // Test energy difference to find minimum number of MC cycles before ground state is reached
        if(cycles > 1){
            test = fabs((Eprev-E)/Eprev);
            if (test < 0.01) countstart = 1;}

        if (countstart == 1 && first){
            minimum_mcs = cycles;
            first = false;
            count = true;
        }

        // Write to file
        if (my_rank == 0) {
            write_to_file(file,cycles,T,average,accepted_configs,n_spins);
        }
    }
    if (my_rank == 0){
        cout << "Min. # cycles " << "\t" << minimum_mcs << endl;
        print_MC(mc, T,average,n_spins);
    }
}


void write_to_file(ofstream &file,int mc, double T,double *average,int accepted_configs,int n_spins){
    // A function that writes expectation values to file

    double norm = 1/((double) (mc));
    double E_exp = average[0]*norm;
    double E_exp2 = average[1]*norm;
    double M_exp2 = average[3]*norm;
    double M_abs = average[4]*norm;
    double Cv = (E_exp2-E_exp*E_exp)/T/T;
    double chi = (M_exp2-M_abs*M_abs)/T;
    double sigmaE = (E_exp2 - E_exp*E_exp)/n_spins/n_spins;

    E_exp /= n_spins*n_spins;
    M_abs /= n_spins*n_spins;
    Cv /= n_spins*n_spins;
    chi /= n_spins*n_spins;

    // Write solution to file
    file << T << "\t" << mc << "\t" << E_exp << "\t" << M_abs << "\t" << Cv << "\t" << chi << "\t" << sigmaE << "\t" << accepted_configs*norm << endl;
}


void print(double E_exp,double M_abs,double Cv,double chi,int N){
    // A function that prints out expectation values to screen

    cout << "L" << "\t" << "\t" << N << endl;
    cout << "<E>/spins" << "\t" << "\t" << fixed << E_exp << endl;
    cout << "<M>/spins" << "\t" << "\t" << M_abs << endl;
    cout << "Cv/spins" << "\t" << "\t" << Cv << endl;
    cout << "chi/spins" << "\t" << "\t" << chi << endl << endl;
}

void print_MC(int mc, double T,double *average,int n_spins){
    // Functions that prints final results of the MC method

    double N = n_spins*n_spins;
    double norm = 1/((double) (mc));
    double E_exp = average[0]*norm;
    double E_exp2 = average[1]*norm;
    double M_exp2 = average[3]*norm;
    double M_abs = average[4]*norm;
    double Cv = (E_exp2-E_exp*E_exp)/T/T;
    double chi = (M_exp2-M_abs*M_abs)/T;

    print(E_exp/N,M_abs/N,Cv/N,chi/N,N);
}
