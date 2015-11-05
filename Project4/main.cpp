/* ISING MODEL IN TWO DIMENSIONS (NO MAGNTIC FIELD)
 *
 * Calculates <E>, <|M|>, Cv, chi as functions of T for a square lattice with 2 possible spin values (+1, -1)
 * using the Metropolis algorithm.
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
<<<<<<< Updated upstream
#include <chrono>
#include <vector>
=======
//#include <chrono>
>>>>>>> Stashed changes
//#include <mpi.h>

using namespace std;

inline int periodic(int i, int limit, int add){return (i+limit+add)%(limit);}
void Metropolis(int,int **,double &,double &,double *,int&);
void initialize(int,int**,double&,double&);
void initialize_random(int n_spins,int **spin_matrix,double& E, double& M);
void ExpectationValues_toFile(double,ofstream&,ofstream&,int,int,int**,bool&,bool &rand,bool &count,double&,double&,double&,double&);
void analytical_cf(double T,double& Z,double& E_exp,double& M_exp,double& Cv,double& chi);

int main()
{
    // CONSTANTS
    double T = 1.0;                                     // Temperature [kT/J]
    int n_spins = 2;                                   // Number of spins
    int N = n_spins*n_spins;                            // Lattice dimensions (square)
    int MC_cycles = 10000;                                // Number of Monte Carlo cycles

    double temp_step = 0.5;                             // Steps in temperature
    double initial_temp = 1.0;                          // Initial temperature
    double final_temp = 3.0;                            // Final temperature


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


    // ANALYTICAL SOLUTION: CLOSED-FORM
    cout << "ANALYTICAL SOLUTION: CLOSED-FORM" << endl;

    double Z_cf,E_exp_cf,M_exp_cf,C_V_cf,chi_cf;
    Z_cf = E_exp_cf = M_exp_cf = C_V_cf = chi_cf = 0.0;

    analytical_cf(T,Z_cf,E_exp_cf,M_exp_cf,C_V_cf,chi_cf);

    cout << "L" << "\t" << N << endl;
    cout << "Z" << "\t" << Z_cf/N << endl;
    cout << "<E>/spins" << "\t" << E_exp_cf/N << endl;
    cout << "<M>/spins" << "\t" << M_exp_cf/N << endl;
    cout << "Cv/spins" << "\t" << C_V_cf/N << endl;
    cout << "chi/spins" << "\t" << chi_cf/N << endl << endl;


    // METROPOLIS ALGORITHM
    cout << "METROPOLIS ALGORITHM" << endl;
    cout << "MC cycles" << "\t" << "\t" << MC_cycles << endl;

    // Expectation values as a function of MC cycles
    bool first = true;
    bool random = true; // True = Random lattice, False = Fixed lattice
    bool count = false;

<<<<<<< Updated upstream
    double E_exp,M_abs,Cv,chi;
    E_exp = M_abs = Cv = chi = 0.0;

    if (random){
        ofstream file_MC("ExpectationValues_MC_random_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
        ofstream file_E("Energy_MC_random_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
        for(int MC = 1; MC<=MC_cycles; MC++)
            ExpectationValues_toFile(T,file_MC,file_E,n_spins,MC,spin_matrix,first,random,count,E_exp,M_abs,Cv,chi);
        file_MC.close();
        file_E.close();
=======
    /*
    if (random){
        ofstream file_MC("ExpectationValues_MC_random_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
        ofstream file_E("Energy_MC_random_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
>>>>>>> Stashed changes
    }
    else{
        ofstream file_MC("ExpectationValues_MC_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
        ofstream file_E("Energy_MC_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
<<<<<<< Updated upstream
        for(int MC = 1; MC<=MC_cycles; MC++)
            ExpectationValues_toFile(T,file_MC,file_E,n_spins,MC,spin_matrix,first,random,count,E_exp,M_abs,Cv,chi);
        file_MC.close();
        file_E.close();
    }

    //ofstream file_MC("ExpectationValues_MC_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
    //ofstream file_E("Energy_MC_" + to_string(n_spins) + "_" + to_string(T) + ".txt");


    cout << "L" << "\t" << "\t" << N << endl;
    cout << "<E>/spins" << "\t" << "\t" << E_exp << endl;
    cout << "<M>/spins" << "\t" << "\t" << M_abs << endl;
    cout << "Cv/spins" << "\t" << "\t" << Cv << endl;
    cout << "chi/spins" << "\t" << "\t" << chi << endl << endl;
=======
    }
    */


    ofstream file_E("Energy_MC_" + to_string(n_spins) + "_" + to_string(T) + ".txt");
    double E_exp,M_abs,Cv,chi;
    E_exp = M_abs = Cv = chi = 0.0;


    ofstream file_MC("ExpectationValues_MC_" + to_string(n_spins) + "_" + to_string(T) + ".txt");

    //for(int MC = 1; MC<=MC_cycles; MC++)
    for(int MC = 1; MC<=MC_cycles; MC+=4)
        ExpectationValues_toFile(T,file_MC,file_E,n_spins,MC,spin_matrix,first,random,count,E_exp,M_abs,Cv,chi);
    file_MC.close();
    file_E.close();
>>>>>>> Stashed changes


    cout << "L" << "\t" << "\t" << N << endl;
    cout << "<E>/spins" << "\t" << "\t" << E_exp << endl;
    cout << "<M>/spins" << "\t" << "\t" << M_abs << endl;
    cout << "Cv/spins" << "\t" << "\t" << Cv << endl;
    cout << "chi/spins" << "\t" << "\t" << chi << endl << endl;

<<<<<<< Updated upstream
    //for(T=initial_temp; T<=final_temp ; T+= temp_step)
        //ExpectationValues_toFile(T,file_T,n_spins,MC_cycles,spin_matrix,first);
    //file_T.close();
=======

    /*
    // Expectation values as a function of temperature variations
    ofstream file_T("ExpectationValues_temp_" + to_string(n_spins) + "_" + to_string(T) + ".txt");

    for(T=initial_temp; T<=final_temp ; T+= temp_step)
        ExpectationValues_toFile(T,file_T,file_E,n_spins,MC_cycles,spin_matrix,first,random,count,E_exp,M_abs,Cv,chi);
    file_T.close();
    */
>>>>>>> Stashed changes


    // Probability P(E_i), plot(?) comparing to sigmaE

    //ofstream fileP("ProbabilityL" + to_string(L) + ".txt"); // File for expectation values
    //fileP << L << "\t" << Probability << "\t" << sigmaE << endl; // Write solution to file
    //fileP.close();

    // Ecpectation values for different L and T. Parallelization!
    //ofstream fileL("ExpectationValuesL" + to_string(n_spins) + ".txt"); // File for expectation values
    //fileL << T << "\t" << E_exp << "\t" << M_exp << "\t" << Cv << "\t" << chi << endl;
    //fileL.close();

    // Estimate critical temperature
    return 0;
}

void analytical_cf(double T,double& Z,double& E_exp,double& M_exp,double& Cv,double& chi){
    Z = 4*cosh(8/T) + 12;
    double Z_inv_cf = 1./Z;
    E_exp = -8*tanh(8/T);
    Cv = (8/T/cosh(8/T))*(8/T/cosh(8/T));
    M_exp = 8*(exp(8/T) + 2)*Z_inv_cf;
    chi = (32*Z_inv_cf*(exp(8/T) + 1) - Z_inv_cf*Z_inv_cf*(8*exp(8/T) + 4)*(8*exp(8/T) + 4))/T;
}

void initialize(int n_spins,int **spin_matrix,double& E, double& M){
  // Setup spin matrix and intial magnetization
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
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0.0,1.0);
    double random;
    for(int x=0;x<n_spins;x++){
        for(int y=0;y<n_spins;y++){
            random = distribution(generator);
            if(random<0.5) spin_matrix[x][y] = -1;
            else spin_matrix[x][y] = 1;
        }
    }
  // Setup spin matrix and intial magnetization
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


void ExpectationValues_toFile(double T,ofstream &file,ofstream &fileE,int n_spins,int mc,int**spin_matrix,bool &first,bool &rand,bool &count,double &E_exp,double &M_abs,double &Cv, double &chi)
{
    double test;

    // Set up array for possible energy changes
    double w[17];
    for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
    for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T);

    // Initialize array for expectation values
    double average[5];

    // Initialize sums
    double M = 0;
    double E = 0;
    int accepted_configs = 0; // Initialize count of accepted configurations
    int countstart = 0;
    vector<double> energy;

    for(int i=0;i<5;i++) average[i] = 0;

    if(rand) initialize_random(n_spins,spin_matrix,E,M);
    else initialize(n_spins,spin_matrix,E,M);

    for(int cycles=1;cycles <= mc;cycles++){
        Metropolis(n_spins,spin_matrix,E,M,w,accepted_configs);
<<<<<<< Updated upstream
        if(count){
            energy.push_back(E); // Adds an element with value E to the end
            //fileE << E/n_spins/n_spins << endl;
=======
        //if(count){
        //    fileE << E/n_spins/n_spins << endl;
>>>>>>> Stashed changes
        }

        // Update expectation values
        //double Eprev = average[0];
        average[0] += E; average[1] += E*E;
        average[2] += M; average[3] += M*M; average[4] += fabs(M);

        //test = fabs((Eprev-average[0])/Eprev);
        //if (test < 0.05) countstart = 1;
    //}

<<<<<<< Updated upstream
    // E_tot = energy.size(); // Number of elements

    if (countstart == 1 && first){
        cout << "Min. # cycles " << "\t" << mc << endl;
        first = false;
        count = true;
    }
=======
    /*if (countstart == 1 && first){
        cout << "Min. # cycles " << "\t" << mc << endl;
        first = false;
        count = true;
    }*/
>>>>>>> Stashed changes

    double norm = 1/((double) (mc));
    E_exp = average[0]*norm;
    double E_exp2 = average[1]*norm;
    double M_exp = average[2]*norm;
    double M_exp2 = average[3]*norm;
    M_abs = average[4]*norm;
    Cv = (E_exp2-E_exp*E_exp)/T/T;
    chi = (M_exp2-M_abs*M_abs)/T;
    double sigmaE = (E_exp2 - E_exp*E_exp)/n_spins/n_spins;

    E_exp /= n_spins*n_spins;
    M_abs /= n_spins*n_spins;
    Cv /= n_spins*n_spins;
    chi /= n_spins*n_spins;

    // Write solution to file
    file << T << "\t" << mc << "\t" << E_exp << "\t" << M_abs << "\t" << Cv << "\t" << chi << "\t" << sigmaE << "\t" << accepted_configs*norm << endl;
}
