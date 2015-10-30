/* ISING MODEL IN TWO DIMENSIONS (NO MAGNTIC FIELD)
 *
 * Calculates <E>, <|M|>, Cv, chi as functions of T for a square lattice with 2 possible spin values
 * using the Metropolis algorithm.
*/

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
//#include <mpi.h>

using namespace std;

inline int periodic(int i, int limit, int add){return (i+limit+add)%(limit);}
void Metropolis(int,int **,double &,double &,double *,int&);
void initialize(int,int**,double&,double&);
void output(int,int,double,double*);
void ExpectationValues_toFile(double,ofstream&,int,int,int**,bool&);
void analytical_cf(double T,double& Z,double& E_exp,double& M_exp,double& Cv,double& chi);
void analytical_sums(double T,double& Z,double& E_exp,double& M_exp,double& Cv,double& chi, int n_spins,int N,int** spin_matrix);


int main()
{
    // CONSTANTS
    double T = 1.0;                                     // Temperature [kT/J]
    int n_spins = 20;                                   // Number of spins
    int N = n_spins*n_spins;                            // Lattice dimensions (square)

    double temp_step = 0.5;                             // Steps in temperature
    double initial_temp = 1.0;                          // Initial temperature
    double final_temp = 3.0;                            // Final temperature

    int MC_cycles = 200;                                // Number of Monte Carlo cycles


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


    // ANALYTICAL SOLUTION: SUMS
    cout << "ANALYTICAL SOLUTION: SUMS" << endl;

    double Z,E_exp,M_exp,C_V,chi;
    Z = E_exp = M_exp = C_V = chi = 0;

    analytical_sums(T,Z,E_exp,M_exp,C_V,chi,n_spins,N,spin_matrix);

    cout << "n_spins" << "\t" << n_spins << endl;
    cout << "Z" << "\t" << Z << endl;
    cout << "<E>" << "\t" << E_exp << endl;
    cout << "<M>" << "\t" << M_exp << endl;
    cout << "Cv" << "\t" << C_V << endl;
    cout << "chi" << "\t" << chi << endl << endl;



    // ANALYTICAL SOLUTION: CLOSED-FORM
    cout << "ANALYTICAL SOLUTION: CLOSED-FORM" << endl;

    double Z_cf,E_exp_cf,M_exp_cf,C_V_cf,chi_cf;
    Z_cf = E_exp_cf = M_exp_cf = C_V_cf = chi_cf = 0;

    analytical_cf(T,Z_cf,E_exp_cf,M_exp_cf,C_V_cf,chi_cf);

    cout << "n_spins" << "\t" << n_spins << endl;
    cout << "Z" << "\t" << Z_cf << endl;
    cout << "<E>" << "\t" << E_exp_cf << endl;
    cout << "<M>" << "\t" << M_exp_cf << endl;
    cout << "Cv" << "\t" << C_V_cf << endl;
    cout << "chi" << "\t" << chi_cf << endl << endl;


    // METROPOLIS ALGORITHM
    cout << "METROPOLIS ALGORITHM" << endl;

    // Expectation values as a function of MC cycles
    ofstream file_MC("ExpectationValues_MC.txt");
    bool first = true;

    for(int MC = 1; MC<=MC_cycles; MC++){
        ExpectationValues_toFile(T,file_MC,n_spins,MC,spin_matrix,first);
    }
    file_MC.close();

    // Expectation values as a function of temperature variations
    ofstream file_T("ExpectationValues_temp.txt");

    for(T=initial_temp; T<=final_temp ; T+= temp_step)
        ExpectationValues_toFile(T,file_T,n_spins,MC_cycles,spin_matrix,first);
    file_T.close();

    // Probability P(E_i), plot comparing to sigmaE
    /*
    int E_count,E_tot;
    metropolis(...,counting=true)
    if(counting){
        E_tot++;
        if(E == E_i) E_count++;
    }
    Probability = E_count/E_tot;

    //ofstream fileP("ProbabilityL" + to_string(L) + ".txt"); // File for expectation values
    //fileP << L << "\t" << Probability << "\t" << sigmaE << endl; // Write solution to file
    //fileP.close();
    */

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
    chi = 32/T*(exp(8/T) + 1)*Z_inv_cf*Z_inv_cf;
}

void analytical_sums(double T,double& Z,double& E_exp,double& M_exp,double& Cv,double& chi,int n_spins,int N,int** spin_matrix){
    // Initialize sums
    double E_exp2 = 0;                          // Expectation value squared, energy <E^2>
    double M_abs = 0;                           // Expectation value, absolute value of net magnetization <|M|>
    double M_exp2 = 0;                          // Expextation value squared, net magnetization <M^2>

    // Calculate energies
    int jm,km;

    double En[N];                               // Energy vector, length N //double *En = new double[N];
    double Ma[N];                               // Magnetization vector, length N //double *Ma = new double[N];

    for(int i=0;i<n_spins;i++){
        En[i] = 0;
        Ma[i] = 0;

        for(int k=0;k<n_spins;k++){
            // Periodic boundary conditions, 1st dimension
            if(k==n_spins-1){km=0;}
            else{km=k+1;}

            for(int j=0;j<n_spins;j++){
                // Periodic boundary conditions, 2nd dimension
                if(j==n_spins-1){jm=0;}
                else{jm=j+1;}

                // Calculate energy and net magnetization
                En[i] -= 2*spin_matrix[k][jm]*spin_matrix[km][j];
                Ma[i] += spin_matrix[k][j];
            }
        }

        Z += exp(-En[i]/T);
        E_exp +=  En[i]*exp(-En[i]/T);
        E_exp2 += En[i]*En[i]*exp(-En[i]/T);
        M_exp +=  Ma[i]*exp(-En[i]/T);
        M_exp2 += Ma[i]*Ma[i]*exp(-En[i]/T);
        M_abs += fabs(Ma[i])*exp(-En[i]/T);
    }

    // Calculate expectation values
    E_exp = E_exp/Z;
    M_exp = M_exp/Z;
    double sigmaE = E_exp2/Z - E_exp*E_exp;
    double sigmaM = M_exp2/Z - M_exp*M_exp;

    Cv = sigmaE/(T*T);                       // Heat capacity
    chi = sigmaM/T;                          // Susceptibility
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


void Metropolis(int n_spins, int **spin_matrix, double& E, double&M, double *w, int& accepted_configs)
{
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

void output(int n_spins,int MCs,double temp,double *average){
    double norm = 1./((double) (MCs));
    double Eaverage = average[0]*norm;
    double E2average = average[1]*norm;
    double Maverage = average[2]*norm;
    double M2average = average[3]*norm;
    double Mabsaverage = average[4]*norm;

    // Expectation values per spin
    double Evariance = (E2average - Eaverage*Eaverage)/n_spins/n_spins;
    double Mvariance = (M2average - Maverage*Maverage)/n_spins/n_spins;
    double M2variance = (M2average - Mabsaverage*Mabsaverage)/n_spins/n_spins;

    // Print results
    cout << "T = "             << "\t" << temp << endl;
    cout << "<E> = "           << "\t" << Eaverage << endl;
    cout << "<E>/spin = "      << "\t" << Eaverage/n_spins/n_spins << endl;
    cout << "Cv = "            << "\t" << Evariance/temp/temp << endl;
    cout << "X = "             << "\t" << M2variance/temp << endl;
    cout << "<M> = "           << "\t" << Maverage << endl;
    cout << "<M>/spin = "      << "\t" << Mabsaverage/n_spins/n_spins << endl << endl;
}

void ExpectationValues_toFile(double T,ofstream &file,int n_spins,int mc,int**spin_matrix,bool &first){
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

    for(int i=0;i<5;i++) average[i] = 0;

    initialize(n_spins,spin_matrix,E,M);

    for(int cycles=1;cycles <= mc;cycles++){
        Metropolis(n_spins,spin_matrix,E,M,w,accepted_configs);

        // Update expectation values
        double Eprev = average[0];
        average[0] += E; average[1] += E*E;
        average[2] += M; average[3] += M*M; average[4] += fabs(M);

        test = fabs((Eprev-average[0])/Eprev);
        if (test < 0.05) countstart = 1;
    }

    if (countstart == 1 && first){
        cout << "Minimum MC cycles = " << "\t" << mc << endl;
        first = false;
    }

    double norm = 1/((double) (mc));
    double E_exp = average[0]*norm;
    double E_exp2 = average[1]*norm;
    double M_exp = average[2]*norm;
    double M_exp2 = average[3]*norm;
    double M_abs = average[4]*norm;
    double Cv = (E_exp2-E_exp*E_exp)/T/T;
    double chi = (M_exp2-M_abs*M_abs)/T;

    E_exp /= n_spins*n_spins;
    M_abs /= n_spins*n_spins;
    Cv /= n_spins*n_spins;
    chi /= n_spins*n_spins;

    // Write solution to file
    file << T << "\t" << mc << "\t" << E_exp << "\t" << M_abs << "\t" << Cv << "\t" << chi << "\t" << accepted_configs*norm << endl;
}
