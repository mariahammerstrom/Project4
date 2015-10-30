/* ISING MODEL IN TWO DIMENSIONS (NO MAGNTIC FIELD)
 *
 * Calculates <E>, <|M|>, Cv, chi as functions of T for a square lattice with 2 possible spin values
 * using the Metropolis algorithm.
*/

#include <iostream>
#include <cmath>
#include <random>
#include <fstream>

using namespace std;

inline int periodic(int i, int limit, int add){return (i+limit+add)%(limit);}
void Metropolis(int,int **,double &,double &,double *,int&);
void initialize(int,double,int**,double&,double&);
void output(int,int,double,double*);


int main()
{
    // CONSTANTS
    double T = 1;                               // Temperature [kT/J]
    double k = 1;                               // Boltzmann's constant
    double J = 1;                               // Coupling constant
    double beta = 1./(k*T);

    int L = 2;                                  // Number of spins
    int N = 2*2;                                // Lattice dimensions (square)
/*

    // ANALYTICAL SOLUTION: SUMS
    // Initialize sums
    double Z = 0;                               // Partition function
    double E_exp = 0;                           // Expectation value, energy <E>
    double E_exp2 = 0;                          // Expectation value squared, energy <E^2>
    double M_exp = 0;                           // Expectation value, net magnetization <M>
    double M_abs = 0;                           // Expectation value, absolute value of net magnetization <|M|>
    double M_exp2 = 0;                          // Expextation value squared, net magnetization <M^2>

    // Set up spin matrix
    int **spin; //int *spin = new int[L];
    spin = new int*[L];
    for(int i=0;i<L;i++)
        spin[i] = new int[L];
    for(int i=0;i<L;i++){
        for(int j=0;j<L;j++){
            spin[i][j] = 1;                     // Spin up: +1, spin down: -1
        }
    }

    // Calculate energies
    int jm,km;

    double En[N];                               // Energy vector, length N //double *En = new double[N];
    double Ma[N];                               // Magnetization vector, length N //double *Ma = new double[N];

    for(int i=0;i<L;i++){
        En[i] = 0;
        Ma[i] = 0;

        for(int k=0;k<L;k++){
            // Periodic boundary conditions, 1st dimension
            if(k==L-1){km=0;}
            else{km=k+1;}

            for(int j=0;j<L;j++){
                // Periodic boundary conditions, 2nd dimension
                if(j==L-1){jm=0;}
                else{jm=j+1;}

                // Calculate energy and net magnetization
                En[i] -= 2*spin[k][jm]*spin[km][j];//+spin[km][j]*spin[k][jm];
                Ma[i] += spin[k][j];
            }
        }

        Z += exp(-En[i]*beta);
        E_exp +=  En[i]*exp(-En[i]*beta);
        E_exp2 += En[i]*En[i]*exp(-En[i]*beta);
        M_exp +=  Ma[i]*exp(-En[i]*beta);
        M_exp2 += Ma[i]*Ma[i]*exp(-En[i]*beta);
        M_abs += fabs(Ma[i])*exp(-En[i]*beta);
    }

    // Calculate expectation values
    E_exp = E_exp/Z;
    M_exp = M_exp/Z;
    double sigmaE = E_exp2/Z - E_exp*E_exp;
    double sigmaM = M_exp2/Z - M_exp*M_exp;

    double Cv = sigmaE/(T*T);                       // Heat capacity
    double chi = sigmaM/T;                          // Susceptibility

    cout << "ANALYTICAL SOLUTION: SUMS" << endl;
    cout << "L" << "\t" << L << endl;
    cout << "Z" << "\t" << Z << endl;
    cout << "<E>" << "\t" << E_exp << endl;
    cout << "<M>" << "\t" << M_exp << endl;
    cout << "Cv" << "\t" << Cv << endl;
    cout << "chi" << "\t" << chi << endl << endl;


    // ANALYTICAL SOLUTION: CLOSED-FORM
    double Z_cf = 4*cosh(8*J*beta) + 12;
    double Z_inv_cf = 1./Z_cf;
    double E_exp_cf = -8*J*tanh(8*J*beta);
    double C_V_cf = k*(8*J*beta/cosh(8*J*beta))*(8*J*beta/cosh(8*J*beta));
    double M_exp_cf = 8*(exp(8*J*beta) + 2)*Z_inv_cf;
    double chi_cf = 32*beta*(exp(8*J*beta) + 1)*Z_inv_cf*Z_inv_cf;

    cout << "ANALYTICAL SOLUTION: CLOSED-FORM" << endl;
    cout << "L" << "\t" << L << endl;
    cout << "Z" << "\t" << Z_cf << endl;
    cout << "<E>" << "\t" << E_exp_cf << endl;
    cout << "<M>" << "\t" << M_exp_cf << endl;
    cout << "Cv" << "\t" << C_V_cf << endl;
    cout << "chi" << "\t" << chi_cf << endl << endl;


    // CLEAR MEMORY
    delete[] spin;
    //delete[] En;
    //delete[] Ma;
*/
/*
    // METROPOLIS ALGORITHM
    cout << "METROPOLIS ALGORITHM" << endl;

    double w[17],average[5];

    double temp_step = 0.5;                             // Steps in temperature
    double initial_temp = 1.0;                          // Initial temperature
    double final_temp = 3.0;                            // Final temperature
    int n_spins = 2;                                    // Number of spins

    int MCs = 5;                                        // No. of Monte Carlo cycles
    int accepted_configs = 0;                           // Initialize count of accepted configurations

    // Setup spin matrix
    int **spin_matrix;
    spin_matrix = new int*[n_spins];

    for(int i=0;i<n_spins;i++)
        spin_matrix[i] = new int[n_spins];
    //spin_matrix = (int**)matrix(n_spins,n_spins,sizeof(int));

    for(int i=0;i<n_spins;i++){
        for(int j=0;j<n_spins;j++){
            spin_matrix[i][j] = 1;
        }
    }


    // Run algortihm
    for(double temp = initial_temp; temp <= final_temp; temp+=temp_step){
        double E = 0;
        double M = 0;

        // Set up array for possible energy changes
        for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
        for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/temp);

        // Initialize array for expectation values
        for(int i=0;i<5;i++) average[i] = 0;
        initialize(n_spins,temp,spin_matrix,E,M);

        // Start Monte Carlo computation
        for(int cycles=1;cycles <= MCs;cycles++){
            Metropolis(n_spins,spin_matrix,E,M,w,accepted_configs);

            // Update expectation values
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
        }

        // Print results
        //output(n_spins,MCs,temp,average);
    }*/


    // Expectation Values as functions of MC cycles
    ofstream fileT1("ExpectationValues_MCsT1.txt"); // Expectation values, T = 1
    ofstream fileT24("ExpectationValues_MCsT24.txt"); // Expectation values, T = 2.4
    L = 20;
    int n_spins = 20;
    int **spin_matrix;
    spin_matrix = new int*[n_spins];
    for(int i=0;i<n_spins;i++)
        spin_matrix[i] = new int[n_spins];
    //spin_matrix = (int**)matrix(n_spins,n_spins,sizeof(int));
    for(int i=0;i<n_spins;i++){
        for(int j=0;j<n_spins;j++){
            spin_matrix[i][j] = 1;
        }
    }
    double w[17],average[5];
    double E_exp,E_exp2,M_exp,M_exp2,M_abs,Cv,chi;
    T= 1.0;
    int mc=10000; // Number of elements MC cycles array
    int MCsa[mc];
    double E = 0; double M = 0;
    int AC;
    // Set up array for possible energy changes
    for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
    for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T);
    // Initialize array for expectation values
    for(int i=0;i<5;i++) average[i] = 0;
    initialize(n_spins,T,spin_matrix,E,M);
    /*for(int i=1;i<=mc;i++){
        E = M = 0;
        AC = 0;
        MCsa[i-1] = i;
        for(int cycles=1;cycles <= MCsa[i-1];cycles++){
            Metropolis(n_spins,spin_matrix,E,M,w,AC);

            // Update expectation values
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
        }
        E_exp = average[0]/n_spins/n_spins;
        E_exp2 = average[1]/n_spins/n_spins;
        M_exp = average[2]/n_spins/n_spins;
        M_exp2 = average[3]/n_spins/n_spins;
        M_abs = average[4]/n_spins/n_spins;
        Cv = (E_exp2-E_exp*E_exp)/T/T;
        chi = (M_exp2-M_abs*M_abs)/T;
        fileT1 << MCsa[i-1] << "\t" << E_exp/((double) MCsa[i-1]) << "\t" << M_exp/((double) MCsa[i-1]) << "\t" << Cv/((double) MCsa[i-1]) << "\t" << chi/((double) MCsa[i-1]) << "\t" << AC/((double) MCsa[i-1]) << endl; // Write solution to file
    }*/

    T = 2.4; //M = E =0;
    int countstart = mc;
    double test;
    // Set up array for possible energy changes
    for(int dE = -8; dE <= 8; dE++) w[dE+8] = 0;
    for(int dE = -8; dE <= 8; dE+=4) w[dE+8] = exp(-dE/T);
    // Initialize array for expectation values
    for(int i=0;i<5;i++) average[i] = 0;
    initialize(n_spins,T,spin_matrix,E,M);
    for(int i=1;i<=mc;i+=100){
        M = E =0;
        AC = 0;
        MCsa[i-1] = i;
        for(int cycles=1;cycles <= MCsa[i-1];cycles++){
            Metropolis(n_spins,spin_matrix,E,M,w,AC);

            // Update expectation values
            double Eprev = average[0];
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
            test = fabs(Eprev-average[0]);
            if (test < 0.05 && i < countstart) countstart = i;
        }
        E_exp = average[0]/n_spins/n_spins;        
        E_exp2 = average[1]/n_spins/n_spins;
        M_exp = average[2]/n_spins/n_spins;
        M_exp2 = average[3]/n_spins/n_spins;
        double M_abs = average[4]/n_spins/n_spins;
        Cv = (E_exp2-E_exp*E_exp)/T/T;
        chi = (M_exp2-M_abs*M_abs)/T;
        fileT24 << MCsa[i-1] << "\t" << E_exp/((double) MCsa[i-1]) << "\t" << M_exp/((double) MCsa[i-1]) << "\t" << Cv/((double) MCsa[i-1]) << "\t" << chi/((double) MCsa[i-1]) << "\t" << AC/((double) MCsa[i-1]) << endl; // Write solution to file
    }
    cout << "Start to count at MC cycles = " << countstart << endl;
    fileT24.close();
    fileT1.close();

    // Accepted conficurations
    //ofstream fileACT("AcceptedConfigurationsT.txt"); // File for expectation values
    //fileACT << T << "\t" << AC << endl; // Write solution to file
    //fileACT.close();


    // Probability
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

void initialize(int n_spins,double temp,int** spin_matrix,double& E,double&M){
    for(int y=0;y<n_spins;y++){
        for(int x=0;x<n_spins;x++){
            if (temp<1.5) spin_matrix[y][x] = 1;
            M += (double) spin_matrix[y][x];
            E -= (double) spin_matrix[y][x]*(spin_matrix[periodic(y,n_spins,-1)][x] + spin_matrix[y][periodic(x,n_spins,-1)]);
        }
    }
}


void Metropolis(int n_spins, int **spin_matrix, double& E, double&M, double *w, int& AC)
{
    default_random_engine generator;                            // start random number generator
    uniform_real_distribution<double> distribution(0.0,1.0);    // pick type of random distribution

    // Loop over all spins
    for(int y=0; y<n_spins*n_spins; y++) {
        //for (int x=0; x<n_spins; x++){

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
            AC += 1;
        }
    }
    //}
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
