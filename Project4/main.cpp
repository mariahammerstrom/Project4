#include <iostream>
#include <cmath>
//#include "lib.h"
#include <random>
#include <fstream>

using namespace std;

inline int periodic(int i, int limit, int add){return (i+limit+add)%(limit);}
void Metropolis(int,int **,double &,double &,double *);
void initialize(int,double,int**,double&,double&);
void output(int,int,double,double*);

int main()
{
    double k = 1; // Boltzmann's constant
    int L,N;
    double T=1; // units kT/J
    double Cv,chi,Z;
    double E_exp,M_exp,E_exp2,M_exp2,sigmaE,sigmaM; // Expectation values

    // Periodic boundary conditions and the Metropolis algorithm

    // Analytical solution
   /* L = 2;
    N = 2*2;
    int jm;
    Z = 0;
    E_exp = 0;
    M_exp = 0;
    E_exp2 = 0;
    M_exp2 = 0;
    int *spin = new int[L];
    double *En = new double[N];
    double *Ma = new double[N];
    for(int i=0;i<N;i++){
        En[i] = 0;
        Ma[i] = 0;
        for(int l=0;l<N;l++){
            spin[l] = 1;
        }
        for(int k=0;k<N;k++){
            if(k==N-1){jm=0;}
            else{jm=k+1;}
            En[i] -= spin[k]*spin[jm];
            Ma[i] += spin[k];
        }
        Z += exp(-En[i]/T);
        E_exp +=  En[i]*exp(-En[i]/T);
        E_exp2 += En[i]*En[i]*exp(-En[i]/T);
        M_exp +=  Ma[i]*exp(-En[i]/T);
        M_exp2 += Ma[i]*Ma[i]*exp(-En[i]/T);
    }
    E_exp = E_exp/Z;
    M_exp = M_exp/Z;
    sigmaE = E_exp2/Z - E_exp*E_exp;
    sigmaM = M_exp2/Z - M_exp*M_exp;
    Cv = sigmaE/(T*T);
    chi = sigmaM/T;
    cout << "L = " << L << endl;
    cout << "Partition function Z = " << Z << endl;
    cout << "<E> = " << E_exp << endl;
    cout << "<M> = " << M_exp << endl;
    cout << "Cv = " << Cv << endl;
    cout << "chi = " << chi << endl;

    delete[] spin;
    delete[] En;
    delete[] Ma;*/

    // b, Ising model; <E>, <|M|>, Cv, chi as functions of T
    int **spin_matrix,n_spins,MCs;
    double w[17],average[5],initial_temp,final_temp,temp_step;
    temp_step = 0.5;
    initial_temp = 1.0;
    final_temp = 3.0;
    MCs = 5;
    n_spins = 2;
    spin_matrix = new int*[n_spins];
    for(int i=0;i<n_spins;i++)
        spin_matrix[i] = new int[n_spins];
    //spin_matrix = (int**)matrix(n_spins,n_spins,sizeof(int));
    for(int i=0;i<n_spins;i++){
        for(int j=0;j<n_spins;j++){
            spin_matrix[i][j] = 1;
            cout << spin_matrix[i][j] << endl;
        }
    }

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
            Metropolis(n_spins,spin_matrix,E,M,w);

            // Update expectation values
            average[0] += E; average[1] += E*E;
            average[2] += M; average[3] += M*M; average[4] += fabs(M);
        }
        // Print results
        output(n_spins,MCs,temp,average);


    }

/*
    // Expectation Values as functions of MC cycles
    L = 20;
    T= 1;
    ofstream fileT1("ExpectationValues_MCsT" + to_string(T) + ".txt"); // File for expectation values
    fileT1 << MCs << "\t" << E_exp << "\t" << M_exp << "\t" << Cv << "\t" << chi << endl; // Write solution to file
    fileT1.close();

    T = 2.4;
    ofstream fileT24("ExpectationValues_MCsT" + to_string(T) + ".txt"); // File for expectation values
    fileT24 << MCs << "\t" << E_exp << "\t" << M_exp << "\t" << Cv << "\t" << chi << endl; // Write solution to file
    fileT24.close();

    // Accepted conficurations
    ofstream fileACMC("AcceptedConfigurationsMC.txt"); // File for expectation values
    fileACMC << MCs << "\t" << AC << endl; // Write solution to file
    fileACMC.close();
    ofstream fileACT("AcceptedConfigurationsT.txt"); // File for expectation values
    fileACT << T << "\t" << AC << endl; // Write solution to file
    fileACT.close();


    // Probability
    ofstream fileP("ProbabilityL" + to_string(L) + ".txt"); // File for expectation values
    fileP << L << Probability << "\t" << sigmaE << endl; // Write solution to file
    fileP.close();


    // Ecpectation values for different L and T. Parallelization!
    ofstream fileL("ExpectationValuesL" + to_string(L) + ".txt"); // File for expectation values
    fileL << T << E_exp << "\t" << M_exp << "\t" << Cv << "\t" << chi << endl;
    fileL.close();
    */

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


void Metropolis(int n_spins, int **spin_matrix, double& E, double&M, double *w)
{
    default_random_engine generator;                   // start random number generator
    uniform_real_distribution<double> distribution(0.0,1.0);   // pick type of random distribution
    // loop over all spins
    for(int y =0; y < n_spins; y++) {
        for (int x= 0; x < n_spins; x++){
            // Find random position
            int ix = (int) (distribution(generator)*(double)n_spins);
            int iy = (int) (distribution(generator)*(double)n_spins);
            int deltaE = 2*spin_matrix[iy][ix]*(spin_matrix[iy][periodic(ix,n_spins,-1)]+ spin_matrix[periodic(iy,n_spins,-1)][ix] + spin_matrix[iy][periodic(ix,n_spins,1)] + spin_matrix[periodic(iy,n_spins,1)][ix]);
            // Here we perform the Metropolis test
            if ( distribution(generator) <= w[deltaE+8] ) {
                spin_matrix[iy][ix] *= -1; // flip one spin and accept new spin config
                // update energy and magnetization
                M += (double) 2*spin_matrix[iy][ix];
                E += (double) deltaE;
            }
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

    // Print
    cout << "T = "             << "\t" << temp << endl;
    cout << "<E> = "              << "\t" << Eaverage << endl;
    cout << "<E>/spin = "      << "\t" << Eaverage/n_spins/n_spins << endl;
    cout << "Cv = "            << "\t" << Evariance/temp/temp << endl;
    cout << "X = "             << "\t" << M2variance/temp << endl;
    cout << "<M> = "              << "\t" << Maverage << endl;
    cout << "<M>/spin = "       << "\t" << Mabsaverage/n_spins/n_spins << endl;
}
