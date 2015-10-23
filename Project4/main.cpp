#include <iostream>
#include <cmath>

using namespace std;

int main()
{
    double k = 1; // Boltzmann's constant
    int L,N;
    double T=1; // units kT/J
    double beta = 1/(k*T);
    double Cv,chi,Z;
    double E_exp,M_exp,E_exp2,M_exp2,sigmaE,sigmaM; // Expectation values
    Z = 0;
    E_exp = 0;
    M_exp = 0;
    E_exp2 = 0;
    M_exp2 = 0;

    // Periodic boundary conditions and the Metropolis algorithm

    // a
    L = 2;
    N = 2*2;
    double *spin = new double[L];
    double *E = new double[N];
    double *M = new double[N];
    for(int i=0;i<N;i++){
        E[i] = 0;
        M[i] = 0;
        for(int l=0;l<L;l++)
            spin[l] = 1;
        for(int k=i;k<L;k++){
            E[i] += spin[k]*spin[k+1];
            M[i] += spin[k];
        }
        Z += exp(-beta*E[i]);
        E_exp +=  E[i]*exp(-beta*E[i]);
        E_exp2 += E[i]*E[i]*exp(-beta*E[i]);
        M_exp +=  M[i]*exp(-beta*E[i]);
        M_exp2 += M[i]*M[i]*exp(-beta*E[i]);
    }
    E_exp = E_exp/Z;
    M_exp = M_exp/Z;
    sigmaE = E_exp2/Z - E_exp*E_exp;
    sigmaM = M_exp2/Z - M_exp*M_exp;
    Cv = sigmaE*beta/T;
    chi = sigmaM*beta;
    cout << "L = " << L << endl;
    cout << "Partition function Z = " << Z << endl;
    cout << "<E> = " << E_exp << endl;
    cout << "<M> = " << M_exp << endl;
    cout << "Cv = " << Cv << endl;
    cout << "chi = " << chi << endl;

    delete[] spin;
    delete[] E;
    delete[] M;

    // b, Ising model; <E>, <|M|>, Cv, chi as functions of T
    // c, <E>, <M> as functions of MC cycles
    // d, probability as a function of E, compare with sigmaE
    // e, plotting and paralization
    // f, estimate critical temperature

    return 0;
}

