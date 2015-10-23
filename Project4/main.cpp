#include <iostream>

using namespace std;

int main()
{
    double k = 1; // Boltzmann's constant
    int L,N;
    double T; // units kT/J
    double J,spin,beta;
    double E,M,Cv,chi,Z;
    double E_exp,M_exp,E_exp2,M_exp2,sigmaE,sigmaM; // Expectation values
    Z = 0;
    E_exp = 0;
    M_exp = 0;
    E_exp2 = 0;
    M_exp2 = 0;
    beta = 1/(k*T);

    // Periodic boundary conditions and the Metropolis algorithm

    // a
    L = 2;
    N = 2*2;
    for(int i=0;i<N;i++){
        for(int k=i-1;k<=i+1;k++){
            for(int j=i-1;j<=i+1;j++)
                E[i] -= J*spin[k]*spin[l];
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

    // b, Ising model; <E>, <|M|>, Cv, chi as functions of T
    // c, <E>, <M> as functions of MC cycles
    // d, probability as a function of E, compare with sigmaE
    // e, plotting and paralization
    // f, estimate critical temperature

    return 0;
}

