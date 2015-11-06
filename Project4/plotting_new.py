"""THE 2-DIMENSIONAL ISING MODEL

A program that plots:
    1) Expectation values,
    2) Number of accepted configurations, and
    3) Number of Monte Carlo cycles
as a function of temperature and number of Monte Carlo cycles.

The program also calculates the probability of the average energy
and compares it to the standard deviation for the average energy."""


import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})


def read_file(filename):
    # Input: filename
    # Output: Arrays of data values
    
    data = np.loadtxt(filename,unpack=True) # Read data
    
    # Model   
    T = data[0]                     # Temperature [kT/J]
    MC_cycles = data[1]             # Number of Monte Carlo cycles

    # Expectation values
    E_avg = data[2]                 # Mean energy
    M_absavg = data[3]              # Mean magnetization (absolute value) 
    C_v = data[4]                   # Specific heat
    X = data[5] 		    # Susceptibility
    E_var = data[6]                 # Variance of energy
    AC = data[7]		    # Accepted configurations
    
    return T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC



def exp_values_MC(L,temps):
    # Input: L = array or list with dimension of lattice(s), temps = array or list of temperatures
    # Output: Plots of expectation values vs. # of Monte Carlo cycles
    
    # Define plots
    plt.figure(1)
    plt.title('Energy expectation value, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle E \rangle$',size=12)
    
    plt.figure(2)
    plt.title('Magnetism expectation value, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle M \rangle$',size=12)
    
    plt.figure(3)
    plt.title('Heat capacity, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$C_V$',size=12)
    
    plt.figure(4)
    plt.title('Susceptibility, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\chi$',size=12)
    
    plt.figure(5)
    plt.title('Accepted configurations vs. MC cycles, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel('Accepted configurations',size=12)    
    
    for i in range(len(temps)):
        # Read file
        filenameMC = '../build-Project4-Desktop_Qt_5_5_0_clang_64bit-Debug/ExpectationValues_MC_%d_%.6f.txt' % (L,temps[i])
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
        
        # Fill plots with values
        plt.figure(1)
        plt.plot(MC_cycles,E_avg,label=r'$T =$ %.2f' % temps[i])
        plt.ylim(-2.2,-1.0)
        plt.legend()
        
        plt.figure(2)
        plt.plot(MC_cycles,M_absavg,label=r'$T =$ %.2f' % temps[i])
        plt.ylim(0.7,1.1)
        plt.legend()
        
        plt.figure(3)
        plt.plot(MC_cycles,C_v,label=r'$T =$ %.2f' % temps[i])
        plt.ylim(-0.5,1.0)
        plt.legend()
        
        plt.figure(4)
        plt.plot(MC_cycles,X,label=r'$T =$ %.2f' % temps[i])
        plt.ylim(-0.1,0.3)
        plt.legend()
        
        plt.figure(5)
        plt.plot(MC_cycles,AC,label=r'$T =$ %.2f' % temps[i])
        plt.ylim(-0.5,2.0)
        plt.legend()
    
    plt.show()
    
    return


def exp_values_T(L):
    # Input: L = array or list of lattice sizes
    # Output: Plots expectation values vs. temperature
    
    # Define plots
    plt.figure(6)
    plt.title('Energy',size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$\langle E \rangle$',size=12)
    
    plt.figure(7)
    plt.title('Magnetism',size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$\langle M \rangle$',size=12)
    
    plt.figure(8)
    plt.title('Heat capacity',size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$C_V$',size=12)
    
    plt.figure(9)
    plt.title('Susceptibility',size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$\chi$',size=12)
    
    plt.figure(10)
    plt.title('Accepted configurations vs. temperature',size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$')
    plt.ylabel('# of accepted configurations',size=12)

    for i in range(len(L)):
        # Read file
        filenameT = '../build-Project4-Desktop_Qt_5_5_0_clang_64bit-Debug/ExpectationValues_temp_%d.txt' % (L[i])
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameT)
        
        # Fill plots with values
        plt.figure(6)
        plt.plot(T,E_avg,label='L = %d' % L[i])
        plt.ylim(-2.5,0)
        plt.legend()
        
        plt.figure(7)
        plt.plot(T,M_absavg,label='L = %d' % L[i])
        plt.ylim(0,1.4)
        plt.legend()
        
        plt.figure(8)
        plt.plot(T,C_v,label='L = %d' % L[i])
        plt.ylim(0,4)
        plt.legend()
        
        plt.figure(9)
        plt.plot(T,X,label='L = %d' % L[i])
        plt.ylim(0,250)
        plt.legend()
        
        plt.figure(10)
        plt.plot(T,AC,label='L = %d' % L[i])
        plt.legend()
    
    plt.show()
    
    return
    
    
def probability(L,temp):
    # Input: L = lattice size, temp = temperature
    # Output: Probability P(E) of a given E and variance of E
    
    filenameMC = '../build-Project4-Desktop_Qt_5_5_0_clang_64bit-Debug/ExpectationValues_MC_%d_%.6f.txt' % (L,temp)
    filenameE = '../build-Project4-Desktop_Qt_5_5_0_clang_64bit-Debug/Energy_MC_%d_%.6f.txt' % (L,temp)
    
    T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
        
    E = np.loadtxt(filenameE,unpack=True)         # Energies
    
    # Round-offs to make counting possible
    precision = 2 # Number of decimal points
    E_avg = np.around(E_avg,decimals=precision)
    E = np.around(E,decimals=precision)

    # Calculate probability
    no_elements = len(E)
    counter = list(E).count(E_avg[-1])
    prob = float(counter)/no_elements
    
    print "Probability: (T = %.2f)" % temp
    print "P(E) = %.2f" % prob
    print "Var(E) = %.2f" % E_var[-1]

    return
    

def main(argv):
    
    # Expectation values vs. # of Monte Carlo cycles
    L = 2
    temps = [1.0, 2.4]
    exp_values_MC(L,temps)
    
    # Expectation values vs. temperature
    L_list = [2,20,40,80]
    #exp_values_T(L_list)
    
    # Calculate probabilty for <E>
    #probability(L,temps[0])

if __name__ == "__main__":
    main(sys.argv[1:]) 