"""THE 2-DIMENSIONAL ISING MODEL

A program that plots:
1) the total number of accepted configurations as a function of the 
total number of Monte Carlo cycles,
2) the total number of accepted configurations as a function of temperature
for a 2-dimensional Ising model using the Metropolis algorithm."""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})


def exp_values(filename):
    """
    Function that plots various expectation values as a function of 
    Monte Carlo cycles.
    File structure: [MC_cycles,E_avg,M_absavg,C_v,X]
    """
    
    data = np.loadtxt(filename,unpack=True)     # Read data
    
    # Model
    # L = data[0]                          # Lattice dimension, L x L
    # T = data[1]                          # Temperature [kT/J]
    MC_cycles = data[0]                  # Number of Monte Carlo cycles

    # Expectation values
    E_avg = data[1]                      # Mean energy
    M_absavg = data[2]                   # Mean magnetization (absolute value) 
    C_v = data[3]                        # Specific heat
    X = data[4]                          # Susceptibility
    
    # Plotting
    plt.figure()
    #plt.plot(MC_cycles,E_avg,label=r'$<E>$')
    #plt.plot(MC_cycles,M_absavg,label=r'$<|M|>$')
    plt.plot(MC_cycles,C_v,label=r'$C_v$')
    plt.plot(MC_cycles,X,label=r'$\chi$')
    #plt.title('Ising model, L = %d, T = %d kT/J' % (L,T))
    plt.xlabel('# of MC cycles',size=14)
    plt.ylabel('Expectation values',size=14)
    plt.legend()
    plt.show()
    
    return 


def accepted_config(filename):
    """
    Function that plots the total no. of accepted configurations as a function 
    of Monte Carlo cycles and of temperature.
    File structure: [total no. of accepted configs,total no. of MC_cycles,T]
    """
    
    data = np.loadtxt(filename,True)      # Read data
    
    accepted_tot = data[0]                # Total no. of accepted configurations
    MC_cycles = data[1]                   # Total no. of Monte Carlo cycles
    T = data[2]                           # Temperature [kT/J]
    
    plt.figure()
    plt.plot(MC_cycles,accepted_tot)
    plt.xlabel('# of MC cycles')
    plt.ylabel('# of accepted configurations')
    plt.title('T = %d' % T)
    
    plt.figure()
    plt.plot(T,accepted_tot)
    plt.xlabel('Temperature [kT/J]')
    plt.ylabel('# of accepted configurations')
    
    plt.show()
    
    return

	
def main(argv):
    filename = 'C:\Users\mariefoss\Documents\GitHub\Project4\Project4\ExpectationValues_MCsT1.txt'
    exp_values(filename)
    # accepted_config(filename)

if __name__ == "__main__":
    main(sys.argv[1:]) 
