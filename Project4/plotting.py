# PROJECT 4: THE 2-DIMENSIONAL ISING MODEL

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})


def exp_values(filename):
    """
    Function that plots various expectation values as a function of 
    Monte Carlo cycles.
    File structure: [L,T,MC_cycles,E_avg,M_absavg,C_v,X]
    """
    
    data = np.loadtxt(filename,True)     # Read data
    
    # Model
    L = data[0]                          # Lattice dimension, L x L
    T = data[1]                          # Temperature [kT/J]
    MC_cycles = data[2]                  # Number of Monte Carlo cycles

    # Expectation values
    E_avg = data[3]                      # Mean energy
    M_absavg = data[4]                   # Mean magnetization (absolute value) 
    C_v = data[5]                        # Specific heat
    X = data[6]                          # Susceptibility
    
    # Plotting
    plt.figure()
    plt.plot(MC_cycles,E_avg,label=r'$<E>$')
    plt.plot(MC_cycles,M_absavg,label=r'$<|M|>$')
    plt.plot(MC_cycles,C_v,label=r'$C_v$')
    plt.plot(MC_cycles,X,label=r'$\chi$')
    plt.title('Ising model, L = %d, T = %d kT/J' % (L,T))
    plt.xlabel('# MC cycles')
    plt.ylabel('Expectation values')
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
    plt.xlabel('# MC cycles')
    plt.ylabel('# accepted configurations')
    
    plt.figure()
    plt.plot(T,accepted_tot)
    plt.xlabel('Temperature [kT/J]')
    plt.ylabel('# accepted configurations')
    
    plt.show()
    
    return


# Call to functions
filename = ''
exp_values(filename)
accepted_config(filename)
