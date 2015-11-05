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


def read_file(filename):
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
    

def probability(filenameMC,filename_E,L):
    data_MC = np.loadtxt(filenameMC,unpack=True)   # Read data
    T = data_MC[0]                                 # Temperature [kT/J]
    MC_cycles = data_MC[1]                         # Number of Monte Carlo cycles
    E_avg = data_MC[2]                             # Expectation value of energy
    E_var = data_MC[6]                             # Variance of energy
        
    E = np.loadtxt(filename_E,unpack=True)         # Energies
    
    # Round-offs to make counting possible
    precision = 2
    E_avg = np.around(E_avg,decimals=precision)
    E = np.around(E,decimals=precision)

    # Calculate probability
    no_elements = len(E)
    counter = list(E).count(E_avg[-1])
    prob = float(counter)/no_elements
    
    print "T = %.2f" % T[0] 
    print "MC cycles = %d" % MC_cycles[-1]
    print
    
    print "Probability:"
    print "P(E) = %.2f" % prob
    print "Var(E) = %.2f" % E_var[-1]

    return


def exp_values_MC(L,temps):
    # L = array or list with dimension of lattice(s)
    # temps = array or list of temperatures
    
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
        filenameMC = 'ExpectationValues_MC_%d_%.6f.txt' % (L,temps[i])
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
        
        # Plot expectation values
        plt.figure(1)
        plt.plot(MC_cycles,E_avg,label=r'$T =$ %.2f' % temps[i])
        plt.legend()
        
        plt.figure(2)
        plt.plot(MC_cycles,M_absavg,label=r'$T =$ %.2f' % temps[i])
        plt.legend()
        
        plt.figure(3)
        plt.plot(MC_cycles,C_v,label=r'$T =$ %.2f' % temps[i])
        plt.legend()
        
        plt.figure(4)
        plt.plot(MC_cycles,X,label=r'$T =$ %.2f' % temps[i])
        plt.legend()
        
        plt.figure(5)
        plt.plot(MC_cycles,AC,label=r'$T =$ %.2f' % temps[i])
        plt.legend()
    
    plt.show()
    
    return

def exp_values_T(L,T):
    # L = lattice size
    # T = temperature
    
    # Define plots
    plt.figure(6)
    plt.title('Energy expectation value, T = %.2f' % T,size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$\langle E \rangle$',size=12)
    
    plt.figure(7)
    plt.title('Magnetism expectation value, T = %.2f' % T,size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$\langle M \rangle$',size=12)
    
    plt.figure(8)
    plt.title('Heat capacity expectation value, T = %.2f' % T,size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$C_V$',size=12)
    
    plt.figure(9)
    plt.title('Susceptibility expectation value, T = %.2f' % T,size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$',size=12)
    plt.ylabel(r'$\chi$',size=12)
    
    plt.figure(10)
    plt.title('Accepted configurations vs. temperature, T = %.2f' % T,size=12)
    plt.xlabel(r'$T \mathrm{[kT/J]}$')
    plt.ylabel('# of accepted configurations',size=12)

    for i in range(len(L)):
        
        # Read file
        filenameT = 'ExpectationValues_temp_%d_%.6f.txt' % (L[i],T)
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameT)
        
        # Plot
        plt.figure(6)
        plt.plot(T,E_avg,label='L = %d' % L[i])
        plt.legend()
        
        plt.figure(7)
        plt.plot(T,M_absavg,label='L = %d' % L[i])
        plt.legend()
        
        plt.figure(8)
        plt.plot(T,C_v,label='L = %d' % L[i])
        plt.legend()
        
        plt.figure(9)
        plt.plot(T,X,label='L = %d' % L[i])
        plt.legend()
        
        plt.figure(10)
        plt.plot(T,AC,label='L = %d' % L[i])
        plt.legend()
    
    plt.show()
    
    return

def main(argv):

    # Change depending on situation!!
    L = 20
    temps = [1.0, 2.4]
    
    #exp_values_MC(L,temps)
    
    L = [2,20]
    exp_values_T(L,temps[0])
    
        

if __name__ == "__main__":
    main(sys.argv[1:]) 