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


def exp_valuesT(filename,L):
    """
    Function that plots various expectation values as a function of temperature.
    File structure: [Temperature,MC_cycles,E_avg,M_absavg,C_v,X, number of accepted]
    """
    
    data = np.loadtxt(filename,unpack=True)     # Read data
    
    # Model   
    T = data[0]                     # Temperature [kT/J]

    # Expectation values
    E_avg = data[2]                 # Mean energy
    M_absavg = data[3]              # Mean magnetization (absolute value) 
    C_v = data[4]                   # Specific heat
    X = data[5] 		    # Susceptibility
    E_var = data[6]                 # Variance of energy
    AC = data[7]		    # Accepted configurations
    
    # Plotting
    plt.figure()
    plt.title('Expectation values vs. temperature, L = %d' % L,size=12)
    plt.plot(T,E_avg,'b',label=r'$<E>$')
    plt.plot(T,M_absavg,'r',label=r'$<|M|>$')
    plt.ylim(-2.2,1.2)
    plt.xlabel(r'Temperature ($kT/J$)')
    plt.ylabel('Expectation values',size=12)
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title('Expectation values vs. temperature, L = %d' % L,size=12)
    plt.plot(T,C_v,'b',label=r'$C_v$')
    plt.plot(T,X,'r',label=r'$\chi$')
    plt.xlabel(r'Temperature ($kT/J$)')
    plt.ylabel('Expectation values',size=12)
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title('Accepted configurations vs. temperature, L = %d' % L,size=12)
    plt.plot(T,AC,label=r'AC')
    plt.xlabel(r'Temperature ($kT/J$)')
    plt.ylabel('# of accepted configurations',size=12)
    plt.show()
    
    return 
	
def exp_valuesMC(filename,L):
    """
    Function that plots various expectation values as a function of 
    Monte Carlo cycles.
    File structure: [Temperature,MC_cycles,E_avg,M_absavg,C_v,X, number of accepted]
    """
    
    data = np.loadtxt(filename,unpack=True)     # Read data
    
    # Model
    T = data[0]                     # Temperature [kT/J]
    MC_cycles = data[1]             # Number of Monte Carlo cycles

    # Expectation values
    E_avg = data[2]                 # Mean energy
    M_absavg = data[3]              # Mean magnetization (absolute value) 
    C_v = data[4]                   # Specific heat
    X = data[5] 					# Susceptibility
    E_var = data[6]                 # Variance of energy
    AC = data[7]					# Accepted configurations
    
    # Plotting
    plt.figure()
    plt.title('Expectation values vs. MC cycles, L = %d, T = %.1f' % (L,T[0]),size=12)
    plt.plot(MC_cycles,E_avg,'b',label=r'$<E>$')
    plt.plot(MC_cycles,M_absavg,'r',label=r'$<|M|>$')
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel('Expectation values',size=12)
    plt.ylim(min(E_avg)-1,max(M_absavg)+1)
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title('Expectation values vs. MC cycles, L = %d, T = %.1f' % (L,T[0]),size=12)
    plt.plot(MC_cycles,C_v,'b',label=r'$Cv$')
    plt.plot(MC_cycles,X,'r',label=r'$\chi$')
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel('Expectation values',size=12)
    plt.legend()
    plt.show()
    
    plt.figure()
    plt.title('Accepted configurations vs. MC cycles, L = %d, T = %.1f' % (L,T[0]),size=12)
    plt.plot(MC_cycles,AC)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel('Accepted configurations',size=12)
    plt.show()
    
    return 
	
def exp_valuesTL(filename2,filename4,filename6,filename8): # Is there a better way to write this??
    """
    Function that plots various expectation values as a function of temperature
    for a given lattice size LxL.
    File structure: [Temperature,E_avg,M_absavg,C_v,X]
    """
    
    data2 = np.loadtxt(filename2,unpack=True)
    data4 = np.loadtxt(filename4,unpack=True)
    data6 = np.loadtxt(filename6,unpack=True)
    data8 = np.loadtxt(filename8,unpack=True)
    
    # Model
    T2 = data2[1]
    T4 = data4[1]
    T6 = data6[1]
    T8 = data8[1]

    # Expectation values
    E_avg2 = data2[2]
    M_absavg2 = data2[3] 
    C_v2 = data2[4]
    X2 = data2[5]
	
    E_avg4 = data4[2]
    M_absavg4 = data4[3] 
    C_v4 = data4[4]
    X4 = data4[5]
	
    E_avg6 = data6[2]
    M_absavg6 = data6[3] 
    C_v6 = data6[4]
    X6 = data6[5]
	
    E_avg8 = data8[2]
    M_absavg8 = data8[3] 
    C_v8 = data8[4]
    X8 = data8[5]
    
    # Plotting, subplot possible
    plt.plot(T2,E_avg2,label=r'$20\times20$')
    plt.plot(T4,E_avg4,label=r'$40\times40$')
    plt.plot(T6,E_avg6,label=r'$60\times60$')
    plt.plot(T8,E_avg8,label=r'$80\times80$')
    plt.xlabel(r'Temperature ($kT/J$)')
    plt.ylabel('Mean energy',size=14)
    plt.legend()
    plt.show()
	
    plt.plot(T2,M_absavg2,label=r'$20\times20$')
    plt.plot(T4,M_absavg4,label=r'$40\times40$')
    plt.plot(T6,M_absavg6,label=r'$60\times60$')
    plt.plot(T8,M_absavg8,label=r'$80\times80$')
    plt.xlabel(r'Temperature ($kT/J$)')
    plt.ylabel('Mean magnetization',size=14)
    plt.legend()
    plt.show()
	
    plt.plot(T2,C_v2,label=r'$20\times20$')
    plt.plot(T4,C_v4,label=r'$40\times40$')
    plt.plot(T6,C_v6,label=r'$60\times60$')
    plt.plot(T8,C_v8,label=r'$80\times80$')
    plt.xlabel(r'Temperature ($kT/J$)')
    plt.ylabel('Specific heat',size=14)
    plt.legend()
    plt.show()
	
    plt.plot(T2,X2,label=r'$20\times20$')
    plt.plot(T4,X4,label=r'$40\times40$')
    plt.plot(T6,X6,label=r'$60\times60$')
    plt.plot(T8,X8,label=r'$80\times80$')
    plt.xlabel(r'Temperature ($kT/J$)')
    plt.ylabel('Susceptibility',size=14)
    plt.legend()
    plt.show()
    
    return 
    

def probability(filename,L):
    data = np.loadtxt(filename,unpack=True)     # Read data
    T = data[0]                                 # Temperature [kT/J]
    MC_cycles = data[1]                         # Number of Monte Carlo cycles
    E_avg = data[2]
    E_var = data[6]                 # Variance of energy

    
    plt.figure()
    plt.title('Probability P(E), L = %d, T = %.1f' % (L,T[0]),size=12)
    plt.hist(E_avg)
    plt.xlabel('Energy values')
    plt.ylabel('Frequency')
    plt.show()
    
    #plt.figure()
    #plt.title('Variance of the energy, L = %d, T = %.1f' % (L,T[0]),size=12)
    #plt.plot(MC_cycles,E_var)
    #plt.xlabel('# of MC cycles',size=12)
    #plt.ylabel(r'$\sigma_E^2$')
    #plt.show
    return
	
def main(argv):
    
    # Expectation values as functions of Monte Carlo cycles
    L = 20                                             # Lattice size
    T = 1.0                                            # Temperature [kT/J]
    
    filenameMC = 'ExpectationValues_MC_%d_%.6f.txt' % (L,T)
    # filenameMC1 = 'C:\Users\mariefoss\Documents\GitHub\Project4\Project4\ExpectationValues_MCsT1.txt'
    #exp_valuesMC(filenameMC,L)
    
    probability(filenameMC,L)
	
    # Expectation values as functions of temperature, a given L
    filenameT = 'ExpectationValues_temp_%d_%.6f.txt' % (L,T)
    #filenameT = 'C:\Users\mariefoss\Documents\GitHub\Project4\Project4\ExpectationValues_temp.txt'
    #exp_valuesT(filenameT,L)
	
	
    # Expectation values as functions of temperature
    L = [20,40,60,80]                                 # Range of lattice sizes
    T = 1.0                                           # Temperature [kT/J]
    filenameTL = np.zeros(len(L))
    
    for i in range(len(L)):
        filenameTL[i] = 'ExpectationValues_temp_%d_%.6f.txt' % (L[i],T)
        
    # exp_valuesT(filenameTL[0],filenameTL[1],filenameTL[2],filenameTL[3])
	

if __name__ == "__main__":
    main(sys.argv[1:]) 
