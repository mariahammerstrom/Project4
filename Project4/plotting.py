"""

THE 2-DIMENSIONAL ISING MODEL

A program that plots:
    1) Expectation values,
    2) Number of accepted configurations, and
    3) Number of Monte Carlo cycles
as a function of temperature and number of Monte Carlo cycles.

The program also calculates:
    1) The probability of the average energy and compares it to the standard 
    deviation for the average energy, and
    2) An estimate for the critical temperature of the system.

"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc
rc('font',**{'family':'serif'})


def read_file(filename):
    # Input: filename for file with structure [T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC]
    # Output: arrays of data values
    
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


    
def exp_values_MC(L,temps,random):
    # Input: L = array or list with dimension of lattice(s), temps = array or list of temperatures
    # Output: plots of expectation values vs. # of Monte Carlo cycles
    
    # Define plots
    plt.figure(1)
    plt.title('Energy expectation value, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle E \rangle$',size=12)
    #plt.xlim(0,20000) # L = 20
    #plt.ylim(-2.1,-1.1) # L = 20
    
    plt.figure(2)
    plt.title('Magnetism expectation value, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle M \rangle$',size=12)
    plt.xlim(0,10000) # L = 20
    plt.ylim(0.2,1.2) # L = 20
    
    plt.figure(3)
    plt.title('Heat capacity, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$C_V$',size=12)
    plt.xlim(0,10000) # L = 20
    plt.ylim(-0.5,2) # L = 20
    
    plt.figure(4)
    plt.title('Susceptibility, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\chi$',size=12)
    plt.xlim(0,30000) # L = 20
    plt.ylim(-1,12) # L = 20
    
    plt.figure(5)
    plt.title('Accepted configurations vs. MC cycles, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel('Accepted configurations',size=12)
    plt.xlim(0,3000)   # L = 20  
    plt.ylim(-10,150) # L = 20
    
    for i in range(len(temps)):
        # Read data
        filenameMC = '../Datasets/ExpectationValues_MC_%d_%.1f_%d.txt' % (L,temps[i],random)
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
        
        if i == 0:
            color = 'b'
        else:
            color = 'r'
        
        # Fill plots with values
        plt.figure(1)
        plt.plot(MC_cycles,E_avg,color,label=r'$T =$ %.2f' % temps[i])
        #coefficients = np.polyfit(MC_cycles,E_avg,6)
        #polynomial = np.poly1d(coefficients)
        #E_avg_smooth = polynomial(MC_cycles)
        #plt.plot(MC_cycles,E_avg_smooth,label=r'$T =$ %.2f' % temps[i])
        #plt.ylim(-2.2,-1.0)
        plt.legend()
        
        plt.figure(2)
        plt.plot(MC_cycles,M_absavg,color,label=r'$T =$ %.2f' % temps[i])
        #coefficients = np.polyfit(MC_cycles,M_absavg,6)
        #polynomial = np.poly1d(coefficients)
        #M_absavg_smooth = polynomial(MC_cycles)
        #plt.plot(MC_cycles,M_absavg_smooth,label=r'$T =$ %.2f' % temps[i])
        #plt.ylim(0.7,1.1)
        plt.legend()
        
        plt.figure(3)
        plt.plot(MC_cycles,C_v,color,label=r'$T =$ %.2f' % temps[i])
        #plt.ylim(-0.5,1.0)
        plt.legend()
        
        plt.figure(4)
        plt.plot(MC_cycles,X,color,label=r'$T =$ %.2f' % temps[i])
        #plt.ylim(-0.1,0.3)
        plt.legend()
        
        plt.figure(5)
        plt.plot(MC_cycles,AC,color,label=r'$T =$ %.2f' % temps[i])
        #plt.ylim(-0.5,2.0)
        plt.legend()
    
    plt.show()
    
    return


def exp_values_T(L_list,random):
    # Input: L_list = array or list of lattice sizes
    # Output: plots expectation values vs. temperature
    
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

    for i in range(len(L_list)):
        # Read file
        filenameT = '../Datasets/ExpectationValues_temp_%d_%d.txt' % (L_list[i],random)
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameT)
        
        # Fill plots with values
        plt.figure(6)
        plt.plot(T,E_avg,label='L = %d' % L_list[i])
        plt.ylim(-2.5,0)
        plt.legend()
        
        plt.figure(7)
        plt.plot(T,M_absavg,label='L = %d' % L_list[i])
        #plt.ylim(0,1.4)
        plt.legend()
        
        plt.figure(8)
        plt.plot(T,C_v,label='L = %d' % L_list[i])
        #plt.ylim(0,4)
        plt.legend(loc=2)
        
        plt.figure(9)
        plt.plot(T,X,label='L = %d' % L_list[i])
        #plt.ylim(0,250)
        plt.legend(loc=2)
        
        plt.figure(10)
        plt.plot(T,AC,label='L = %d' % L_list[i])
        plt.legend(loc=2)
    
    plt.show()
    
    return

def exp_rand_MC(L,temp,random_list):
    # Input: L = lattice size, temp = temperature, random = 0 (non-random spin matrix) or 1 (random spin matrix)
    # Ouput: plot of expectation values <E> and <M> for random and non-random spin matrix
    
    # Define plots
    plt.figure(11)
    plt.title('Energy, L = %d, T = %.1f' % (L,temp),size=12)
    plt.xlim(0,2000) # L = 20
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle E \rangle$',size=12)
    
    plt.figure(12)
    plt.title('Magnetism, L = %d, T = %.1f' % (L,temp),size=12)
    plt.xlim(0,35000) # L = 20
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle M \rangle$',size=12)
    
    for i in range(len(random_list)):
        
        # Read data
        filenameMC = '../Datasets/ExpectationValues_MC_%d_%.1f_%d.txt' % (L,temp,random_list[i])
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
        
        if i == 0:
            plot_label = 'Random start configuration'
            color = 'b'
        else:
            plot_label = 'Ground state as start'
            color = 'r'
        
        # Fill plots with values
        plt.figure(11)
        plt.plot(MC_cycles,E_avg,color,label=plot_label)
        plt.legend()
        
        plt.figure(12)
        plt.plot(MC_cycles,M_absavg,color,label=plot_label)
        plt.legend()
    
    plt.show()
        
    return
    
    
def probability(L,temp,random):
    # Input: L = lattice size, temp = temperature
    # Output: probability P(E) of a given E and variance of E
    
    filenameMC = '../Datasets/ExpectationValues_MC_%d_%.1f_%d.txt' % (L,temp,random)
    filenameE = '../Datasets/Energy_MC_%d_%.1f_%d.txt' % (L,temp,random)
    
    # Read data
    T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
    E = np.loadtxt(filenameE,unpack=True)
    
    # Round-off to make counting possible
    precision = 2 # Number of decimal points
    E_avg = np.around(E_avg,decimals=precision)
    E = np.around(E,decimals=precision)
    E = E[3000:-1] # Only include energies after most likely state is reached

    # Calculate probability
    no_elements = len(E)
    counter = list(E).count(E_avg[-1])
    prob = float(counter)/no_elements
    
    print "Probability: (T = %.2f)" % temp
    print "E_avg = %.2f" % E_avg[-1]
    print "E = %.2f" % E[-1]
    print "P(E) = %.2f" % prob
    print "Var(E) = %.2f" % E_var[-1]
    
    #print len(E[3000-2:-1]),len(MC_cycles[20+3000:-1])

    return


def T_C_estimate(L_list,random):
    # Input: L = array or list of lattice sizes
    # Output: T_C = numerical estimate of critical temperature
    
    print "Critical temperature:"
    
    nu = 1
    a = 5
    T_C_exact = 2.269 # Exact critical temperature [kT/J]
        
    print "Exact = \t %.3f" % T_C_exact
    
    T_C_estimate = np.zeros(len(L_list))
    T_max = np.zeros(len(L_list))
    a_estimate = np.zeros(len(L_list))
    
    print "\nBEFORE:"
    
    for i in range(len(L_list)):
        #filenameT = '../Datasets/ExpectationValues_temp_%d_%d.txt' % (L_list[i],random)
        filenameT = 'C:\Users\mariefoss\Documents\GitHub\Project4\Datasets\ExpectationValues_temp_%d_%d.txt' % (L_list[i],random)
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameT) # Read data
        
        # Find critical temperature from data set:
        T_max[i] = T[np.where(X == max(X))]
        print "L = ", L_list[i], ", T_max = ", T_max[i]
        
        # Make estimate for T(L -> inf) with a = 1:
        T_C_estimate[i] = T_max[i] - a*L_list[i]**(-1./nu)
        
        # Find value of "a" to make the solution equal to the exact one:
        a_estimate[i] = (T_max[i] - T_C_exact)*L_list[i]
        
<<<<<<< HEAD
        
        print "Numerical = \t %.3f \t L = %d \t a = %.2f \t T_max = %.5f" % (T_C_estimate[i],L_list[i],a_estimate[i],T_max[i])
=======
        #print "Numerical = \t %.3f \t L = %d" % (T_C_estimate[i],L_list[i])
>>>>>>> origin/master
    
    
    print "\nAFTER:"
    a_avg = np.average(a_estimate) # Average "a"-values
    #a_avg = 0.75
    
    # Calculate estimate for T(L -> inf) with a = a_avg:
    for i in range(len(L_list)):
        print "Numerical = \t %.3f \t L = %d" % (T_max[i] - a_avg*L_list[i]**(-1./nu),L_list[i])
    
    print "\na_avg = \t %.2f" % a_avg
        
    return
    

def timeusage():
    # Plot time usage
    L_sizes = np.array([10,20,40,60,80,100])
    time = np.array([57,232,923,2886,7320,31582]) # [s]
    plt.figure()
    plt.plot(L_sizes,time,'k',label='Data')
    plt.plot(L_sizes,L_sizes**2,'b--',label=r'$L^2$')
    plt.plot(L_sizes,np.exp(L_sizes/9.5),'r--',label=r'$e^{L/9.5}$')
    plt.xlabel('Lattice size',size=12)
    plt.ylabel('Time [s]')
    plt.legend()
    plt.show()
    return
    

def main(argv): # Change these values according to run!
    
    # Is the spin matrix initially random? True = 1, false = 0.
    random = 0                                 
    random_list = [0,1]
    
    # How large is the lattice?
<<<<<<< HEAD
    L = 20
=======
    L = 10
>>>>>>> origin/master
    L_list = [10,20,40,60,80,100]
    
    # What is the temperature?
    temp = 2.4
    temp_list = [1.0,2.4]
    
    # Run calculations:
<<<<<<< HEAD
    #exp_values_MC(L,temp_list,random)      # Expectation values vs. # of Monte Carlo cycles
    #exp_rand_MC(L,temp,random_list)        # Expectation values vs. # of Monte Carlo cycles: Different initial spin matrix!
    #exp_values_T(L_list,random)            # Expectation values vs. temperature
    probability(L,temp,random)             # Calculate probabilty for <E>
    #T_C_estimate(L_list,random)            # Calculate estimate of T_crit
    #timeusage()                         # Plot time usage to run code (based on dataset)
=======
    #exp_values_MC(L,temps,random)         # Expectation values vs. # of Monte Carlo cycles
    #exp_rand(L,temp,random_list)          # Expectation values vs. # of Monte Carlo cycles: Different initial spin matrix!
    #exp_values_T(L_list,random)          # Expectation values vs. temperature
    #probability(L,temp,random)           # Calculate probabilty for <E>
    T_C_estimate(L_list,random)          # Calculate estimate of T_crit
    
    #filename = '../Datasets/ExpectationValues_MC_20_1.0_1.txt'
    #print np.loadtxt(filename)
>>>>>>> origin/master
    

if __name__ == "__main__":
    main(sys.argv[1:]) 