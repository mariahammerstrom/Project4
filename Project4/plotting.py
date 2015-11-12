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
    X = data[5] 		   	 		# Susceptibility
    E_var = data[6]                 # Variance of energy
    AC = data[7]		    		# Accepted configurations
    
    return T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC


    
def exp_values_MC(L,temps,random):
    # Input: L = array or list with dimension of lattice(s), temps = array or list of temperatures,
    # random = 0 (non-random initial spin matrix) or 1 (random initial spin matrix)
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
        plt.legend()
        
        plt.figure(2)
        plt.plot(MC_cycles,M_absavg,color,label=r'$T =$ %.2f' % temps[i])
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
    # Input: L = lattice size, temp = temperature, random = 0 (non-random initial spin matrix) or 1 (random initial spin matrix)
    # Ouput: plot of expectation values <E> and <M> for random and non-random spin matrix
    
    # Define plots
    plt.figure(11)
    plt.title('Energy, L = %d, T = %.1f' % (L,temp),size=12)
    plt.xlim(0,2000) # L = 20
    #plt.ylim(-2.1,-0.6) # T = 1.0
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle E \rangle$',size=12)
    
    plt.figure(12)
    plt.title('Magnetism, L = %d, T = %.1f' % (L,temp),size=12)
    plt.xlim(0,35000) # L = 20
    #plt.ylim(0.2,1.2) # T = 1.0
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\langle M \rangle$',size=12)
    
    for i in range(len(random_list)):
        
        # Read data
        filenameMC = '../Datasets/ExpectationValues_MC_%d_%.1f_%d.txt' % (L,temp,random_list[i])
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
        
        if i == 0:
            plot_label = 'Ground state as start'
            color = 'b'
        else:
            plot_label = 'Random start configuration'
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
    
    
def probability(L,temp,random,cutoff):
    # Input: L = lattice size, temp = temperature, 
    # random = 0 (non-random initial spin matrix) or 1 (random initial spin matrix),
    # cutoff = no. of MC cycles to ignore.
    # Output: probability P(E) of a given E and variance of E.
    
    filenameMC = '../Datasets/ExpectationValues_MC_%d_%.1f_%d.txt' % (L,temp,random)
    filenameE = '../Datasets/Energy_MC_%d_%.1f_%d.txt' % (L,temp,random)
    
    # Read data
    T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameMC)
    E = np.loadtxt(filenameE,unpack=True)
    
    # Round-off to make counting possible
    precision = 2 # Number of decimal points
    E_avg = np.around(E_avg,decimals=precision)
    E = np.around(E,decimals=precision)
    E = E[cutoff:-1] # Only include energies after most likely state is reached

    # Calculate probability
    no_elements = len(E)
    counter = list(E).count(E_avg[-1])
    prob = float(counter)/no_elements
    
    print "Probability: (T = %.2f)" % temp
    print "E_avg = %.2f" % E_avg[-1]
    print "E = %.2f" % E[-1]
    print "P(E) = %.2f" % prob
    print "Var(E) = %.2f" % E_var[-1]
    
    return


def T_C_estimate(L_list,random):
    # Input: L = array or list of lattice sizes, random = 0 (non-random initial spin matrix) or 1 (random initial spin matrix)
    # Output: T_C = numerical estimate of critical temperature
    
    print "CRITICAL TEMPERATURE:"
    nu = 1
    a = 0.25
    T_C_exact = 2.269 # Exact critical temperature [kT/J]
    T_max = np.zeros(len(L_list))
        
    print "Exact = \t %.3f" % T_C_exact
    print "a = \t \t %.2f \n" % a
    
    for i in range(len(L_list)):
        filenameT = '../Datasets/ExpectationValues_temp_%d_%d.txt' % (L_list[i],random)
        T,MC_cycles,E_avg,M_absavg,C_v,X,E_var,AC = read_file(filenameT) # Read data
        
        # Find critical temperature from data set:
        T_max[i] = T[np.where(X == max(X))]
        T_estimate = T_max[i] - a*L_list[i]**(-1./nu)
        rel_error = (T_estimate - T_C_exact)/T_C_exact
        print "L = %d: \t T_max = %.3f \t T_estimate = %.3f \t Rel.error = %.2e" % (L_list[i],T_max[i],T_estimate,rel_error)
            
    return
    

def timeusage():
    # Plot time usage for different lattice sizes L
    
    L_sizes = np.array([10,20,40,60,80,100])          # lattice sizes
    time = np.array([57,232,923,2886,7320,31582])     # time usage [s]
    
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
    L = 20
    L_list = [80,100,120]
    
    # What is the temperature?
    temp = 2.4
    temp_list = [1.0,2.4]
    yes = 1
	
    while yes==1:
        print " "
        print "What do you want to do?"
        print "1 = Plot expectation values against MC cycles;"
        print "2 = Plot expectation values against MC cycles, random vs. fixed ground state;"
        print "3 = Plot expectation values against temperature;"
        print "4 = Calculate probability for <E>;"
        print "5 = Calculate estimate of critical temperature;"
        print "6 = Plot time usage;"
        number = raw_input("Choose a number: ")
        option = int(number)

        if option == 1:
	        exp_values_MC(L,temp_list,random)      # Expectation values vs. number of Monte Carlo cycles
        elif option == 2:
            exp_rand_MC(L,temp,random_list)        # Expectation values vs. number of Monte Carlo cycles: Random vs. non-random
        elif option == 3:
            exp_values_T(L_list,random)            # Expectation values vs. temperature
        elif option == 4:
            probability(L,temp,random,cutoff=3000)  # Calculate probabilty for <E>
        elif option == 5:
            T_C_estimate(L_list,random)            # Calculate estimate of T_crit based on analytical value of the constant a
        elif option == 6:
            timeusage()                            # Plot time usage to run code (based on dataset)
        else:
            print "Invalid choice. Aborting."
        new = raw_input("Want to try again? Yes(1) or No(0): ")
        yes = int(new)
        if yes==0:
            print "Bye."
        else:
            yes=1
		
		
if __name__ == "__main__":
    main(sys.argv[1:]) 