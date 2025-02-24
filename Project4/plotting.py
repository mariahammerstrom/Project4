"""
THE 2-DIMENSIONAL ISING MODEL

A program that plots:
    1) Expectation values,
    2) Number of accepted configurations, and
    3) Number of Monte Carlo cycles
as a function of temperature and number of Monte Carlo cycles.
The program also plots the time usage of the program based on a given data set.

The program calculates:
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
    #plt.xlim(0,10000) # L = 20
    #plt.ylim(0.2,1.2) # L = 20
    
    plt.figure(3)
    plt.title('Heat capacity, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$C_V$',size=12)
    #plt.xlim(0,10000) # L = 20
    #plt.ylim(-0.5,2) # L = 20
    
    plt.figure(4)
    plt.title('Susceptibility, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel(r'$\chi$',size=12)
    #plt.xlim(0,30000) # L = 20
    #plt.ylim(-1,12) # L = 20
    
    plt.figure(5)
    plt.title('Accepted configurations vs. MC cycles, L = %d' % L,size=12)
    plt.xlabel('# of MC cycles',size=12)
    plt.ylabel('Accepted configurations',size=12)
    #plt.xlim(0,3000)   # L = 20  
    #plt.ylim(-10,150) # L = 20
    
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
    E = E[cutoff:-1]
    
    # Plot probability distribution
    results, edges = np.histogram(E, normed=True)
    binWidth = edges[1] - edges[0]
    plt.bar(edges[:-1], results*binWidth, binWidth)
    plt.xlabel(r'$E$')
    plt.ylabel(r'$P(E)$')
    plt.title(r'Probability distribution, $L$ = %d, $T$ = %.2f, $<E>$ = %.2f' % (L,temp,E_avg[-1]),size=12)
    plt.show()
    
    print "T = %.2f" % temp
    print "E_avg = %.2f" % E_avg[-1]
    print "Var(E) = %.2f" % E_var[-1]
    
    print "\nSum(P(E)) = ", sum(results*binWidth) # Check that the probabilities sum to 1.0
    
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
	

def string2array(list):
    list = list.split(',')
    size = len(list)
    array = np.zeros(size)
    for i in range(size):
	    array[i] = float(list[i])
    return array
    

def main(argv):

    # Command line options
    yes = 1
    random_list = [0,1]
    temp_list = [1.0,2.4]
    
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

        if option == 1: # Expectation values vs. number of Monte Carlo cycles
            l = raw_input("Which L do you want to plot for? Choose 10,20,40,60,80 or 100: ")
            L = int(l)
            rand = raw_input("Do you want to start with a random ground state? Yes(1), No(0): ")
            random = int(rand)
            exp_values_MC(L,temp_list,random)      
        elif option == 2: # Expectation values vs. number of Monte Carlo cycles: Random vs. non-random
            l = raw_input("Which L do you want to plot for? Choose 10,20,40,60,80 or 100: ")
            L = int(l)
            temperature = raw_input("T=1.0 or T=2.4? ")
            temp = float(temperature)
            exp_rand_MC(L,temp,random_list)        
        elif option == 3: # Expectation values vs. temperature
            list = raw_input("Which L do you want to plot for? Choose one or more from 10,20,40,60,80,100. Separate by commas: ")
            L_list = string2array(list)
            rand = raw_input("Do you want to start with a random ground state? Yes(1), No(0): ")
            random = int(rand)
            exp_values_T(L_list,random)            
        elif option == 4: # Calculate probabilty for <E>
            l = raw_input("Which L do you want to plot for? Choose 10,20,40,60,80 or 100: ")
            L = int(l)
            temperature = raw_input("T=1.0 or T=2.4? ")
            temp = float(temperature)
            rand = raw_input("Do you want to start with a random ground state? Yes(1), No(0): ")
            random = int(rand)
            probability(L,temp,random,cutoff=100)  
        elif option == 5: # Calculate estimate of T_crit based on analytical value of the constant a
            list = raw_input("Which L do you want to estimate T_C for? Choose one or more from 10,20,40,60,80,100. Separate by commas: ")
            L_list = string2array(list)
            rand = raw_input("Do you want to start with a random ground state? Yes(1), No(0): ")
            random = int(rand)
            T_C_estimate(L_list,random)            
        elif option == 6: # Plot time usage to run code (based on dataset)
            timeusage()                            
        else:
            print "Invalid choice. Aborting."
        new = raw_input("Want to try again? Yes(1) or No(0): ")
        yes = int(new)
        if yes==0:
            print "Bye!"
        else:
            yes=1
		
		
if __name__ == "__main__":
    main(sys.argv[1:]) 