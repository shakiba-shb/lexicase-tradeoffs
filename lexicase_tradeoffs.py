import my_module
import numpy as np
import random
from random import choices
from numpy.random import rand
import time
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import bootstrap
import pandas as pd
import seaborn as sns
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


#parameters
max_loops = 10000
dim = [5] #dimension of genomes
G = [500] #number of generations
S = [500] #population size 
n_iters = 30 #number of iterations
damp = [1] #dampening factor
MU = [0.002] #mutation rate
p_thresh = 0.5 #Threshold probability of survival of a genome

def fitness_function(x, damp):
    #Function used for building phenotypes based on genotypes.
    #Genotype is a vector of dim integers between 0 and 4.
    #Phenotype is a vector where each number in the original vector has all of the other numbers subtracted from it.
    #Fitness is the highest value in the phenotype.
    x_copy = list(np.copy(x))
    s = sum(x_copy)
    for i in range(len(x_copy)):
        x_copy[i] = x_copy[i] - s/damp + x_copy[i]/damp
    return x_copy

def mutants(pop, mu, m):
    #Gets a population and returns all its m mutants with mutation rate mu.
    new_pop = []
    for p in pop:
        p_copy = p.copy()
        if (rand() < mu):
            I = choices(list(np.arange(0, d)), k=m)
            for i in I:           
                if (random.random() < 0.5):
                    p_copy[i] = p_copy[i] + 1
                else:
                    p_copy[i] = p_copy[i] - 1
                if(p_copy[i] > 4): #cap genes at 4
                    p_copy[i] = 4
                if(p_copy[i] < 0): #cap genes at 0
                    p_copy[i] = 0
        new_pop.append(p_copy)
    return new_pop

def four_counter(x):
    #Counts number of 4s in genome (4 is the maximum value a genome can take)
    counter = 0
    for i in x:
        if (i == 4):
            counter = counter+1
    return counter

N_loops = []
error_min = []
error_max = []
opt = []
df = pd.DataFrame(columns = ['G', 'S', 'dim', 'damp', 'MU', 'mean_fail_loop', 'p_fail'])

for g in G:
    for ss in S:
        for d in dim:
            for dd in damp:
                for mu in MU:
                    print("mu = ", mu, "G = ", g, "S = ", ss, 'damp = ', dd, 'dim = ', d)
                    time.sleep(3)
                    n_loops = [] #keeps number of loops until termination in each iteration
                    opt_tracker = [] #keeps track of number of optimums found

                    for it in range(n_iters):  
                        #print("it = ", it)
                        terminate = False
                        initial_pop = [[0 for i in range (d)]] #initial population
                        counter = 0 #indicates termination of while loop
                        opt_counter = 0
                        prob = [] #keeps probabilities of genomes
                        
                        while (counter < max_loops):
                            
                            muts = mutants(initial_pop, mu, 1) #create mutants of current genotypes
                            new_pop = initial_pop + muts #add mutants to population
                            for i in range(len(new_pop)): #remove similar genotypes
                                new_pop[i] = tuple(new_pop[i])
                            new_pop = set(new_pop)   
                            new_pop = [list(ele) for ele in new_pop]

                            phenotypes = [] 
                            for p in new_pop: #create phenotypes from genotypes based on defines fitness function
                                phenotypes.append(list(fitness_function(p, dd)))

                            prob = my_module.LexicaseFitness(phenotypes) #calculate probablity of being selected by lexicase selection
                            P_survival = list((np.ones(len(prob)) - (np.ones(len(prob)) - prob)**ss)**g) #probability of survival based on equation(3)

                            survivors = []
                            for p in range(len(new_pop)): #find survivors based on probability of survival
                                if (P_survival[p] >= p_thresh):
                                    survivors.append(new_pop[p])
                            
                            if (survivors == []): #define what happens if no individual survives
                                print(" No survivors at iteration ", it, "at loop ", counter, "for mu = ", mu, "G = ", g, "S = ", ss)
                                for p in range(len(new_pop)):
                                    if (P_survival[p] >= np.mean(prob)):
                                        survivors.append(new_pop[p])            
                            
                            for s in survivors:
                                #Look for optimums in the population 
                                if (four_counter(s) == 1 and np.sum(s) == 4):
                                    print("Optimum found at loop", counter)
                                    terminate = True
                                                
                                #if (four_counter(s) >= 2):
                                    #print("found 2 fours at loop ", counter)
                                    #terminate = True

                                #Look for local optima where population gets stuck
                                if (counter == max_loops - 1):
                                    if((len(survivors) == 1) and (not np.any(s))):
                                        print("stuck at all zeros")
                                        terminate = True
                            
                            #print("survivors = ", survivors)
                            #time.sleep(0.01)

                            if (terminate == True):
                                break

                            initial_pop = survivors #set current survivers as initial population of the next loop   
                            counter = counter + 1
                
                        n_loops.append(counter)
                        opt_tracker.append(opt_counter)
                        print("iter = ", it, "n_loops = ", n_loops[it])

                    p_fail = 1 - ((sum(1 for i in n_loops if i < max_loops - 1))/n_iters) #probability of failure
                    N_loops.append(np.average(n_loops))
                    #opt.append(np.average(opt_tracker))
                    print('p_fail = ', p_fail)
                    new_row = {'G': g, 'S':ss, 'dim':d, 'damp':dd, 'MU': (mu*d)/10, 'mean_fail_loop':np.mean(n_loops), 'p_fail':p_fail}
                    df = df.append(new_row, ignore_index = True)

                    data = (n_loops,)    
                    bootstrap_ci = bootstrap(data, np.mean, confidence_level=0.95, random_state=1, method='percentile')
                    error_min.append(bootstrap_ci.confidence_interval[0])
                    error_max.append(bootstrap_ci.confidence_interval[1])


df.to_csv("out.csv", index=False)

##Use following code to create heatmaps of any 2 parameters:
# df = df.pivot('param1', 'param2', 'p_fail')
# plt.figure("heat map of parameters")
# plt.title("probabilty of failure")
# ax = sns.heatmap(df, cmap = 'coolwarm')
# ax.invert_yaxis()
# plt.savefig("heatmap.png")

##Use following code to plot number of loops to fail for each mutation rate:
#error = [error_min, error_max]
#plt.figure("number of loops to fail")
#N_loops = np.array(N_loops)
#error = np.array(error)
##plt.errorbar(MU, N_loops, yerr = np.abs(N_loops - error) , fmt = 'o', c = 'red', capsize = 5)
#plt.title("number of loops to fail")
#plt.xlabel("mutation rate")
#plt.ylabel("number of loops")  
#plt.savefig("mutation rate vs number of loops.png")
