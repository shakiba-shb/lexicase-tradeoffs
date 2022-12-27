import example
import numpy as np
import random
from random import choices
from numpy.random import rand
import matplotlib.pyplot as plt
import time
import scipy.stats as st

#parameters
dim = 15 #length of genomes
G = 1000 #number of generations
S = 100 #population size 
n_iters = 10
damp = 2
#MU = [0.05, 0.3, 0.7]
MU = [0.01, 0.05, 0.1, 0.3, 0.5, 0.7]
p_thresh = 0.5
from numpy.random import rand

def fitness_function(x, damp):
    #Genotype is a vector of 5 integers between 0 and 5.
    #Phenotype is a vector where each number in the original vector has all of the other numbers subtracted from it.
    #Fitness is the highest value in the phenotype.
    x_copy = np.copy(x)
    s = sum(x_copy)
    for i in range(len(x_copy)):
        x_copy[i] = x_copy[i] - s/damp + x_copy[i]
    return x_copy

def mutants(pop, mu, m):
    # Gets a population and returns all its m mutants with mutation rate mu.
    new_pop = []
    for p in pop:
        p_copy = p.copy()
        if (rand() < mu):
            I = choices(list(np.arange(0, dim)), k=m)
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
    # Counts number of 4s in genome
    counter = 0
    for i in x:
        if (i == 4):
            counter = counter+1
    return counter

N_loops = []
std = []
error_min = []
error_max = []
opt = []

for mu in MU:
    print("mu = ", mu)
    n_loops = []
    opt_tracker = []
    zeros = []
    for it in range(n_iters):  
        #print("it = ", it)
        terminate = False
        #initial_pop = [[random.randint(1,3) for i in range(dim)]]
        initial_pop = [[0 for i in range (dim)]]
        counter = 0 # indicates termination of while loop
        opt_counter = 0
        prob = []
        while (counter < 10000):
            print("Loop" , counter, "iter", it, "mu", mu, flush = True)

            muts = mutants(initial_pop, mu, 2)
            new_pop = initial_pop + muts
            for i in range(len(new_pop)):
                new_pop[i] = tuple(new_pop[i])
            new_pop = set(new_pop)   
            new_pop = [list(ele) for ele in new_pop]

            phenotypes = []
            for p in new_pop:
                phenotypes.append(list(fitness_function(p, damp)))

            prob = example.LexicaseFitness(phenotypes)
                #df = {"mu": mu,"iter": it, "loop": counter, "genome":temp[p].geno, "time":tac-tic, "prob":prob[p]}

            P_survival = (np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G

            survivors = []
            for p in range(len(new_pop)):
                if (P_survival[p] >= p_thresh):
                    survivors.append(new_pop[p])
                        
            for s in survivors:
                if (four_counter(s) == 1 and np.sum(s) == 4):
                    opt_counter = opt_counter + 1
                    #print("Optimum found at mu: ", mu, "iter: ", it, "loop: ", counter)
                if (four_counter(s) >= 2):
                    terminate = True
                            
            if (survivors == []):
                print(" No survivors at iteration ", it, "at loop ", counter, "for mu = ", mu)
            
            print ("#survivors = ", len(survivors))
            print("survivors = ", survivors)
            print("prob = ", prob)
            #time.sleep(0.05)
            initial_pop = survivors    
            counter = counter + 1
            
            if (terminate == True):
                break    

            #print("loop ", counter, "time = ", tac0 - tic0 )  
 
        n_loops.append(counter)
        opt_tracker.append(opt_counter)
        #print("n_loops = ", n_loops)


    #error_min.append(list(st.t.interval(alpha=0.95, df=len(n_loops)-1, loc=np.mean(n_loops), scale=st.sem(n_loops)))[0])
    #error_max.append(list(st.t.interval(alpha=0.95, df=len(n_loops)-1, loc=np.mean(n_loops), scale=st.sem(n_loops)))[1])
    std.append(np.std(n_loops))

    N_loops.append(np.average(n_loops))
    opt.append(np.average(opt_tracker))

print("N_loops = ", N_loops)
print("opt = ", opt)
#error = [error_min, error_max]

plt.figure("number of loops to fail")
plt.plot(MU, N_loops)
plt.errorbar(MU, N_loops, yerr = np.multiply(0.7, std), fmt = 'o', c = 'red', capsize = 10)
plt.title("number of loops to fail")
plt.xlabel("mutation rate")
plt.ylabel("number of loops")  
plt.savefig("mutation rate vs number of loops.png")