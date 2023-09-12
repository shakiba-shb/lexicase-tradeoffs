import ec_ecology_toolbox as eco
import numpy as np
import random
import time
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.stats import bootstrap
import pandas as pd
import seaborn as sns
import uuid
import json
import sys
import os

def fitness_function(x, damp):
    """
    Fitness function based on Equation (1)

    Parameters:
    - x (list): genotype
    - damp (float): dampening factor
    
    Returns:
    - y (list): phenotype
    """
    x = np.array(x)
    s = np.sum(x)
    return x - s / damp + x / damp

def mutants(pop, mu, m, per_genome = True):
    """
    Generates mutants from a population.

    Parameters:
    - pop (list of lists): input population
    - mu (float): mutation rate
    - m (int): The number of mutations that could happen (single mutation or more).
    - per_genome (bool): if True, returns per_genome mutants otherwise returns per-site mutants.
    """
    new_pop = []
    for indv in pop:
        mutations_count = 0  # Track the number of mutations
        mutated_indices = set()  # Track the indices that have been mutated
        indv = indv.copy()
        
        if per_genome:
            # Per-genome mutation
            for _ in range(m):
                gene_index = random.randint(0, len(indv) - 1)
                if random.random() < mu:
                    if gene_index not in mutated_indices:
                        indv[gene_index] += random.choice([-1, 1])
                        indv[gene_index] = max(0, min(indv[gene_index], 4))
                        mutated_indices.add(gene_index)
        else:
            # Per-site mutation
            order = random.sample(list(np.arange(0, len(indv))), len(indv))
            for _ in range(m):
                for gene_index in order:
                    if (random.random()) < mu and (gene_index not in mutated_indices):
                        indv[gene_index] += random.choice([-1, 1])
                        indv[gene_index] = max(0, min(indv[gene_index], 4))
                        mutated_indices.add(gene_index)
                        break

        new_pop.append(indv)
    return new_pop

def experiment (G = None, S = None, Dim = None, Damp = None, MU = None, Seed = None, epsilon = None, max_loops = None , p_thresh = None, rdir = "" ):

    runid = uuid.uuid4() 
    print("G = ", G, "S = ", S, "Dim = ", Dim, 'Damp = ', Damp, 'MU = ', MU, 'Seed = ', Seed, 'max_loops =', max_loops, 'epsilon =', epsilon)
    terminate = False
    initial_pop = [[0 for i in range (Dim)]] #initial population
    counter = 0 #indicates termination of while loop
    prob = [] #keeps probabilities of genomes

    while (counter < max_loops):

        muts = mutants(initial_pop, MU, 1) #create mutants of current genotypes in the population
        new_pop = initial_pop + muts #add mutants to current population
        for i in range(len(new_pop)): #remove similar genotypes
            new_pop[i] = tuple(new_pop[i])
        new_pop = set(new_pop)   
        new_pop = [list(ele) for ele in new_pop]

        phenotypes = [] 
        for genome in new_pop: #create phenotypes from genotypes based on defined fitness function
            phenotypes.append(list(fitness_function(genome, Damp)))

        prob = eco.LexicaseFitness(phenotypes, epsilon) #calculate probablity of being selected by lexicase selection for all phenotypes
        P_survival = list((np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G) #calculate probability of survival based on equation(?) for all phenotypes
        
        survivors = []
        for pheno in range(len(new_pop)): #find survivors based on probability of survival
            if (P_survival[pheno] >= p_thresh):
                survivors.append(new_pop[pheno])
        
        if (survivors == []): #define what happens if no individual survives, look at average probabilities instead of p_thresh
            print(" No survivors at seed ", Seed, "loop ", counter, "for mu = ", MU, "G = ", G, "S = ", S, "Dim = ", Dim)
            for pheno in range(len(new_pop)): 
                if (P_survival[pheno] >= np.mean(P_survival)):
                    survivors.append(new_pop[pheno])

        for s in survivors:
        #Look for optimums in the population 
            if (s.count(4) == 1 and np.sum(s) == 4):
                #print("Optimum found at loop", counter)
                terminate = True

            if (counter == max_loops - 1): #Check if we're stuck at all-zeros
                if((len(survivors) == 1) and (not np.any(s))):
                    print("stuck at all zeros")
                    terminate = True

        if (terminate == True):
            break

        initial_pop = survivors #set current survivers as initial population of the next loop   
        counter = counter + 1

    new_row = {'G': G, 'S':S, 'Dim':Dim, 'Damp':Damp, 'MU': MU, 'Seed': int(Seed), 'fail_loop':counter}
    
    filename = rdir + f'/runid-{runid}.json'
    
    with open(filename, 'w') as of:
        json.dump(new_row, of)
    
    return new_row

import fire
if __name__ == '__main__':
    fire.Fire(experiment)