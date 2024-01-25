from ec_ecology_toolbox import selection_probabilities as eco
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
import math
import csv

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
        indv_copy = indv.copy()
        
        if per_genome:
            # Per-genome mutation
            for _ in range(m):
                gene_index = random.randint(0, len(indv_copy) - 1)
                if random.random() < mu:
                    if gene_index not in mutated_indices:
                        indv_copy[gene_index] += random.choice([-1, 1])
                        indv_copy[gene_index] = max(0, min(indv_copy[gene_index], 4))
                        mutated_indices.add(gene_index)
        else:
            # Per-site mutation
            order = random.sample(list(np.arange(0, len(indv_copy))), len(indv_copy))
            for _ in range(m):
                for gene_index in order:
                    if (random.random()) < mu and (gene_index not in mutated_indices):
                        indv_copy[gene_index] += random.choice([-1, 1])
                        indv_copy[gene_index] = max(0, min(indv_copy[gene_index], 4))
                        mutated_indices.add(gene_index)
                        break

        new_pop.append(indv_copy)
    return new_pop

def all_populations_identical(population_dict):
    # Check if there are no generations in the dictionary
    if not population_dict:
        print("pop_history empty")
        return False

    # Get the first generation's population
    first_generation_population = population_dict.get(min(population_dict.keys()))

    # Compare all other generations' populations to the first generation
    for generation, population in population_dict.items():
        if population != first_generation_population:
            return False  # Populations are not identical

    # If all populations are identical
    return True

def experiment (S = None, G = None, Dim = None, Damp = None, MU = None, Seed = None, epsilon = None, max_loops = None , rdir = "" ):

    runid = uuid.uuid4()
    print("G = ", G, "S = ", S, "Dim = ", Dim, 'Damp = ', Damp, 'MU = ', MU, 'Seed = ', Seed, 'max_loops =', max_loops, 'epsilon =', epsilon)
    terminate = False
    initial_pop = [[0 for i in range (Dim)]] #initial population
    counter = 0 #indicates termination of while loop
    fail_reason = ''
    last_pop = []
    pop_history = dict()
    p_thresh = 0.5

    while (counter < max_loops):

        muts = mutants(initial_pop, MU, 1) #create mutants of current genotypes in the population
        new_pop = initial_pop + muts #add mutants to current population
        new_pop = list(set(tuple(genotype) for genotype in new_pop)) #remove similar genotypes
        new_pop = [list(genotype) for genotype in new_pop]

        phenotypes = [] 
        for genome in new_pop: #create phenotypes from genotypes based on defined fitness function
            phenotypes.append(list(fitness_function(genome, Damp)))
        prob = eco.LexicaseFitness(phenotypes, epsilon) #calculate probablity of being selected by lexicase selection for all phenotypes
        P_survival = list((np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G) #calculate probability of survival based on equation(?) for all phenotypes

        for p in prob: 
            assert 0 <= p <= 1, f"invalid probability value: prob= {prob}, phenotypes= {phenotypes}"
            
        survivors = []
        for pheno in range(len(new_pop)): #find survivors based on probability of survival
            if (P_survival[pheno] >= p_thresh):
                survivors.append(new_pop[pheno])
        
        if (survivors == []): #define what happens if no individual survives
            print(" No survivors at seed ", Seed, "loop ", counter, "for mu = ", MU, "G = ", G, "S = ", S, "Dim = ", Dim)
            #print(len(new_pop))
            # for pheno in range(len(new_pop)): 
            #     if (P_survival[pheno] >= np.mean(P_survival)):
            #         survivors.append(new_pop[pheno])
            # if (survivors == []):
            #     for pheno in range(len(new_pop)):
            #         if(math.isclose(P_survival[pheno],np.mean(P_survival))):
            #             survivors.append(new_pop[pheno])
            x = int(1/(1 - (1 - p_thresh**(1/G))**(1/S)))
            sample = random.sample(new_pop, x)
            phenotypes = [] 
            for genome in sample: 
                phenotypes.append(list(fitness_function(genome, Damp)))
            prob = eco.LexicaseFitness(phenotypes, epsilon) 
            P_survival = list((np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G) 
            
            for pheno in range(len(sample)): 
                if (P_survival[pheno] >= p_thresh):
                    survivors.append(new_pop[pheno])

        for s in survivors:
        #Look for optimums in the population 
            if (s.count(4) == 1 and np.sum(s) == 4):
                #print("Optimum found at loop", counter)
                print("optimum found at loop ", counter)
                terminate = True

        # if (counter == max_loops - 1): #Check if we're stuck at all-zeros
        #     if((len(survivors) == 1) and (not np.any(survivors[0]))):
        #         print("stuck at all zeros")
        #         terminate = True

        # if (counter > max_loops - 1000):
        #     pop_history[counter] = survivors
        
        if (terminate == True):
            break

        initial_pop = survivors #set current survivers as initial population of the next loop   
        if (survivors == []):
            print("no one is going to the next loop")
        counter = counter + 1

    #if (counter == max_loops):
        # print("History checking...")
        # if all_populations_identical(pop_history):
        #     fail_reason = 'stuck'
        # else:
        #     fail_reason = 'searching'
        # last_pop = list(pop_history.values())[-1]
    last_pop = survivors

    new_row = {'G': G, 'S':S, 'Dim':Dim, 'Damp':Damp, 'MU': MU, 'Seed': int(Seed), 'max_loops': max_loops, 'epsilon': epsilon, 'fail_loop':counter, 'last_pop': last_pop}
    
    filename = rdir + f'/runid-{runid}.json'
    
    with open(filename, 'w') as of:
        json.dump(new_row, of)
    
    return new_row

# import fire
# if __name__ == '__main__':
#     fire.Fire(experiment)

import argparse
if __name__ == '__main__':

    # parse command line arguments
    parser = argparse.ArgumentParser(
        description="Run a single experiment with specified variables.", add_help=False)
    parser.add_argument('-S', action='store', type=int, default=100,
                        help='Population size')
    parser.add_argument('-G', action='store', type=int, default=100,
                        help='G')                
    parser.add_argument('-Dim', action='store', type=int, default=5,
                        help='Dimension size')
    parser.add_argument('-Damp', action='store', type=int, default=1,
                        help='Dampening factor')
    parser.add_argument('-MU', action='store', type=float, default=0.01,
                        help='mutation rate')
    parser.add_argument('-Seed', action='store', type=int, default=42,
                        help='seed')                    
    parser.add_argument('-epsilon', action='store', type=float, default=0,
                        help='epsilon')
    parser.add_argument('-max_loops', action='store', type=int, default=1000,
                        help='maximum number of loops')    
    parser.add_argument('-rdir', action='store', default='results', type=str,
                        help='Name of save file')                                   
    
    args = parser.parse_args()
    experiment( args.S, args.G, args.Dim, args.Damp, args.MU, args.Seed, args.epsilon, args.max_loops, args.rdir )