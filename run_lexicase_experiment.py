import ec_ecology_toolbox as eco
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
import uuid
import json
import sys
import os
import json_numpy

def fitness_function(x, damp):
    """
    Fitness function based on Equation (1)

    Parameters:
    - x (list): genotype
    - damp (float): dampening factor
    
    Returns:
    - y (list): phenotype
    """
    y = [x_i - sum(x_j/damp for x_j in x if x_j != x_i) for x_i in x]
    return y

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
            while mutations_count < m:
                if random.random() < mu:
                    gene_index = random.randint(0, len(indv) - 1)
                    if gene_index not in mutated_indices:
                        indv[gene_index] += random.choice([-1, 1])
                        indv[gene_index] = max(0, min(indv[gene_index], 4))
                        mutations_count += 1
                        mutated_indices.add(gene_index)
        else:
            # Per-site mutation
            order = random.sample(list(np.arange(0, len(indv))), len(indv))
            while mutations_count < m:
                for gene_index in order:
                    if random.random() < mu:
                        indv[gene_index] += random.choice([-1, 1])
                        indv[gene_index] = max(0, min(indv[gene_index], 4))
                        mutations_count += 1
                        break

        new_pop.append(indv)
    return new_pop

def experiment (G = None, S = None, Dim = None, Damp = None, MU = None, Seed = None):

    runid = uuid.uuid4()
    print("G = ", G, "S = ", S, "Dim = ", Dim, 'Damp = ', Damp, 'MU = ', MU, 'Seed = ', Seed)

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

        prob = eco.LexicaseFitness(phenotypes) #calculate probablity of being selected by lexicase selection for all phenotypes
        P_survival = list((np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G) #calculate probability of survival based on equation(?) for all phenotypes
        
        survivors = []
        for pheno in range(len(new_pop)): #find survivors based on probability of survival
            if (P_survival[pheno] >= p_thresh):
                survivors.append(new_pop[pheno])
        
        if (survivors == []): #define what happens if no individual survives, look at average probabilities instead of p_thresh
            print(" No survivors at seed ", Seed, "loop ", counter, "for mu = ", MU, "G = ", G, "S = ", S)
            for pheno in range(len(new_pop)): 
                if (P_survival[pheno] >= np.mean(P_survival)):
                    survivors.append(new_pop[pheno])

        for s in survivors:
        #Look for optimums in the population 
            if (s.count(4) == 1 and np.sum(s) == 4):
                print("Optimum found at loop", counter)
                terminate = True

            if (counter == max_loops - 1): #Check if we're stuck at all-zeros
                if((len(survivors) == 1) and (not np.any(s))):
                    print("stuck at all zeros")
                    terminate = True

        if (terminate == True):
            break

        initial_pop = survivors #set current survivers as initial population of the next loop   
        counter = counter + 1

    #new_row = {'G': [G], 'S':[S], 'Dim':[Dim], 'Damp':[Damp], 'MU': [MU], 'Seed': [Seed], 'fail_loop':[counter]}
    #res = pd.concat([res, pd.DataFrame([new_row])], ignore_index=True)
    new_row = {'G': G, 'S':S, 'Dim':Dim, 'Damp':Damp, 'MU': MU, 'Seed': int(Seed), 'fail_loop':counter}
    
    filename = rdir + f'/runid-{runid}.json'
    
    with open(filename, 'w') as of:
        json.dump(new_row, of)
    
    return new_row


################################################################################
# Save data
################################################################################
if len(sys.argv) > 1:
    rdir = sys.argv[1]
else:
    from datetime import datetime
    now = datetime.now()
    dt_string = now.strftime("%Y-%m-%d_%H-%M%S")
    rdir = f'results_{dt_string}'
    print('rdir: ',rdir)
os.makedirs(rdir,exist_ok=False)

################################################################################
# run experiment
################################################################################
from pqdm.processes import pqdm
from tqdm import tqdm
import itertools as it

n_iters=10
N_JOBS=64
args=[]

G = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
S = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
Dim = [5]
Damp = [1]
MU = [0.1]
seeds = np.random.randint(1,2**15,n_iters) 
max_loops = 1000
p_thresh = 0.5
results=pd.DataFrame(columns = ['G', 'S', 'Dim', 'Damp', 'MU', 'fail_loop', 'Seed'])

# construct a list of arguments
for g,s,dim,damp,mu,seed in it.product(G, S, Dim, Damp, MU, seeds):
    args.append(
        {'G':g,
         'S':s,
         'Dim':dim,
         'Damp':damp,
         'MU': mu, 
         'Seed': seed
        }
    )
# run the experiment
print(len(args),'experiments')
DEBUG=False
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    if DEBUG:
        #################
        # serial version for debugging
        results=[]
        for a in tqdm(args):
            results.append(experiment(**a)) 
            #results = pd.concat([results, pd.DataFrame(experiment(**a))], ignore_index=True)
        #results.to_csv('~/lexicase-tradeoffs/out.csv', index = False)
    else:
        #################
        # run in parallel. set n_jobs to the number of parallel jobs you want to run 
        #res= pqdm(args, experiment, n_jobs=min(len(args),N_JOBS), argument_type='kwargs') 
        #for r in res:
        #    results = pd.concat(([results, r]), ignore_index=True)
        #results.to_csv('~/lexicase-tradeoffs/out.csv', index = False)

        results= pqdm(args, experiment, n_jobs=min(len(args),N_JOBS), argument_type='kwargs') 
        results = [r for r in results if r != None]

print('done!')