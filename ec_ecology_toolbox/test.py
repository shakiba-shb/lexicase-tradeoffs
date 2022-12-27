import example
import numpy as np

def fitness_function(x):
    #Genotype is a vector of 5 integers between 0 and 5.
    #Phenotype is a vector where each number in the original vector has all of the other numbers subtracted from it.
    #Fitness is the highest value in the phenotype.
    x_copy = np.copy(x)
    s = sum(x_copy)
    for i in range(len(x_copy)):
        x_copy[i] = x_copy[i] - s/2 + x_copy[i]
    return x_copy

pop = [[0, 1, 0, 0, 0], [2, 1, 0, 0, 0 ]]

phenotypes = []
for p in pop:
    phenotypes.append(list(fitness_function(p)))

print(phenotypes)
#print(fitness_function(phenotypes))
#print(example.LexicaseFitness(phenotypes))