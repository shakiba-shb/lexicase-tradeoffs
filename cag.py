from ec_ecology_toolbox import community_assembly_graph
from ec_ecology_toolbox import selection_probabilities as eco
from ec_ecology_toolbox.community_assembly_graph.example_nodes import Set_Node
import numpy as np
import random
import json
import os
import pandas as pd
from itertools import combinations
from math import ceil
from copy import deepcopy

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

def contains_optimum(cag):
    opt = False
    for n in cag.keys():
        for i in n:
            i = list(i)
            if (i.count(4) == 1 and np.sum(i) == 4):
                opt = True
    if opt:
        return True

def adjacent_nodes(n):
    adjacent_nodes = set()
    for i in range(len(n)):
        for adj in (-1, 1):
            temp_node = deepcopy(n)
            temp_node[i] = temp_node[i] + adj
            temp_node[i] = max(0, min(temp_node[i], 4))
            adjacent_nodes.add(tuple(temp_node))
    if tuple(n) in adjacent_nodes:
        adjacent_nodes.discard(tuple(n))
    return adjacent_nodes

# df = pd.DataFrame()
# directory = '/home/shakiba/lexicase-tradeoffs/S_Dim_ep0_fig7'

# for filename in os.listdir(directory):
#     if filename.endswith('.json'):
#         file_path = os.path.join(directory, filename)
#         with open(file_path, 'r') as f:
#             data = json.load(f)
#             data['last_pop'] = [tuple(l) for l in data['last_pop']]
#             df = pd.concat([df, pd.DataFrame(data, index=[0])], ignore_index=True)

def my_custom_user_func(pop):
    result = dict()
    adj = set()
    for gene in pop:
        adj.update(adjacent_nodes(list(gene)))

    for gene in adj:
        node = Set_Node(pop)
        #node = deepcopy(pop)
        #node = [list(i) for i in node]
        members = [list(i) for i in node.members]
        members.append(list(gene))
        #members = [[0, 1, 0], [2, 1, 0], [0, 0, 1], [0, 2, 1], [1, 0, 1], [1, 1, 0], [2, 1, 2], [1, 2, 1], [0, 2, 0], [0, 0, 0], [1, 1, 2], [1, 0, 0], [2, 0, 0], [2, 2, 0], [0, 1, 1], [2, 1, 1], [1, 2, 0], [1, 3, 1], [3, 1, 0], [1, 0, 2], [1, 1, 1]]

        phenotypes = []
        for i in members:
            phenotypes.append(list(fitness_function(i, damp=1)))
        prob = eco.LexicaseFitness(phenotypes, epsilon)
        P_survival = list((np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G)
        #P_survival = [0.3]*len(prob)

        survivors = []
        for m, p in zip(members, P_survival):
            if p > t:
                survivors.append(m)

        if survivors == []:
            x = int(1/(1 - (1 - t**(1/G))**(1/S)))
            all_combinations = list(combinations(members, x))
            for c in all_combinations:
                new_node = frozenset(tuple(i) for i in c)
                if new_node != node.members:
                    if new_node in result:
                        result[new_node] += mu/(len(adj)*len(all_combinations))
                    else:
                        result[new_node] = mu/(len(adj)*len(all_combinations))
                        #result[tuple(tuple(i) for i in all_combinations)] = (mu/(len(adj)*len(all_combinations)))
        else:           
            new_node = frozenset(tuple(i) for i in survivors)
            if new_node != node.members:
                if new_node in result:
                    result[new_node] += mu/len(adj)
                else:
                    result[new_node] = mu/len(adj)
            
    if len(result) == 0:
        node.is_sink = True
    else:
        node.is_sink = False
    return result


last_pop = [[1, 1, 0], [0, 2, 0], [0, 1, 1]]
#last_pop  = [[0,0,0]]
S = 100
G = 100
t = 0.5
mu = 0.1
epsilon = 0
max_graph_size = 5
last_pop = frozenset(tuple(i) for i in last_pop)
#print(my_custom_user_func(last_pop))


#print(adjacent_nodes(start_node))
#CAG = community_assembly_graph.CAG(1, custom_user_func, 3, False)
#print(adjacent_nodes(start_node))
# for start_node in last_pop:
#     start_node = tuple(start_node)
CAG = community_assembly_graph.CAG(last_pop, my_custom_user_func, max_graph_size, False)
print("CAG = ",CAG)
print("len(CAG) = ", len(CAG))

if contains_optimum(CAG):
    print("optimum is reachable")
else:
    if(len(CAG) < max_graph_size):
        print("optimum is not reachable")

    elif(len(CAG) >= max_graph_size):
        if contains_optimum(CAG):
            print("optimum is reachable")
        else:
            print("optimum isn't easily reachable")