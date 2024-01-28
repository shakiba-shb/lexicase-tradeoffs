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
        for i in n.members:
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

def my_custom_user_func(pop):
    result = dict()
    adj = set()
    for gene in pop.members:
        adj.update(adjacent_nodes(list(gene)))

    for gene in adj:
        node = deepcopy(pop)
        #node = deepcopy(pop)
        #node = [list(i) for i in node]
        members = [list(i) for i in node.members]
        members.append(list(gene))
        members = list(set(tuple(m) for m in members)) #remove similar genotypes
        members = [list(m) for m in members]
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
                new_node = Set_Node([tuple(i) for i in c])
                if new_node != node:
                    if new_node in result:
                        result[new_node] += mu/(len(adj)*len(all_combinations))
                    else:
                        result[new_node] = mu/(len(adj)*len(all_combinations))
                        #result[tuple(tuple(i) for i in all_combinations)] = (mu/(len(adj)*len(all_combinations)))
        else:           
            new_node = Set_Node([tuple(i) for i in survivors])
            if new_node != node:
                if new_node in result:
                    result[new_node] += mu/len(adj)
                else:
                    result[new_node] = mu/len(adj)
            
    if len(result) == 0:
        pop.is_sink = True
    else:
        pop.is_sink = False
    return result

directory = '/home/shakiba/lexicase-tradeoffs/results'
df = pd.DataFrame(columns=['G', 'S', 'Dim', 'Damp', 'MU', 'Seed', 'max_loops', 'epsilon', 'fail_loop', 'cag_result', 'cag_size'])
for filename in os.listdir(directory):
    if filename.endswith('.json'):
        file_path = os.path.join(directory, filename)
        with open(file_path, 'r') as f:
            data = json.load(f)

        #last_pop = [[0,0,0,0,0]]
        last_pop  = data['last_pop']
        S = data['S']
        G = data['G']
        mu = data['MU']
        epsilon = data['epsilon'] 
        dim = data['Dim']
        damp = data['Damp']
        seed = data['Seed']
        max_loops = data['max_loops']
        fail_loop = data['fail_loop']
        t = 0.5
        cag_result = None
        cag_size = None
        #print(my_custom_user_func(last_pop))
        if(fail_loop == max_loops):
            #Optimum was not found, check CAG
            #assert len(last_pop) != []
            max_graph_size = 100
            last_pop = Set_Node([tuple(i) for i in last_pop])
            CAG = community_assembly_graph.CAG(last_pop, my_custom_user_func, max_graph_size, False)
            print('S', S, 'G', G, 'Dim',dim, 'Damp',damp, 'MU', mu, 'Seed', int(seed))
            edge_df = pd.DataFrame(columns = ['from', 'to', 'edge_weight'])
            node_df = pd.DataFrame(columns = ['node', 'is_sink'])
            for key, value in CAG.items():
                new_row = {'node': str(key.members), 'is_sink': key.is_sink}
                node_df = pd.concat([node_df, pd.DataFrame(new_row, index=[0])], ignore_index=True)
                for k, v in value.items():
                    new_row = {'from': str(key.members), 'to': str(k.members), 'edge_weight': v}
                    edge_df = pd.concat([edge_df, pd.DataFrame(new_row, index=[0])], ignore_index=True)

            filename = '/'.join([directory, f"edge_S-{S}_G-{G}_Dim-{dim}_Damp-{damp}_MU-{mu}_Seed-{seed}_eps-{epsilon}.csv"])
            edge_df.to_csv(filename, index=False)
            filename = '/'.join([directory, f"node_S-{S}_G-{G}_Dim-{dim}_Damp-{damp}_MU-{mu}_Seed-{seed}_eps-{epsilon}.csv"])
            node_df.to_csv(filename, index=False)
            
            print("CAG = ",CAG)
            print("len(CAG) = ", len(CAG))

            if contains_optimum(CAG):
                cag_result = 'optimum is reachable'
                cag_size = len(CAG)
            else:
                if(len(CAG) < max_graph_size):
                    cag_result = "optimum is not reachable"
                    cag_size = len(CAG)

                elif(len(CAG) >= max_graph_size):
                    cag_result = "optimum isn't easily reachable"
                    cag_size = len(CAG)

        new_row = {'S': S, 'G':G, 'Dim':dim, 'Damp':damp, 'MU': mu, 'Seed': int(seed), 'max_loops': max_loops, 'epsilon': epsilon, 'fail_loop':fail_loop, 'cag_result':cag_result, 'cag_size':cag_size}
        df = pd.concat([df, pd.DataFrame(new_row, index=[0])], ignore_index=True)
filename = '/'.join([directory, "all_data.csv"])
df.to_csv(filename, index=False)