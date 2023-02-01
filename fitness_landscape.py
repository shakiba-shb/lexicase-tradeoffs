from leixcase_equation import *
import itertools
import networkx as nx
import matplotlib.pyplot as plt
import math

dim = 3
L = 5
damp = 2
G = 1 #number of generations
S = 10 #population size 
p_thresh = .0001

class community_info:
    def __init__(self, id, transitions_to):
        self.id = id
        self.transitions_to = transitions_to


all_genotypes = list(itertools.product(*([[i for i in range(L)]]*dim)))
all_phenotypes = [fitness_function(i, damp) for i in all_genotypes]
print(all_genotypes)
# stable = set()
n_genotypes = len(all_genotypes)
n_communities = 1 << n_genotypes
print(n_communities)
all_communities = {}

# adj = [[0 for i in range(n_genotypes)] for j in range(n_genotypes)]
adj = {i: [] for i in range(n_genotypes)}
for i in range(n_genotypes):
    for j in range(i+1, n_genotypes):
        if math.dist(all_genotypes[i], all_genotypes[j]) == 1:
            adj[i].append(j)
            adj[j].append(i)


def pos_to_dec(pos):
    result = 0
    for val in pos:
        result += 1<<val
    return result


def dec_to_pos(dec):
    result = []
    for i in range(n_genotypes):
        if dec & (1<<i):
            result.append(i)

    return result


def analyze_community(i):
    indices = [j for j in range(n_genotypes) if i & (1<<j)]
    population = [all_phenotypes[j] for j in indices]
    if population == []:
        all_communities[i] = community_info(i,i)
        return
    # print("seg fault?", population)
    prob = example.LexicaseFitness(population)
    # print("no seg fault")
    # print(/"")
    P_survival = (np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G
    survivors = [indices[j] for j in range(len(indices)) if P_survival[j] >= p_thresh]
    while survivors != indices:
        # print("loopign", indices, survivors)
        indices = survivors
        population = [all_phenotypes[j] for j in indices]
        if population == []:
            all_communities[i] = community_info(i,i)
            return
        # print("seg fault?", population)
        prob = example.LexicaseFitness(population)
        # print("no seg fault")
        # print(/"")
        P_survival = (np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G
        survivors = [indices[j] for j in range(len(indices)) if P_survival[j] >= p_thresh]
        
    id = pos_to_dec(survivors)
    # print(bin(i), bin(id), prob, P_survival, survivors, p_thresh)
    all_communities[i] = community_info(i, id)
    # if i == id:
    #     stable.add(i)

# for i in range(1, n_communities):
#     indices = [j for j in range(n_genotypes) if i & (1<<j)]
#     population = [all_phenotypes[j] for j in indices]
#     prob = example.LexicaseFitness(population)
#     P_survival = (np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G

#     survivors = [indices[j] for j in range(len(indices)) if P_survival[j] >= p_thresh]
#     id = pos_to_dec(survivors)
#     all_communities[i] = community_info(i, id)
#     if i == id:
#         stable.add(i)


graph = nx.DiGraph()
# G.add_nodes_from(stable)
# print(stable)

unexplored = set([1])
analyze_community(1)
print(adj)

# count = 0
while unexplored:
    # count += 1
    print(len(unexplored))
    curr = unexplored.pop()

    adj_genotypes = set()
    for i in dec_to_pos(curr):
        adj_genotypes.update(adj[i])
            
    for a in adj_genotypes:
        if not curr & (1 << a):
            new_comm = curr + (1<<a)
            # print(curr, i, a, new_comm)
            if new_comm not in all_communities:
                analyze_community(new_comm)
                analyze_community(all_communities[new_comm].transitions_to)
                # assert(all_communities[new_comm].transitions_to == all_communities[all_communities[new_comm].transitions_to].transitions_to)
                # print("new", all_communities[new_comm].transitions_to)
                unexplored.add(all_communities[new_comm].transitions_to)

            result = all_communities[new_comm].transitions_to
            if result != curr:
                graph.add_edge(curr, result, l=a) 

# for comm in stable:
#     adj_comm = set()
#     for i in dec_to_pos(comm):
#         adj_comm.update(adj[i])
#     adj_comm -= set(dec_to_pos(comm))
#     # genotypes = [all_genotypes[i] for i in dec_to_pos(comm)]
#     for i in adj_comm:
#         # if not comm & (1 << i):
#             # new_genotype = all_genotypes[comm + (1<<i)]

#         result = all_communities[comm + (1<<i)].transitions_to
#         if result != comm:
#             G.add_edge(comm, result, label=i) 

o_degrees = graph.out_degree() #Dict with Node ID, Degree
i_degrees = graph.in_degree() #Dict with Node ID, Degree
nodes = graph.nodes()
n_color = np.asarray([2*(o_degrees[n] == 0) + 1*(i_degrees[n] == 0) for n in nodes])
# print(nodes, n_color)
print([n for n in nodes if o_degrees[n] == 0])  
# pos = nx.spring_layout(graph)
# nx.draw(graph, pos, with_labels=False, font_weight='bold', node_color=n_color)
# nx.draw_networkx_edge_labels(graph, pos)
# plt.savefig("test.png")
