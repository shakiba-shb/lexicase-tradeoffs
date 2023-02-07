from leixcase_equation import *
import itertools
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import math
import collections
import random

dim = 5
L = 5
damp = 1
G = 500 #number of generations
S = 100 #population size 
p_thresh = .0001
# start = [(4,4,4,0,0,0), (0,0,0,4,4,4)]
start = [(4,0,0,0,0), (0,4,0,0,0), (0,0,4,0,0), (0,0,0,4,0), (0,0,0,0,4)]

# Adapted from https://stackoverflow.com/questions/9326851/generating-all-permutations-of-2-ones-and-3-zeroes-with-itertools
def combiner(zeros=3, ones=2, first_val = 0, second_val = 1):
    for indices in itertools.combinations(range(zeros+ones), ones):
        item = [first_val] * (zeros+ones)
        for index in indices:
            item[index] = second_val
        yield item


class community_info:
    def __init__(self, id, transitions_to):
        self.id = id
        self.transitions_to = transitions_to

    def __repr__(self) -> str:
        return str([all_genotypes[i] for i in self.id]) + "->" + str([all_genotypes[i] for i in self.transitions_to])



all_genotypes = list(set(list(itertools.product([i for i in range(L)], repeat=dim))))

# all_genotypes = list(combiner(dim,0)) + \
#                 list(combiner(dim-1,1)) + \
#                 list(combiner(dim-2,2)) + \
#                 list(combiner(dim-3,3)) + \
#                 list(combiner(dim-4,4)) + \
#                 list(combiner(dim-5,5)) + \
#                 list(combiner(dim-1,1,0,2)) + \
#                 list(combiner(dim-2,2,0,2)) + \
#                 list(set(itertools.permutations([0]*(dim-2)+[1,2]))) + \
#                 list(set(itertools.permutations([0]*(dim-3) + [1,1,2]))) + \
#                 list(set(itertools.permutations([0]*(dim-3)+ [1,2,2]))) 

# all_genotypes = [[0]*dim] + [[0]*i + [1] + [0]*(dim-i-1) for i in range(dim)] + [[0]*i + [2] + [0]*(dim-i-1) for i in range(dim)] + \
#                 [[1] + [0]*i + [1] + [0]*(dim-i-2) for i in range(dim-1)] + [[0]*i + [1] + [0]*(dim-i-2) + [1] for i in range(dim-1)] + \
#                 [[2] + [0]*i + [1] + [0]*(dim-i-2) for i in range(dim-1)] + [[0]*i + [1] + [0]*(dim-i-2) + [2] for i in range(dim-1)] + \
#                 [[1] + [0]*i + [2] + [0]*(dim-i-2) for i in range(dim-1)] + [[0]*i + [2] + [0]*(dim-i-2) + [1] for i in range(dim-1)] + \
#                 [[2] + [0]*i + [2] + [0]*(dim-i-2) for i in range(dim-1)]

# all_genotypes = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 1), (0, 0, 0, 1, 0), (0, 0, 0, 1, 1), (0, 0, 1, 0, 0), (0, 0, 1, 0, 1), (0, 0, 1, 1, 0), (0, 0, 1, 1, 1), 
#                  (0, 1, 0, 0, 0), (0, 1, 0, 0, 1), (0, 1, 0, 1, 0), (0, 1, 0, 1, 1), (0, 1, 1, 0, 0), (0, 1, 1, 0, 1), (0, 1, 1, 1, 0), (0, 1, 1, 1, 1), 
#                  (1, 0, 0, 0, 0), (1, 0, 0, 0, 1), (1, 0, 0, 1, 0), (1, 0, 0, 1, 1), (1, 0, 1, 0, 0), (1, 0, 1, 0, 1), (1, 0, 1, 1, 0), (1, 0, 1, 1, 1), 
#                  (1, 1, 0, 0, 0), (1, 1, 0, 0, 1), (1, 1, 0, 1, 0), (1, 1, 0, 1, 1), (1, 1, 1, 0, 0), (1, 1, 1, 0, 1), (1, 1, 1, 1, 0), (1, 1, 1, 1, 1)]

# all_genotypes = [[0,0,0,0,0], [1,0,0,0,0], [0,1,0,0,0], [0,0,1,0,0], [0,0,0,1,0], [0,0,0,0,1], 
#                  [0,0,0,2,0], [0,0,0,0,2], [2,0,0,0,0], [0,2,0,0,0], [0,0,2,0,0],
#                  [1,1,0,0,0], [1,0,0,1,0], [1,0,0,0,1], [0,1,0,1,0], [0,1,0,0,1],
#                  [1,0,1,0,0], [0,1,1,0,0], [0,0,1,1,0], [0,0,0,1,1], [0,0,1,0,1], 
#                  [2,1,0,0,0], [1,2,0,0,0], [1,0,0,0,2], [1,0,0,2,0], [1,0,2,0,0], 
#                  [2,0,0,0,1], [0,2,0,0,1], [0,0,0,2,1], [0,0,2,0,1], [0,0,0,1,2], 
#                  [2,2,0,0,0], [2,0,1,0,0], [1,1,1,0,0], [1,0,1,1,0], [1,0,1,0,1], [0,0,1,1,1]]
                 #[0,1,1,1,0],
                 #[0,0,1,1,2], [0,0,1,2,2], [0,0,2,1,2], [2,0,1,1,0]]#,
#                 [0,0,3,2,2], [0,0,3,1,2]]
# all_genotypes = list(itertools.product(*([[i for i in range(L)]]*dim)))

start = [all_genotypes.index(i) for i in start]

all_phenotypes = [fitness_function(i, damp) for i in all_genotypes]
# print([(i, all_genotypes[i]) for i in range(len(all_genotypes))])
print(start, [all_phenotypes[i] for i in start])
# exit()
# stable = set()
n_genotypes = len(all_genotypes)
# n_communities = 1 << n_genotypes
# print(n_communities)
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
    # print()
    # print("analyzing", i)
    indices = list(i)
    population = [all_phenotypes[j] for j in indices]
    # print(population)
    # exit()
    if population == []:
        all_communities[i] = community_info(i,frozenset([0]))
        return
    # print("seg fault?", population)
    prob = example.LexicaseFitness(population)
    # print("no seg fault")
    # print(/"")
    P_survival = (np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G
    survivors = [indices[j] for j in range(len(indices)) if P_survival[j] >= p_thresh]
    # print()
    # print([all_genotypes[j] for j in indices], P_survival)
    # print(survivors)
    while survivors != indices:
        # print("loopign", indices, survivors)
        indices = survivors
        population = [all_phenotypes[j] for j in indices]
        if population == []:
            all_communities[i] = community_info(i,frozenset([0]))
            return
        # print("seg fault?", population)
        prob = example.LexicaseFitness(population)
        # print("no seg fault")
        # print(/"")
        P_survival = (np.ones(len(prob)) - (np.ones(len(prob)) - prob)**S)**G
        survivors = [indices[j] for j in range(len(indices)) if P_survival[j] >= p_thresh]
        
    # id = pos_to_dec(survivors)
    # print(i, prob, P_survival, survivors, p_thresh)
    return community_info(i, frozenset(survivors))
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

unexplored = collections.OrderedDict()
unexplored[frozenset(start)] = True
res = analyze_community(frozenset(start))
all_communities[frozenset(start)] = res
graph.add_node(frozenset(start), pos = (sum([i for j in [all_genotypes[i] for i in start] for i in j ]) + (random.random()-.5),len(start)+ (random.random()-.5)))

# print(adj)

count = 0
while unexplored and count < 500:
    count += 1
    print(len(unexplored))
    curr = unexplored.popitem(False)[0]

    # assert(graph.has_node(curr))

    adj_genotypes = set()
    for i in curr:
        # print(i, adj[i])
        adj_genotypes.update(adj[i])
    
    # print(curr, all_communities[curr])
    assert(all_communities[curr].transitions_to == curr)
    for a in adj_genotypes:

        if a not in curr:
            new_comm = curr.union(frozenset([a]))
            # print(curr, a, new_comm)
            if new_comm not in all_communities:
                comm_info = analyze_community(new_comm)
                if comm_info.transitions_to not in all_communities:
                    unexplored[comm_info.transitions_to] = True
                    if not graph.has_node(comm_info.transitions_to):
                        graph.add_node(comm_info.transitions_to, pos = (sum([i for j in [all_genotypes[i] for i in comm_info.transitions_to] for i in j ]) + (random.random()-.5), len(comm_info.transitions_to) + (random.random()-.5)))
                all_communities[new_comm] = comm_info
                if comm_info.transitions_to != new_comm:
                    trans_info = analyze_community(all_communities[new_comm].transitions_to)
                    all_communities[comm_info.transitions_to] = trans_info
                # print(new_comm, all_communities[new_comm].transitions_to)
                # assert(all_communities[new_comm].transitions_to == all_communities[all_communities[new_comm].transitions_to].transitions_to)
                # print("new", all_communities[new_comm].transitions_to)

            result = all_communities[new_comm].transitions_to
            if not graph.has_node(result):
                graph.add_node(result, pos = (sum([i for j in [all_genotypes[i] for i in result] for i in j ]) + (random.random()-.5),len(result) + (random.random()-.5)))                
            if result != curr:
                graph.add_edge(curr, result, l=a)

print(all_communities.values())

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

o_degrees = graph.out_degree() # Dict with Node ID, Degree
i_degrees = graph.in_degree() # Dict with Node ID, Degree
nodes = graph.nodes()
graph.remove_nodes_from([n for n in nodes if o_degrees[n] == 0])
nodes = graph.nodes()

# n_color = np.asarray([2*(o_degrees[n] == 0) + 1*(i_degrees[n] == 0) for n in nodes])

# print(nodes, n_color)
# print(all_genotypes)
print(set([tuple(all_genotypes[i]) for n in nodes for i in n if o_degrees[n] == 0]))  
print([[all_genotypes[i] for i in n] for n in nodes if o_degrees[n] == 0])  

fig, ax = plt.subplots(figsize=(20, 20), dpi=250)

nx.relabel_nodes(graph, {n:str([i for i in n]) for n in nodes}, copy=False)
pr = nx.pagerank(graph, alpha=0.9)
max_pr = max(pr.values())

pos=nx.get_node_attributes(graph,'pos')
# print(pos)
for key in pos.keys():
    new_x = pos[key][0]
    new_y = math.log(pr[key]/max_pr)
    pos[key] = {"pos":(new_x, new_y)}
nx.set_node_attributes(graph, pos)
pos=nx.get_node_attributes(graph,'pos')
nodes = graph.nodes()

# print(pr)

n_color = np.asarray([max([all_genotypes[i].count(L-1) for i in eval(n) ]) for n in nodes])
# pos = nx.spring_layout(graph)
drawn_nodes = nx.draw_networkx_nodes(graph,pos, node_color=n_color, ax=ax)
nx.draw_networkx_edges(graph, pos, ax=ax)
# print(drawn_nodes)
# cmap = plt.cm.plasma
# pc = mpl.collections.PatchCollection(drawn_nodes, cmap=cmap)
# pc.set_array(list(range(max(n_color))))

ax.set_axis_off()
plt.colorbar(drawn_nodes)
plt.savefig("test.png")
