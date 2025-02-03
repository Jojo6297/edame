import numpy as np
import pandas as pd
import networkx as nx
import random
from collections import Counter
from itertools import *
import seaborn as sns
import copy

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
#rcParams
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['figure.dpi'] = 300

import sys
print(f'\n_____________________________________________________________________________________________________________________________________________________________________\n')

# FUNCTIONS #
def operation1(mat1, mat2): #mat2 is a list
    '''to generate attractors'''
    o_count = 0
    new_list = []
    while True:
        simulated_matrix = [None]*len(mat2)
        mat = list(np.dot(mat1, mat2))
        mat2 = list(mat2)
        new_list.append(mat2)
        for i in range(len(mat)):
            if mat[i] > 0: #  and mat2[i] > 0:
                simulated_matrix[i] = 1
            elif mat[i] == 0:
                simulated_matrix[i] = mat2[i]
            else:
                simulated_matrix[i] = 0
        o_count += 1
        if mat2 == simulated_matrix:                                                                                                                    # for asymptotic attractor
            break                    
        elif simulated_matrix in new_list:                                                                                                              # for cycles
            #print ("Cycle Exists")
            break                          
        mat2 = simulated_matrix
        #print(o_count)
    return (simulated_matrix)

def entropy_temp(vec):
    t1 = len(vec)
    ones_count = Counter(vec)[1]
    p1 = ones_count/float(t1)
    p0 = 1 - p1
    return (p1)

#defining "AND" operation
def op(pair):
    pair = list(pair)
    return(pair[0]*pair[1])

# This function accepts a matrix, a function and how many times you want to repeat the operations
def shift_score(twovec, op, rep = 5):
    x1 = np.array(list(map(op, twovec)))
    x2 = entropy_temp(x1)
    p3 = entropy_temp(np.transpose(twovec)[0])
    p4 = entropy_temp(np.transpose(twovec)[1])
    x4 = (p3*p4)
    x5 = np.sqrt(((p3*p4)*(1-(p3*p4)))/(len(twovec)))
    if (x5 == 0):
        return 0
    return ((x2 - x4) / x5)

# create network from positive, negative edges
def gen_network(node_num, pos_edges, neg_edges):
    Adj_upper = pd.DataFrame(np.zeros((node_num, node_num)))
    for i in pos_edges:
        Adj_upper.loc[i[0]][i[1]] = 1
    for j in neg_edges:
        Adj_upper.loc[j[0]][j[1]] = -1
    return (Adj_upper + np.transpose(Adj_upper))                                                                                                        # Adjacency of the ESABO network #

# Generating unique attractors for a network
def gen_attr(sample_num, df_phyla, netw):
    final_matrices = []
    for i in range(sample_num):
        a2 = df_phyla[df_phyla.columns[i]]
        M2 = operation1(netw, a2)
        final_matrices.append(M2)
    uniq_attractor = [list(x) for x in set(tuple(i) for i in final_matrices)]
    uniq_attractor = np.array(uniq_attractor)
    return (uniq_attractor)

# Finding Jaccard index
def Jaccard(li1, li2):
    '''takes lists of tuples'''
    intersection = set.intersection(set(li1), set(li2))
    union = set.union(set(li1), set(li2))
    return(len(intersection)/len(union))

def mod2(pair):
    pair = list(pair)
    if pair[0] == pair[1]:
        return 0
    else:
        return 1

def hamming(a, b):
    '''takes individual attr'''
    return (sum([mod2([a[i], b[i]]) for i in range(len(a))]))


def attr_difference(setA, setB):
    '''For each attractor in setA find its closest attractors in SetB (by hamming)'''
    hamm_dist_dict = {}
    for item1 in setA:
        hamm_dist = []
        for item2 in setB:
            hamm_dist.append(hamming(item1, item2))
        hamm_dist_dict[item1] = [setB[i] for i in range(len(hamm_dist)) if hamm_dist[i] == min(hamm_dist) and hamm_dist[i] != 0]

    set_A_B = []
    for key,val in hamm_dist_dict.items():
        set_A_B.append([tuple([k-v1 for k,v1 in zip(key,v)]) for v in val])                                                                             # Find (cD - A*)

    return(set_A_B)


############## Compare the original and obtained attractors to find the switching events ##############
def sorted_nodes_array(attr, attr_astr, species_num):
    cD = list(set(attr) - set(attr_astr))                                                                                                               # Find destroyed attractors #
    cC = list(set(attr_astr) - set(attr))                                                                                                               # Find created attractors #
    a_astr = copy.deepcopy(attr_astr)                                                                                                                                  # New network attractors (A*) #


    # For each attractor in cD find its closest attractors in A* (by hamming)
    cD_a_astr = attr_difference(cD, a_astr)                                                                                                             # Find (cD - A*)
    # For each attractor in cC find its closest attractors in A (by hamming)
    a_cC = attr_difference(cC, attr)                                                                                                                    # Find (cC - A)
    # For each attractor in cD find its closest attractors in cC (by hamming)
    cD_cC = attr_difference(cD, cC)                                                                                                                     # Find (cD - cC)


    # For each matrix, stack the attractors together and generate an array with absolute sums over the rows.
    a_vec = pd.DataFrame([list(item) for sublist in cD_a_astr for item in sublist]).abs().sum(axis = 0).tolist()
    b_vec = pd.DataFrame([list(item) for sublist in a_cC for item in sublist]).abs().sum(axis = 0).tolist()
    c_vec = pd.DataFrame([list(item) for sublist in cD_cC for item in sublist]).abs().sum(axis = 0).tolist()
    new_list_vec = [item for item in [a_vec, b_vec, c_vec] if len(item) > 0]                                                                            # Keep non zero vectors

    # Sum the three arrays and create the final "SORTED NODES ARRAY (S_{N})"
    plist = pd.DataFrame(new_list_vec).abs().sum(axis = 0).tolist()
    switch_linedict = dict(zip(range(species_num), plist))
    return_dict = dict(sorted(switch_linedict.items(), key = lambda x:x[1], reverse = True))
    return (return_dict)


def netw_to_edges(g):
    temp_upper = pd.DataFrame(np.triu(g, k = 0))
    mask1 = (temp_upper == 1)
    mask2 = (temp_upper == -1)
    stacked1 = mask1.stack()
    stacked2 = mask2.stack()
    posE = stacked1[stacked1].index.tolist()                                                                                                            # get list of positive edges
    negE = stacked2[stacked2].index.tolist()                                                                                                            # get list of negative edges
    return(posE, negE)

# Possible edge combinations
P_c = list(set([i for i in permutations([0,0,1,1,-1,-1], 2)]))

seed1 = 921

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

np.random.seed(seed1)
no_of_species = 10
tot_edges = float(20)
m = int(tot_edges//2)

print(f'\n\nWe start with an original network (G) with {int(no_of_species)} nodes and {int(tot_edges)} edges.')

phyla_df = pd.DataFrame(product([0,1],repeat=no_of_species), columns= [i for i in range(no_of_species)]).T
no_of_samples = phyla_df.shape[1]

############## Generating a random network to be used as an input for inference algorithm ##############
while True:
    # generate a random initial connected network
    G = nx.gnm_random_graph(no_of_species, tot_edges, seed = seed1, directed=False)
    if nx.is_connected(G) == True:
        print(f'\nG is a connected network')
        break

Adj = nx.to_numpy_matrix(G)
Adj_upper = pd.DataFrame(np.triu(Adj, k = 0))
mask = (Adj_upper == 1)
stacked = mask.stack()
index_list = stacked[stacked].index.tolist()  # get list of positive edges

random.seed(seed1)
neg_index = sorted(random.sample(index_list, k = m))
for i,j in neg_index:
    Adj_upper.loc[i, j] = -1                                                                                                                            # Randomly assigning -1 to tot_edges/2 edges #
Int_mat = Adj_upper + np.transpose(Adj_upper)                                                                                                           # Int_mat: final matrix with m/2 +ve and -ve int. #
pos_index = sorted([item for item in index_list if item not in neg_index])

############## Find  attractors of the network based on thte majority vore update rules (refer to equation 1 in manuscript) ##############
uniq_attractor = gen_attr(no_of_samples, phyla_df, Int_mat)
A = [tuple(j) for j in uniq_attractor]                                                                                                                  # Attractors of the original network #

##############  Attempt to infer the initial network using ESABO method ##############
scores = []                                                                                                                                             # To save shift scores #
Edge_list = []
score_sign_list = []

for i in range(no_of_species):
    for j in range(no_of_species):
        if j > i:
            Edge_list.append((i ,j))
            scores.append(shift_score(uniq_attractor[:, [i,j]], op))                                                                                    # Calculating the shift scores (z-scores) #

score_sign_list = [[i,j] for i,j in zip(Edge_list, scores)]
score_sign_list = sorted(score_sign_list, key=lambda x: x[1], reverse = True)

m_positive_zscores_edges = [i[0] for i in score_sign_list[:m]]                                                                                          # Select edges with highest positive scores ## could be done by choosing a positive threshold, we chose top m edges for simplicity #
m_negative_zscores_edges = [i[0] for i in score_sign_list[-m:]]                                                                                         # Select edges with highest negative scores ## could be done by choosing a negative threshold, we chose bottom m edges for simplicity #

Rec_pos_J = Jaccard(pos_index, m_positive_zscores_edges)                                                                                                # Find overlap between original positive edges and infered positive edges #
Rec_neg_J = Jaccard(neg_index, m_negative_zscores_edges)                                                                                                # Find overlap between original negative edges and infered negative edges #
Rec_J_index = (Rec_pos_J + Rec_neg_J) / 2.0                                                                                                             # Find overlap between original and infered edges #
print(f'\nFirst we find attractors of the original network and then find the ESABO network using these attractors.\n')

print(f"Jaccard index between Original network and Network generated by ESABO = {Rec_J_index}")

# Save the ESABO network
Rec_Int_mat2 = gen_network(no_of_species, m_positive_zscores_edges, m_negative_zscores_edges)           # Adjacency of the ESABO network #

# Generating unique attractors for ESABO implementation
Rec_uniq_attractor2 = gen_attr(no_of_samples, phyla_df, Rec_Int_mat2)
Rec_li_uniq_attr2 = [tuple(j) for j in Rec_uniq_attractor2]                                                                                             # Attractors of the ESABO network #

J_Attr = Jaccard(A, Rec_li_uniq_attr2)                                                                                                                  # Find overlap between original and infered edges #
print(f"Jaccard index between Attractors of original network and Attractors of ESABO network = {J_Attr}\n")

############## Compare the original and ESABO attractors to find the switching events ##############
############## Obtain the sorted nodes array ##############
EdgeChange_dict = sorted_nodes_array(A, Rec_li_uniq_attr2, no_of_species)
EdgeChange_dict_val = copy.deepcopy(EdgeChange_dict)


print(f'Then obtain the sorted nodes array:')
print(f"Nodes (keys) with their respective peak heights (values):\n{EdgeChange_dict_val}")
print(f"Sorted nodes array:\n{list(EdgeChange_dict_val.keys())}")
print(f"\n")

pos_new_original = [item for item in m_positive_zscores_edges if item not in pos_index]
pos_original_new = [item for item in pos_index if item not in m_positive_zscores_edges]
neg_new_original = [item for item in m_negative_zscores_edges if item not in neg_index]
neg_original_new = [item for item in neg_index if item not in m_negative_zscores_edges]

print(f'Edge differences between the original and ESABO network:')
print(f'Original positive edges not in ESABO network -> {pos_original_new}')
print(f'Original negative edges not in ESABO network -> {neg_original_new}')
print(f'ESABO positive edges not in original network -> {pos_new_original}')
print(f'ESABO negative edges not in original network -> {neg_new_original}\n\n')

# Save the networks and attractors from original and ESABO implementation
original_netw = [G, pos_index, neg_index, A]
esabo_netw = [m_positive_zscores_edges, m_negative_zscores_edges, Rec_li_uniq_attr2, Rec_J_index, J_Attr]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################### Creating indices of nodes to iterate over ######################
cont_counting = [[j for j in range(i)] for i in range(2,no_of_species)]
allnodes = []
for item in cont_counting:
    allnodes.append([item+[p] for p in range(len(item),no_of_species)])

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(f"\'\'\'RUNNING THE EDAME ALGORITHM\'\'\' ")
pic_length = 10                                                                                                                                         # We keep looking for a solution which takes less than 7 iterations, to plot the figure afterwards.
while pic_length > 8:
    ## Initialize values and list objects
    random.seed(None)
    init_n = 3                                                                                                                                          # Minimum number of nodes considered at once
    n = init_n                                                                                                                                          # Counter for number of nodes
    JI_t = 0                                                                                                                                            # Initialize temporary Jaccard index => 0
    iter_count = 0                                                                                                                                      # Counter for number of iterations
    fail = 0                                                                                                                                            # Counter for failures
    p_count = 0                                                                                                                                         # Counter for printing

    temp_G = Rec_Int_mat2.copy(deep = True)                                                                                                             # To update new networks with JI > previous JI
    JI_A = Jaccard(Rec_li_uniq_attr2, A)                                                                                                                # To update attractor JI

    network_info = [esabo_netw]                                                                                                                         # To save all networks

    # In case algorithm fails, save some variables to use again
    J_val = JI_A
    EdgeChange_dict_val = copy.deepcopy(EdgeChange_dict)
    tempfail = 0                                                                                                                                        # Counter to trigger the adjustment of n when the algorithm fails
    faillimit = 5                                                                                                                                       # Threshold to stop the algorithm if it fails a lot of times

    print(f"Initial Sorted nodes array: {list(EdgeChange_dict_val.keys())}")                                                                            # Sorted nodes array

    ## Run the loop to iteratively update the ESABO network by comparing attractors of the temporary and the new network generated by changing edges, 
    # that are indicated by the sorted nodes array generated after every update
    while n >= init_n  and n <= no_of_species and fail < faillimit:                                                                                     # Loop should run until 3 <= n <= total number of nodes, and before the faillimit is achieved
        node_count = len(allnodes[n-3][0])                                                                                                              # number of nodes to be considered
        print(f'\nNumber of nodes considered: {node_count}')
        
        # Iterate over items in allnodes[0], then in allnodes[1], and so on
        for item__1 in allnodes[n-3]:                                                          
            if len(list(EdgeChange_dict.keys())) > 0:                                                                                                   # If sorted nodes array exists
                check_nodes = list(map(list(EdgeChange_dict.keys()).__getitem__, item__1))                                                              # Create a list of nodes from sorted node array present at the indices in item__1
            else:
                check_nodes = list(EdgeChange_dict.keys())[:n]

            check_edge = [i for i in combinations(check_nodes, 2)]                                                                                      # Combinations of nodes to create a list of possible edges to be modified
            edge_comb = [i for i in combinations(check_edge, 2)]                                                                                        # Combinations of edges to create pairs of edges to work with

            #print(f"Position of the last node = {node_count}")
            print(f"N_list: {check_nodes}")
            if p_count == 0:
                print(f'\nE_list and E_pair printed once for explaination: ')
                print(f"E_list: {check_edge}")
                print(f"E_pair: {edge_comb}\n")
                p_count += 1

            if fail > 0:
                edge_comb = random.sample(edge_comb, len(edge_comb))                                                                                    # If the algorithm fails once, edge pairs are selected at random

            for item in edge_comb:
                number = 0
                count = 0
                P_c = random.sample(P_c, len(P_c))                                                                                                      # Randomly picking up possible edge value combinations

                for comb in P_c:
                    G_t = temp_G.copy(deep = True)                                                                                                      # store a temporary network

                    # check for each comb in P_c, if any edge pair in "E_pair" has values == P_c[0] and P_c[1]
                    if G_t[item[0][0]][item[0][1]] == comb[0] and G_t[item[1][0]][item[1][1]] == comb[1]:
                        #number += 1            

                        # If E_pair fits a P_c value -> Run a loop on P_c and give this edge pair new values
                        for diffedge in [o for o in P_c if o != comb]:
                            G_t[item[0][0]][item[0][1]] = diffedge[0]
                            G_t[item[0][1]][item[0][0]] = diffedge[0]
                            G_t[item[1][0]][item[1][1]] = diffedge[1]
                            G_t[item[1][1]][item[1][0]] = diffedge[1]

                            # with the modified edges in the network, create new attractors
                            temp_uniq_attractor = gen_attr(no_of_samples, phyla_df, G_t)
                            A_t = [tuple(j) for j in temp_uniq_attractor]

                            JI_t = Jaccard(A_t, A)                                                                                                      # Jaccard index between modified network attractors and previous network attractors
                            
                            if number > 5:                                                                                                              # If a single edge pair can't be modified to any comb in P_c after 5 tries, try the next edge pair
                                break

                            if JI_t < JI_A:                                                                                                             # If new JI < previous attractor JI (JI_A)
                                continue

                            
                            if JI_t > JI_A and JI_t != 1:                                                                                               # If (new JI > previous attractor JI) but (new JI != 1)
                                iter_count += 1
                                print(f'Iteration --> {iter_count}')
                                print(f"Attractor Jaccard index increase -> {JI_A} -> {JI_t}")
                                #print(f'JI_prev = {JI_A}, JI_new = {JI_t}')
                                print(f'Edge pair modified: {item}')

                                # Save new network information for the figure later on
                                iter_posedge, iter_negedge = netw_to_edges(G_t)
                                pos_Jac = Jaccard([tuple(sorted(i)) for i in pos_index], [tuple(sorted(i)) for i in iter_posedge])
                                neg_Jac = Jaccard([tuple(sorted(i)) for i in neg_index], [tuple(sorted(i)) for i in iter_negedge])
                                Jac = (pos_Jac + neg_Jac) / 2.0

                                network_info.append([iter_posedge, iter_negedge, A_t, Jac, JI_t])
                                temp_G = G_t.copy(deep = True)                                                                                          # Updated network

                                ############## Compare the original and updated network attractors to find the switching events ##############
                                ############## Obtain sorted nodes array ##############
                                EdgeChange_dict = sorted_nodes_array(A, A_t, no_of_species)
                                print(f"Updated sorted nodes array:\n{list(EdgeChange_dict.keys())}")
                                break

                            if JI_t == 1:
                                iter_count += 1
                                print(f'Iteration --> {iter_count}')
                                print(f"Attractor Jaccard index increase -> {JI_A} -> {JI_t}")
                                #print(f'JI_prev = {JI_A}, JI_new = {JI_t}')
                                print(f'Edge pair modified: {item}')
                                print(f'Network detected\n')

                                # Save new network information for the figure later on
                                iter_posedge, iter_negedge = netw_to_edges(G_t)
                                pos_Jac = Jaccard([tuple(sorted(i)) for i in pos_index], [tuple(sorted(i)) for i in iter_posedge])
                                neg_Jac = Jaccard([tuple(sorted(i)) for i in neg_index], [tuple(sorted(i)) for i in iter_negedge])
                                Jac = (pos_Jac + neg_Jac) / 2.0
                                network_info.append([iter_posedge, iter_negedge, A_t, Jac, JI_t])

                                if iter_count > 8:
                                    print(f"\nThe algorithm found the original network, but took > 8 iterations, therefore, we run it again to get the result in <= 8 iterations to plot a clean figure. It might take some time.\n")
                                    pic_length = iter_count

                                    JI_A = J_val
                                    JI_t = 0
                                    iter_count = 0
                                    EdgeChange_dict = copy.deepcopy(EdgeChange_dict_val)
                                    temp_G = Rec_Int_mat2.copy(deep = True)
                                    network_info = [esabo_netw]

                                else:
                                    pic_length = iter_count

                                break

                        # Break or continue the other loops based on the comparison of previous and new JI.
                        if JI_t > JI_A and JI_t != 1:
                            break
                        elif JI_t == 1:
                            break
                        else:
                            continue

                    if number == 0:
                        count += 1
                    if count > init_n:
                        break
                    
                # Break or continue the other loops based on the comparison of previous and new JI.
                if JI_t > JI_A and JI_t != 1:
                    break
                elif JI_t == 1:
                    break
                else:
                    continue

            # Break or continue the other loops based on the comparison of previous and new JI.
            if JI_t > JI_A and JI_t != 1:
                break
            elif JI_t == 1:
                break
            else:
                node_count += 1                                                                                                                         # If the network was not updated, increment node_count
                continue
        
        # Break or continue the other loops based on the comparison of previous and new JI.
        if JI_t > JI_A and JI_t != 1:
            n = init_n                                                                                                                                  # If the network was updated, start the node counter from the start (n_init)
            JI_A = JI_t
            
        elif JI_t == 1:
            break

        # If we have exhausted all the nodes, the algorithm fails. We have to start again. Due to the randomness in edge pair selection and P_c selection, the trajectory of edge modificstions can change.
        elif len(allnodes[n-3][0]) == no_of_species:
            fail += 1
            tempfail += 1

            if fail == 1:
                print(f'\nAlgorithm failed once')
            else:
                print(f'\nAlgorithm failed -> {fail} times')

            if tempfail + init_n <= no_of_species:                                                                                                      # Tempfail is a counter for failure, if algorithm fails, start the algorithm with n += 1
                n = tempfail + init_n
            else:
                tempfail = 2
                n = tempfail + init_n                                                                                                                   # After (tempfail + init_n) > no_of_species, reset tempfail = 2
            print(f'Starting again with n --> {n}')

            # Setting all the variables to initial values to retart the process
            JI_A = J_val
            JI_t = 0
            iter_count = 0
            EdgeChange_dict = copy.deepcopy(EdgeChange_dict_val)
            temp_G = Rec_Int_mat2.copy(deep = True)
            network_info = [esabo_netw]
            continue
        else:
            n += 1                                                                                                                                      # If (new JI < previous attractor JI) increment n
            continue

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

######################### Let us plot the figure to visualize the network evolution ######################
all_pos_edge = [original_netw[1]] + [i[0] for i in network_info]
all_neg_edge = [original_netw[2]] + [i[1] for i in network_info]
all_attr = [original_netw[3]] + [i[2] for i in network_info]
all_Jindex = [i[3] for i in network_info]
all_Jattr = [i[4] for i in network_info]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

len_A = [len(i) for i in all_attr]
same_attrlen_A = [len(i) for i in all_attr]
#print(f'All attractors ---> {len_A}')
same_attr_as_original = [sorted(list(set.intersection(set(all_attr[i]), set(all_attr[0])))) for i in range(len(all_attr))]
len_S = [len(i) for i in same_attr_as_original]
#print(f'Attr present in both sets ---> {len_S}')
different_attr_than_original  = [sorted(list((Counter(all_attr[i]) - Counter(all_attr[0])).elements())) for i in range(len(all_attr))]
len_D = [len(i) for i in different_attr_than_original]
#print(f'Attr present in this network but not in original ---> {len_D}')
different_attr_than_individuals  = [sorted(list((Counter(all_attr[0]) - Counter(all_attr[i])).elements())) for i in range(len(all_attr))]
len_D2 = [len(i) for i in different_attr_than_individuals]
#print(f'Attractors present in original network but not in this ---> {len_D2}')
Attr_combined = [pd.DataFrame(same_attr_as_original[i]+different_attr_than_original[i]) for i in range(len(all_attr))]
same_attr = copy.deepcopy(Attr_combined)
diff_attr = copy.deepcopy(Attr_combined)

for j in range(len(same_attr)):
    same_attr[j].loc[same_attr[j].index < len_S[j]] = [False for i in range(no_of_species)]
    same_attr[j].loc[same_attr[j].index >= len_S[j]] = [True for i in range(no_of_species)]
    same_attr[j] = same_attr[j].to_numpy()

    diff_attr[j].loc[diff_attr[j].index < len_S[j]] = [True for i in range(no_of_species)]
    diff_attr[j].loc[diff_attr[j].index >= len_S[j]] = [False for i in range(no_of_species)]
    diff_attr[j] = diff_attr[j].to_numpy()

max_miss = max([len(i) for i in different_attr_than_individuals])
len_M = [len(item) for item in different_attr_than_individuals]
missing_attr = [pd.DataFrame(item + (max_miss - len(item))*[tuple(1 for i in range(no_of_species))]) for item in different_attr_than_individuals]

init = copy.deepcopy(missing_attr)
ones_s = copy.deepcopy(missing_attr)
for j in range(len(init)):
    init[j].loc[init[j].index < len_D2[j]] = [False for i in range(no_of_species)]
    init[j].loc[init[j].index >= len_D2[j]] = [True for i in range(no_of_species)]
    init[j] = init[j].to_numpy()

    ones_s[j].loc[ones_s[j].index < len_D2[j]] = [True for i in range(no_of_species)]
    ones_s[j].loc[ones_s[j].index >= len_D2[j]] = [False for i in range(no_of_species)]
    ones_s[j] = ones_s[j].to_numpy()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Defining color maps
Greens_map = sns.light_palette((125, 90, 55), input="husl", n_colors= 10)
Blacks_r_map = sns.light_palette((  0,  0,  0), input="husl", n_colors= 10, reverse = True)
Blacks_map = sns.light_palette((  0,  0,  0), input="husl", n_colors= 10, reverse = False)
Greys_r_map  = sns.light_palette((249,  0, 50), input="husl", n_colors= 10, reverse = True)
Greys_map  = sns.light_palette((249,  0, 50), input="husl", n_colors= 10, reverse = False)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## The final figure

print(f"Plotting the figure!!")
fig = plt.figure(figsize = (len(network_info)+3, 6))

ax__ = plt.subplot2grid((23,len(all_pos_edge)), (0,0), rowspan=11, colspan= 1)
nodes = mlines.Line2D([], [], color='k', linestyle = '', linewidth = 1,
                         marker='o', markersize=3, 
                        label = f'Node')
blue_line = mlines.Line2D([], [], color= 'blue', linestyle = '-', linewidth = 1,
                        label='Positive edges')
red_line = mlines.Line2D([], [], color= 'red', linestyle = '-', linewidth = 1,
                        label='Negative edges')
black_line = mlines.Line2D([], [], color= Blacks_r_map[0], linestyle = '--', linewidth = 1,
                        label = f'Edges in this\nnetwork, but\n not in ' + r'$G_0$')
silver_line = mlines.Line2D([], [], color=Greys_r_map[0], linestyle = '--', linewidth = 1,
                        label = f'Edges in ' + r'$G_0$'+ f',\nbut not in\nthis network')
lg = ax__.legend(handles=[nodes, blue_line, red_line, black_line, silver_line], fontsize = 6, loc = 'upper right', labelspacing = .4, frameon = False) #, bbox_to_anchor=(.01,.3, .2, .6))
lg.set_title(f'Network Properties',prop={'size':6.5, 'weight':'bold'})
ax_ll = plt.gca().add_artist(lg)

green_patch = mpatches.Patch(color = Greens_map[9], label = f'Attractors of the\nnetwork')
black_patch = mpatches.Patch(color = Blacks_r_map[0], label = f'Attractors of this\nnetwork,\nbut not in ' + r'$A_0$')
silver_patch = mpatches.Patch(color = Greys_r_map[0], label = f'Attractors in\n' + r'$A_0$'+ f', but not\nof this network')
lg2 = ax__.legend(handles=[green_patch, black_patch, silver_patch], fontsize = 6, loc = 'lower right', labelspacing = .4, frameon = False) #, bbox_to_anchor=(.01,.3, .2, .6))
lg2.set_title(f'Attractor Properties',prop={'size':6.5, 'weight':'bold'})
plt.axis('off')

ax_ = plt.subplot2grid((23,len(all_pos_edge)), (0,1), rowspan=10, colspan= len(all_pos_edge)-1)
plt.axis([-.5, len(all_pos_edge)-1.5, np.min([all_Jindex+ all_Jattr])-0.05, 1.05])
plt.plot(all_Jindex, 'x--', label = 'Network Overlap', linewidth = 1, ms = 4.5)
plt.plot(all_Jattr, 'x--', label = 'Attractor Overlap', linewidth = 1, ms = 4.5)
ax_.legend(fontsize = 6.5, loc = 'best')
plt.ylabel('Overlap (Jaccard Index)', fontsize = 6.5)
plt.xlabel('Number of Iterations', fontsize = 7)
plt.xticks([i for i in range(len(all_pos_edge) - 1)], fontsize = 6)
plt.yticks(fontsize = 5)
plt.xticks(fontsize = 6)
ax_.yaxis.set_label_coords(-(2/(len(network_info)*10)), 0.5)

Orig_g = original_netw[0]
pos = nx.shell_layout(Orig_g)
for j in range(len(all_pos_edge)):
    ax1 = plt.subplot2grid((23,len(all_pos_edge)), (12,j), rowspan=4)
    orig_pos = [i for i in all_pos_edge[0] if i not in all_pos_edge[j]]
    orig_neg = [i for i in all_neg_edge[0] if i not in all_neg_edge[j]]
    edges1 = all_pos_edge[j] + all_neg_edge[j] + orig_pos + orig_neg
    g = nx.Graph()
    g = nx.from_edgelist(edges1)
    edge_col = [Blacks_r_map[0] if edge not in all_neg_edge[0]+all_pos_edge[0] else Greys_r_map[0] if edge not in all_neg_edge[j]+all_pos_edge[j] else 'red' if edge in all_neg_edge[j] else 'blue' for edge in [tuple(sorted(item)) for item in g.edges()]]
    edge_style = ['dashed' if edge not in all_neg_edge[0]+all_pos_edge[0] else 'dashed' if edge not in all_neg_edge[j]+all_pos_edge[j] else 'solid' for edge in [tuple(sorted(item)) for item in g.edges()]]
    edge_width = [1 if edge not in all_neg_edge[0]+all_pos_edge[0] else 1 if edge not in all_neg_edge[j]+all_pos_edge[j] else .9 for edge in [tuple(sorted(item)) for item in g.edges()]]

    nx.draw_networkx_nodes(g,pos,
                            node_color = 'k',
                            node_size = 1,
                            alpha = 1,
                            ax = ax1
                            )

    nx.draw_networkx_edges(g,pos,
                            width = edge_width,
                            alpha = .9,
                            edge_color = edge_col,
                            style = edge_style,
                            ax = ax1)
    if j == 0:
        ax1.set_title(r'Original Network ($G_{0}$)', fontsize = 6, va = 'top')
    if j == 1:
        ax1.set_title(r'ESABO Network ($G^{*}$)', fontsize = 6, va = 'top')

for j in range(len(all_pos_edge)):
    ax2 = plt.subplot2grid((23,len(all_pos_edge)), (17,j), rowspan=4)
    sns.heatmap(Attr_combined[j], mask = same_attr[j], cmap = Greens_map, cbar = False, xticklabels= False, yticklabels=False)
    sns.heatmap(Attr_combined[j], mask = diff_attr[j], cmap = Blacks_map, cbar = False, xticklabels= False, yticklabels=False)
    
    if j == 0:
        ax2.set_title(r'Original Attractors ($A_{0}$)', fontsize = 6, va = 'top')

    if j == 1:
        ax2.set_title(r'ESABO Attractors ($A^{*}$)', fontsize = 6, va = 'top')

for j in range(len(all_pos_edge)):
    ax3 = plt.subplot2grid((23,len(all_pos_edge)), (21,j), rowspan=2)
    sns.heatmap(missing_attr[j], mask = init[j], cmap = Greys_map, cbar = False, xticklabels= False, yticklabels=False)
    sns.heatmap(missing_attr[j], mask = ones_s[j], cmap = 'binary', cbar = False, xticklabels= False, yticklabels=False)

plt.figtext(.12, .139, f'Missing attractors'+r'$\rightarrow$', fontsize = 7.5)

plt.subplots_adjust(wspace = 0.3)
plt.savefig('Network_Evolution_from_ESABO_to_Original.png', dpi = 300, facecolor = 'white', bbox_inches = 'tight')
plt.close()