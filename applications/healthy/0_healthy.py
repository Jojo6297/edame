# Importing necessary libraries for data handling, calculations, and operations
import numpy as np
import pandas as pd
import os
import random
import collections
from collections import Counter
from itertools import *
import copy
import sys

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# FUNCTIONS

def abd(df):
    '''Transforms the DataFrame into a binary abundance matrix, where each non-zero entry becomes 1.'''
    df = pd.DataFrame(df)
    df = df.fillna(0)  # Replace NaN values with 0
    for ind in df.index:
        for col in df.columns:
            if df.loc[ind, col] != 0:
                df.loc[ind, col] = 1  # Convert all non-zero entries to 1
    return df

def operation1(mat1, mat2):
    '''Generates attractors by iterating matrix multiplication until a stable or cyclic state is found.'''
    o_count = 0
    new_list = []  # Stores states to detect cycles
    while True:
        simulated_matrix = [None] * len(mat2)
        mat = list(np.dot(mat1, mat2))  # Matrix multiplication
        mat2 = list(mat2)
        new_list.append(mat2)
        # Evaluate elements to create binary states
        for i in range(len(mat)):
            if mat[i] > 0:
                simulated_matrix[i] = 1
            elif mat[i] == 0:
                simulated_matrix[i] = mat2[i]
            else:
                simulated_matrix[i] = 0
        o_count += 1
        # Stop iteration if an attractor (stable or cyclic) is found
        if mat2 == simulated_matrix or simulated_matrix in new_list:
            break
        mat2 = simulated_matrix
    return simulated_matrix

def entropy_temp(vec):
    '''Calculates the proportion of 1s in a binary vector.'''
    t1 = len(vec)
    ones_count = collections.Counter(vec)[1]
    p1 = ones_count / float(t1)
    return p1

# Define "AND" operation for vector pairs
def op(pair):
    pair = list(pair)
    return pair[0] * pair[1]

# This function computes a shift score to assess binary similarity, given a function and repetitions
def shift_score(twovec, op, rep=5):
    x1 = np.array(list(map(op, twovec)))
    x2 = entropy_temp(x1)
    p3 = entropy_temp(np.transpose(twovec)[0])
    p4 = entropy_temp(np.transpose(twovec)[1])
    x4 = (p3 * p4)
    x5 = np.sqrt(((p3 * p4) * (1 - (p3 * p4))) / (len(twovec)))
    if x5 == 0:
        return 0
    return (x2 - x4) / x5

# Generates a network adjacency matrix based on positive and negative edges
def gen_network(node_num, pos_edges, neg_edges):
    Adj_upper = pd.DataFrame(np.zeros((node_num, node_num)))
    for i in pos_edges:
        Adj_upper.loc[i[0]][i[1]] = 1
    for j in neg_edges:
        Adj_upper.loc[j[0]][j[1]] = -1
    return Adj_upper + np.transpose(Adj_upper)

# Generates unique attractors for a network by applying matrix operations
def gen_attr(sample_num, df_phyla, netw):
    final_matrices = []
    for i in range(sample_num):
        a2 = df_phyla[df_phyla.columns[i]]
        M2 = operation1(netw, a2)
        final_matrices.append(M2)
    uniq_attractor = [list(x) for x in set(tuple(i) for i in final_matrices)]
    return np.array(uniq_attractor)

# Calculates the Jaccard similarity index between two sets of tuples
def Jaccard(li1, li2):
    intersection = set.intersection(set(li1), set(li2))
    union = set.union(set(li1), set(li2))
    return len(intersection) / len(union)

# Define XOR operation for binary pairs
def mod2(pair):
    return 0 if pair[0] == pair[1] else 1

# Calculates the Hamming distance between two binary attractors
def hamming(a, b):
    return sum([mod2([a[i], b[i]]) for i in range(len(a))])

# Determines differences in attractors based on Hamming distance
def attr_difference(setA, setB):
    hamm_dist_dict = {}
    for item1 in setA:
        hamm_dist = [hamming(item1, item2) for item2 in setB]
        hamm_dist_dict[item1] = [setB[i] for i in range(len(hamm_dist)) if hamm_dist[i] == min(hamm_dist) and hamm_dist[i] != 0]

    set_A_B = []
    for key, val in hamm_dist_dict.items():
        set_A_B.append([tuple([k - v1 for k, v1 in zip(key, v)]) for v in val])

    return set_A_B

# Detects switching events based on node array sorting
def sorted_nodes_array(attr, attr_astr, species_num):
    cD = list(set(attr) - set(attr_astr))
    cC = list(set(attr_astr) - set(attr))
    a_astr = copy.deepcopy(attr_astr)

    cD_a_astr = attr_difference(cD, a_astr)
    a_cC = attr_difference(cC, attr)
    cD_cC = attr_difference(cD, cC)

    a_vec = pd.DataFrame([list(item) for sublist in cD_a_astr for item in sublist]).abs().sum(axis=0).tolist()
    b_vec = pd.DataFrame([list(item) for sublist in a_cC for item in sublist]).abs().sum(axis=0).tolist()
    c_vec = pd.DataFrame([list(item) for sublist in cD_cC for item in sublist]).abs().sum(axis=0).tolist()
    new_list_vec = [item for item in [a_vec, b_vec, c_vec] if len(item) > 0]

    plist = pd.DataFrame(new_list_vec).abs().sum(axis=0).tolist()
    switch_linedict = dict(zip(range(species_num), plist))
    return dict(sorted(switch_linedict.items(), key=lambda x: x[1], reverse=True))

def netw_to_edges(g):
    temp_upper = pd.DataFrame(np.triu(g, k=0))
    mask1 = (temp_upper == 1)
    mask2 = (temp_upper == -1)
    posE = mask1.stack()[mask1].index.tolist()
    negE = mask2.stack()[mask2].index.tolist()
    return posE, negE

P_c = list(set([i for i in permutations([0, 0, 1, 1, -1, -1], 2)]))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

edge = sys.argv[1]

if edge == 'yesedge':
    edge_p = int(sys.argv[2])
    edge_n = int(sys.argv[3])

if edge == 'noedge':
    edge_p = sys.argv[2]
    edge_n = sys.argv[3]

exp_num = sys.argv[4]

dir = os.getcwd()
pdf1 = pd.read_csv(f"{dir}/{exp_num}_healthy.relative.abundance.csv", index_col='phyla')

abd1 = abd(pdf1)
abd1 = abd1.astype('Int64')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

labels__ = dict(zip(range(abd1.shape[0]), abd1.index))

final_matrices = [list(abd1[item]) for item in abd1.columns]
uniq_attractor = [list(x) for x in set(tuple(i) for i in final_matrices)]
uniq_attractor = np.array(uniq_attractor)
a = [tuple(j) for j in uniq_attractor]

## Uncomment: If want to export initial attractors
pd.DataFrame(a).to_csv(f'Exp_{exp_num}_healthy_InitialAttractors.csv')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

no_of_species = abd1.shape[0]

if no_of_species <= 10:     # exhaustively enumerate all possible binary states
    phyla_df = pd.DataFrame(product([0,1],repeat=no_of_species), columns= [i for i in range(no_of_species)]).T
    no_of_samples = phyla_df.shape[1]
else:                       # If the number of species is higher than 10, choose random 1000 binary initial states
    phyla_df = pd.DataFrame(np.random.randint(0,2,size=(no_of_species, 1000)))
    no_of_samples = phyla_df.shape[1]

cont_counting = [[j for j in range(i)] for i in range(2,no_of_species)]
allnodes = []
for item in cont_counting:
    allnodes.append([item+[p] for p in range(len(item),no_of_species)])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Calculating the shift scores
scores = []                                            # to save shift scores #
Edge_list = []
score_sign_list = []

for i in range(no_of_species):
    for j in range(no_of_species):
        if j > i:
            scores.append(shift_score(uniq_attractor[:, [i,j]], op))
            Edge_list.append((i ,j))

score_sign_list = [[i,j] for i,j in zip(scores, Edge_list)]
score_sign_list = sorted(score_sign_list, key=lambda x: x[0], reverse = True)

neg_items = [i[0] for i in score_sign_list if i[0] < 0]
pos_items = [i[0] for i in score_sign_list if i[0] > 0]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## If we choose "noedge" option, then run ESABO iteratively for combinations of negative and positive thresholds 
## to find the one with highest attractor JI.
if edge == 'noedge':
    pos_threshold = [i[0] for i in score_sign_list if i[0] > 0][1:] + [0]
    neg_threshold = [0] + [i[0] for i in score_sign_list if i[0] < 0][:-1]
    
    thresholdlist = []
    for item1 in neg_threshold:
        for item2 in pos_threshold:
            thresholdlist.append((item1, item2))

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    A_Jacs = dict()
    len_attr = []
    for comb in thresholdlist:
        m_positive_zscores_edges = [i[1] for i in score_sign_list if i[0] > comb[1]]
        m_negative_zscores_edges = [i[1] for i in score_sign_list if i[0] < comb[0]]

        Int_mat2 = gen_network(no_of_species, m_positive_zscores_edges, m_negative_zscores_edges)
        Rec_final_matrices2 = gen_attr(no_of_samples, phyla_df, Int_mat2)
        A_astr = [tuple(j) for j in Rec_final_matrices2]

        J_Attr = Jaccard(a, A_astr)
        A_Jacs[comb] = [J_Attr, abs(len(A_astr)-len(a))]

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    A_Jacs1 =  dict(sorted(A_Jacs.items(), key = lambda x:x[1][0], reverse = True))
    negth = list(A_Jacs1.keys())[0][0]
    posth = list(A_Jacs1.keys())[0][1]

    m_positive_zscores_edges = [i[1] for i in score_sign_list if i[0] > posth]
    m_negative_zscores_edges = [i[1] for i in score_sign_list if i[0] < negth]

    if len(m_positive_zscores_edges) == 0:
        m_positive_zscores_edges = [i[1] for i in score_sign_list[:1]]
    if len(m_negative_zscores_edges) == 0:
        m_negative_zscores_edges = [j[1] for j in score_sign_list[-1:]]

else:
    if len(pos_items) == 0:
        m_positive_zscores_edges = []
    elif edge_p > len(pos_items):
        edge_p_copy = len(pos_items)
        m_positive_zscores_edges = [i[1] for i in score_sign_list[:edge_p_copy]]
    else:
        m_positive_zscores_edges = [i[1] for i in score_sign_list[:edge_p]]

    if len(neg_items) == 0:
        m_negative_zscores_edges = []
    elif edge_n > len(neg_items):
        edge_n_copy = len(neg_items)
        m_negative_zscores_edges = [j[1] for j in score_sign_list[-edge_n_copy:]]
    else:
        m_negative_zscores_edges = [j[1] for j in score_sign_list[-edge_n:]]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## ESABO Network

Int_mat2 = gen_network(no_of_species, m_positive_zscores_edges, m_negative_zscores_edges)

## Uncomment: If want to export ESABO network
Int_mat2.to_csv(f'Exp_{exp_num}_{edge}_p{edge_p}_n{edge_n}_healthy_ESABONetworkIntMat.csv')

## ________________________________________________ Generating unique attractors _________________________________________________

Rec_final_matrices2 = gen_attr(no_of_samples, phyla_df, Int_mat2)
A_astr = [tuple(j) for j in Rec_final_matrices2]
J_Attr = Jaccard(a, A_astr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

EdgeChange_dict = sorted_nodes_array(a, A_astr, no_of_species)
EdgeChange_dict_val = copy.deepcopy(EdgeChange_dict)
print(EdgeChange_dict_val)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Run EDAME algorithm
networkobtained_li = []
Attrobtained_li = []

init_n = 3
n = init_n

numnodes_li = []
edge_detected_li = []
sign_li = []
Js_list = []

new_Int_mat2 = Int_mat2.copy(deep = True)
Temp_orig_J = Jaccard(A_astr, a)
J_val = Temp_orig_J

J = 0
fail = 0
tempfail = 0
iter_count = 0
faillimit = int(sys.argv[5])

save_ = [Temp_orig_J]
save_li = []
Edgechanges = []

while n >= init_n  and n < no_of_species+1 and fail < faillimit:
    node_count = len(allnodes[n-3][0])
    for item__1 in allnodes[n-3]:
        if len(list(EdgeChange_dict.keys())) > 0:
            check_nodes = list(map(list(EdgeChange_dict.keys()).__getitem__, item__1))
        else:
            check_nodes = list(EdgeChange_dict.keys())[:n]

        print(check_nodes)
        check_edge = [i for i in combinations(check_nodes, 2)]
        edge_comb = [i for i in combinations(check_edge, 2)]
        if fail > 0:
            edge_comb = random.sample(edge_comb, len(edge_comb))
        for item in edge_comb:
            number = 0
            count = 0
            
            P_c = random.sample(P_c, len(P_c))
            for combo in P_c:
                G_t = new_Int_mat2.copy(deep = True)
                if G_t[item[0][0]][item[0][1]] == combo[0] and G_t[item[1][0]][item[1][1]] == combo[1]:
                    number += 1

                    for diffedge in [o for o in P_c if o != combo]:
                        G_t[item[0][0]][item[0][1]] = diffedge[0]
                        G_t[item[0][1]][item[0][0]] = diffedge[0]
                        G_t[item[1][0]][item[1][1]] = diffedge[1]
                        G_t[item[1][1]][item[1][0]] = diffedge[1]


                        if no_of_species > 10:
                            phyla_df = pd.DataFrame(np.random.randint(0,2,size=(no_of_species, 1000)))
                            no_of_samples = phyla_df.shape[1]

                        # with the modified edges in the network, create new attractors
                        temp_uniq_attractor = gen_attr(no_of_samples, phyla_df, G_t)
                        A_t = [tuple(j) for j in temp_uniq_attractor]

                        J = Jaccard(A_t, a)

                        if number > 5:
                            break

                        if J < Temp_orig_J:
                            continue

                        if J > Temp_orig_J and J != 1:
                            print(Temp_orig_J, J)

                            save_.append(J)
                            networkobtained = G_t
                            Attrobtained = A_t
                            
                            Js_list.append(J)

                            edge_detected_li.append(item)
                            numnodes_li.append(node_count)
                            sign_li.append(diffedge)
                            iter_count += 1
                            print(f'itercount --------> {iter_count}')
                            new_Int_mat2 = G_t.copy(deep = True)

                            EdgeChange_dict = sorted_nodes_array(a, A_t, no_of_species)
                            #print(EdgeChange_dict)
                            break
                        
                        if J == 1:
                            print(f'Network detected')
                            save_.append(J)
                            networkobtained = G_t
                            Attrobtained = A_t
                            J_attr_end = J

                            edge_detected_li.append(item)
                            numnodes_li.append(node_count)
                            sign_li.append(diffedge)
                            iter_count += 1
                            print(f'itercount --------> {iter_count}')

                            Dict1 = dict(zip(edge_detected_li, sign_li))
                            Edge_detected = [[[k[0], v[0]],[k[1], v[1]]] for k,v in Dict1.items()]
                            
                            num_of_nodes_used = numnodes_li
                            edge_detected = Edge_detected
                            sign = sign_li
                            IterCount = iter_count
                            Failed = fail

                            networkobtained_li.append(networkobtained)
                            Attrobtained_li.append(Attrobtained)
                            save_li.append(save_)
                            
                            break

                    if J > Temp_orig_J and J != 1:
                        break
                    elif J == 1:
                        break
                    else:
                        continue
                if number == 0:
                    count += 1
                if count > init_n:
                    break
            if J > Temp_orig_J and J != 1:
                break
            elif J == 1:
                break
            else:
                continue
        if J > Temp_orig_J and J!= 1:
            n = init_n
            break
        elif J == 1:
            break
        else:
            node_count += 1
            continue
    
    if J > Temp_orig_J and J!= 1:
        Temp_orig_J = J
        
    elif J == 1:
        break

    elif len(allnodes[n-3][0]) == no_of_species:
        networkobtained_li.append(networkobtained)
        Attrobtained_li.append(Attrobtained)
        save_li.append(save_)
        
        Dict1 = dict(zip(edge_detected_li, sign_li))
        Edge_detected = [[[k[0], v[0]],[k[1], v[1]]] for k,v in Dict1.items()]
        Edgechanges.append(Edge_detected)

        fail += 1
        tempfail += 1
        if fail == 1:
            print(f'failed Once')
        else:
            print(f'failed -> {fail} times')    
        
        numnodes_li = []
        edge_detected_li = []
        sign_li = []
        iter_count = 0

        Temp_orig_J = J_val
        J = 0
        EdgeChange_dict = copy.deepcopy(EdgeChange_dict_val)
        new_Int_mat2 = Int_mat2.copy(deep = True)

        if tempfail + init_n <= no_of_species:
            n = tempfail + init_n # failes once, start with 4 shuffled nodes
        else:
            tempfail = 2
            n = tempfail + init_n
        print(f'after fail n ======> {n}')

        save_ = [J_val]
        networkobtained = []
        Attrobtained = []
        continue
    else:
        n += 1
        continue

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Select all output networks which gave the same end attractor Jaccard index
max_val = max({max(i) for i in save_li})
test_max_num = [i for i, li in enumerate(save_li) if li[-1] == max_val]
print(f'Runs with the maximum attractor JI -> {test_max_num}')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for x in test_max_num:

    NetwObIntMat = networkobtained_li[x]
    AttrOb = Attrobtained_li[x]
    values = save_li[x]

    NetwObIntMat.to_csv(f'Exp_{exp_num}_{edge}_p{edge_p}_n{edge_n}_{x}_healthy_FinalNetworkIntMat.csv')

    ## Uncomment: If want to export attractors of the Final network
    pd.DataFrame(AttrOb).to_csv(f'Exp_{exp_num}_{edge}_p{edge_p}_n{edge_n}_{x}_healthy_FinalAttractors.csv')

    ## Export the iteration scores
    pd.DataFrame(values).to_csv(f'Exp_{exp_num}_{edge}_p{edge_p}_n{edge_n}_{x}_healthy_Iterations.csv', index = False)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~