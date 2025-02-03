import numpy as np
import pandas as pd
import os
from collections import Counter
from itertools import *
import glob
import shutil
from natsort import natsorted
import sys

###################################################################################################################

# Finding Jaccard index
def Jaccard(li1, li2):
    '''takes lists of tuples'''
    intersection = set.intersection(set(li1), set(li2))
    union = set.union(set(li1), set(li2))
    return(len(intersection)/len(union))

###################################################################################################################

condd = sys.argv[1]
inp_dir = f'{os.getcwd()}/{condd}'
out_dir = f'{os.getcwd()}/combined'
#os.makedirs(out_dir, exist_ok = True)

###################################################################################################################

Exps = [i.split("\\")[-1].split(f"_{condd}")[0] for i in glob.glob(f'{inp_dir}/*relative.abundance.csv')]
p_edges = [4, 5, 6, 'noedge']
n_edges = [2, 'noedge']

###################################################################################################################

for exp in Exps:
    for edge_p in p_edges:
        for edge_n in n_edges:

            files = natsorted(glob.glob(f'{inp_dir}/Exp_{exp}_*_p{edge_p}_n{edge_n}_*_FinalNetworkIntMat.csv'))
            Iter_files = natsorted(glob.glob(f'{inp_dir}/Exp_{exp}_*_p{edge_p}_n{edge_n}_*_Iterations.csv'))
            #print(files)
            #print(Iter_files)
            
            if len(files) == 1:
                copy_files = glob.glob(f'{inp_dir}/Exp_{exp}_*_p{edge_p}_n{edge_n}_*_FinalNetworkIntMat.csv')
                for f in copy_files:
                    shutil.copy2(f, out_dir)
                
                iter_copy_files = glob.glob(f'{inp_dir}/Exp_{exp}_*_p{edge_p}_n{edge_n}_*_Iterations.csv')
                for f2 in iter_copy_files:
                    shutil.copy2(f2, out_dir)

            if len(files) > 1:
                dfs = []
                for f2 in files:
                    df = pd.read_csv(f2, index_col=0, header=0)
                    dfs.append(df)

                df_count = pd.DataFrame()
                for n1,i in enumerate(dfs):
                    count = 0
                    for n2,j in enumerate(dfs):
                        df_count.loc[n1,n2] = (i == j).all().all()
                np.fill_diagonal(df_count.values, np.nan)

                #print(df_count)

                if True not in [item for sublist in df_count.values for item in sublist]:

                    net_edges = {}
                    for en, net in enumerate(dfs):
                        pos_e = []
                        neg_e = []
                        for r in net.index:
                            for c in net.columns:
                                if int(r) < int(c):
                                    if net.loc[r,c] == 1:
                                        pos_e.append((str(r), c))
                                    if net.loc[r,c] == -1:
                                        neg_e.append((str(r), c))
                        net_edges[en] = [pos_e, neg_e]

                    J_df = pd.DataFrame(columns = list(net_edges.keys()), index = list(net_edges.keys()))
                    for k1,v1 in net_edges.items():
                        for k2, v2 in net_edges.items():
                            if k1 != k2:
                                J_df.loc[k1, k2] = (Jaccard(v1[0], v2[0]) + Jaccard(v1[1], v2[1])) / 2


                    max_val = J_df.max(axis = 1).max()
                    #print(f'Max --> {max_val}')
                    pop = [[l,m] for l in J_df.index for m in J_df.columns if (J_df.loc[l, m] == max_val)][0][0]

                    x = '_'.join(files[pop].split('\\')[-1].split('_')[:-2])

                    #print(x)
                
                else:
                    same_dict = {}
                    for ind in df_count.index:
                        same_dict[ind] = Counter(df_count.loc[ind, :])[True]

                    c_ounter = [k for k,v in same_dict.items() if v == max(same_dict.values())][0]
                    x = '_'.join(files[c_ounter].split('\\')[-1].split('_')[:-2])

                    #print(x)

                copy_files_ = glob.glob(f'{inp_dir}/{x}_*_FinalNetworkIntMat.csv')
                for f_ in copy_files_:
                    shutil.copy2(f_, out_dir)

                iter_copy_files = glob.glob(f'{inp_dir}/{x}_*_Iterations.csv')
                for f2 in iter_copy_files:
                    shutil.copy2(f2, out_dir)