import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from itertools import *
import seaborn as sns
import glob

from matplotlib.colors import LogNorm, Normalize

#####################################################################################################################

# Finding Jaccard index
def Jaccard(li1, li2):
    '''takes lists of tuples'''
    intersection = set.intersection(set(li1), set(li2))
    union = set.union(set(li1), set(li2))
    return(len(intersection)/len(union))

#####################################################################################################################

specie_df = pd.read_csv(f'MatchSpecies.csv')

dir = os.getcwd()
out_dir = f"{dir}/Heatmaps"
os.makedirs(out_dir, exist_ok=True)

#####################################################################################################################

edge = 'yesedge'
p_edges = [4, 5, 6]
n_edges = [2]

for edge_p in p_edges:
    for edge_n in n_edges:

        healthy_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*healthy_FinalNetworkIntMat.csv')[0] for i in range(1,5)]
        T2D_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*T2D_FinalNetworkIntMat.csv')[0] for i in range(1,5)]
        CRC_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*CRC_FinalNetworkIntMat.csv')[0] for i in range(1,5)]
        IBD_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*IBD_FinalNetworkIntMat.csv')[0] for i in range(1,5)]
        
        ####################################################
        healthy_FinalIntMat_dfs = []
        for FinI_file in healthy_FinalInt_files:
            Fini_df = pd.read_csv(FinI_file, index_col = 0)
            Fini_df.columns = specie_df['Species']
            Fini_df.index = Fini_df.columns
            healthy_FinalIntMat_dfs.append(Fini_df)

        T2D_FinalIntMat_dfs = []
        for FinI_file in T2D_FinalInt_files:
            Fini_df = pd.read_csv(FinI_file, index_col = 0)
            Fini_df.columns = specie_df['Species']
            Fini_df.index = Fini_df.columns
            T2D_FinalIntMat_dfs.append(Fini_df)

        CRC_FinalIntMat_dfs = []
        for FinI_file in CRC_FinalInt_files:
            Fini_df = pd.read_csv(FinI_file, index_col = 0)
            Fini_df.columns = specie_df['Species']
            Fini_df.index = Fini_df.columns
            CRC_FinalIntMat_dfs.append(Fini_df)

        IBD_FinalIntMat_dfs = []
        for FinI_file in IBD_FinalInt_files:
            Fini_df = pd.read_csv(FinI_file, index_col = 0)
            Fini_df.columns = specie_df['Species']
            Fini_df.index = Fini_df.columns
            IBD_FinalIntMat_dfs.append(Fini_df)

        #################################################### 
        healthy_main_df = pd.DataFrame()
        for i in range(len(healthy_FinalIntMat_dfs)):
            df = healthy_FinalIntMat_dfs[i]
            df.index.name = None; df.columns.name = None
            df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
            df = df.stack().reset_index()
            df.columns = ['Species1','Species2',f'healthy{i}']
            if i == 0:
                healthy_main_df = df
            else:
                healthy_main_df = healthy_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

        T2D_main_df = pd.DataFrame()
        for i in range(len(T2D_FinalIntMat_dfs)):
            df = T2D_FinalIntMat_dfs[i]
            df.index.name = None; df.columns.name = None
            df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
            df = df.stack().reset_index()
            df.columns = ['Species1','Species2',f'T2D{i}']
            if i == 0:
                T2D_main_df = df
            else:
                T2D_main_df = T2D_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

        CRC_main_df = pd.DataFrame()
        for i in range(len(CRC_FinalIntMat_dfs)):
            df = CRC_FinalIntMat_dfs[i]
            df.index.name = None; df.columns.name = None
            df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
            df = df.stack().reset_index()
            df.columns = ['Species1','Species2',f'CRC{i}']
            if i == 0:
                CRC_main_df = df
            else:
                CRC_main_df = CRC_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

        IBD_main_df = pd.DataFrame()
        for i in range(len(IBD_FinalIntMat_dfs)):
            df = IBD_FinalIntMat_dfs[i]
            df.index.name = None; df.columns.name = None
            df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
            df = df.stack().reset_index()
            df.columns = ['Species1','Species2',f'IBD{i}']
            if i == 0:
                IBD_main_df = df
            else:
                IBD_main_df = IBD_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

        #################################################### Dataframe
        main_df = healthy_main_df.merge(T2D_main_df, on=['Species1', 'Species2'], how = 'outer').merge(CRC_main_df, on=['Species1', 'Species2'], how = 'outer').merge(IBD_main_df, on=['Species1', 'Species2'], how = 'outer')
        main_df['Species'] = main_df['Species1']+'-'+main_df['Species2']
        del main_df['Species1']
        del main_df['Species2']
        main_df.index = main_df['Species']
        del main_df['Species']

        df_m = main_df.copy(deep = True)

        #################################################### JI
        pos_indices = {}
        neg_indices = {}
        for col in df_m.columns:
            pos1 = list(df_m.loc[df_m[col] == 1].index)
            pos_indices[col] = [(i.split('-')[0],i.split('-')[1]) for i in pos1]
            neg1 = list(df_m.loc[df_m[col] == -1].index)
            neg_indices[col] = [(i.split('-')[0],i.split('-')[1]) for i in neg1]

        JI_df = pd.DataFrame(columns=main_df.columns, index=main_df.columns)
        for k1 in pos_indices.keys():
            for k2 in pos_indices.keys():
                if k1 != k2:    
                    pos_JI = Jaccard(pos_indices[k1], pos_indices[k2])
                    neg_JI = Jaccard(neg_indices[k1], neg_indices[k2])
                    JI_df.loc[k1][k2] = (pos_JI + neg_JI)/2
        JI_df = JI_df.apply(pd.to_numeric)
        np.fill_diagonal(JI_df.values, np.nan)

        #################################################### Figure
        sns.heatmap(JI_df, cmap = 'RdYlBu', linewidths=0.5, annot=False, center = 0.2) #norm=LogNorm()
        plt.axvline(x=4, color='black', linestyle='--')
        plt.axvline(x=8, color='black', linestyle='--')
        plt.axvline(x=12, color='black', linestyle='--')
        plt.axhline(y=4, color='black', linestyle='--')
        plt.axhline(y=8, color='black', linestyle='--')
        plt.axhline(y=12, color='black', linestyle='--')
        plt.title(f"EDAME_p{edge_p}-n{edge_n}")
        plt.savefig(f'{out_dir}/Heatmap_JIdf_p{edge_p}_n{edge_n}.png', dpi = 300, bbox_inches = 'tight')
        plt.close()

        #################################################### Save similarity
        JI_df.to_csv(f'{out_dir}/EDAME_JIdf_p{edge_p}_n{edge_n}.csv')

#####################################################################################################################

edge = 'noedge'
edge_p = 'noedge'
edge_n = 'noedge'

healthy_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*healthy_FinalNetworkIntMat.csv')[0] for i in range(1,5)]
T2D_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*T2D_FinalNetworkIntMat.csv')[0] for i in range(1,5)]
CRC_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*CRC_FinalNetworkIntMat.csv')[0] for i in range(1,5)]
IBD_FinalInt_files = [glob.glob(dir + f'/Exp_{i}_*{edge}_p{edge_p}_n{edge_n}*IBD_FinalNetworkIntMat.csv')[0] for i in range(1,5)]

####################################################
healthy_FinalIntMat_dfs = []
for FinI_file in healthy_FinalInt_files:
    Fini_df = pd.read_csv(FinI_file, index_col = 0)
    Fini_df.columns = specie_df['Species']
    Fini_df.index = Fini_df.columns
    healthy_FinalIntMat_dfs.append(Fini_df)

T2D_FinalIntMat_dfs = []
for FinI_file in T2D_FinalInt_files:
    Fini_df = pd.read_csv(FinI_file, index_col = 0)
    Fini_df.columns = specie_df['Species']
    Fini_df.index = Fini_df.columns
    T2D_FinalIntMat_dfs.append(Fini_df)

CRC_FinalIntMat_dfs = []
for FinI_file in CRC_FinalInt_files:
    Fini_df = pd.read_csv(FinI_file, index_col = 0)
    Fini_df.columns = specie_df['Species']
    Fini_df.index = Fini_df.columns
    CRC_FinalIntMat_dfs.append(Fini_df)

IBD_FinalIntMat_dfs = []
for FinI_file in IBD_FinalInt_files:
    Fini_df = pd.read_csv(FinI_file, index_col = 0)
    Fini_df.columns = specie_df['Species']
    Fini_df.index = Fini_df.columns
    IBD_FinalIntMat_dfs.append(Fini_df)

#################################################### 
healthy_main_df = pd.DataFrame()
for i in range(len(healthy_FinalIntMat_dfs)):
    df = healthy_FinalIntMat_dfs[i]
    df.index.name = None; df.columns.name = None
    df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
    df = df.stack().reset_index()
    df.columns = ['Species1','Species2',f'healthy{i}']
    if i == 0:
        healthy_main_df = df
    else:
        healthy_main_df = healthy_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

T2D_main_df = pd.DataFrame()
for i in range(len(T2D_FinalIntMat_dfs)):
    df = T2D_FinalIntMat_dfs[i]
    df.index.name = None; df.columns.name = None
    df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
    df = df.stack().reset_index()
    df.columns = ['Species1','Species2',f'T2D{i}']
    if i == 0:
        T2D_main_df = df
    else:
        T2D_main_df = T2D_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

CRC_main_df = pd.DataFrame()
for i in range(len(CRC_FinalIntMat_dfs)):
    df = CRC_FinalIntMat_dfs[i]
    df.index.name = None; df.columns.name = None
    df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
    df = df.stack().reset_index()
    df.columns = ['Species1','Species2',f'CRC{i}']
    if i == 0:
        CRC_main_df = df
    else:
        CRC_main_df = CRC_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

IBD_main_df = pd.DataFrame()
for i in range(len(IBD_FinalIntMat_dfs)):
    df = IBD_FinalIntMat_dfs[i]
    df.index.name = None; df.columns.name = None
    df = df.where(np.triu(np.ones(df.shape), +1).astype(bool))
    df = df.stack().reset_index()
    df.columns = ['Species1','Species2',f'IBD{i}']
    if i == 0:
        IBD_main_df = df
    else:
        IBD_main_df = IBD_main_df.merge(df, on=['Species1', 'Species2'], how = 'outer')

#################################################### Dataframe
main_df = healthy_main_df.merge(T2D_main_df, on=['Species1', 'Species2'], how = 'outer').merge(CRC_main_df, on=['Species1', 'Species2'], how = 'outer').merge(IBD_main_df, on=['Species1', 'Species2'], how = 'outer')
main_df['Species'] = main_df['Species1']+'-'+main_df['Species2']
del main_df['Species1']
del main_df['Species2']
main_df.index = main_df['Species']
del main_df['Species']

df_m = main_df.copy(deep = True)

#################################################### JI
pos_indices = {}
neg_indices = {}
for col in df_m.columns:
    pos1 = list(df_m.loc[df_m[col] == 1].index)
    pos_indices[col] = [(i.split('-')[0],i.split('-')[1]) for i in pos1]
    neg1 = list(df_m.loc[df_m[col] == -1].index)
    neg_indices[col] = [(i.split('-')[0],i.split('-')[1]) for i in neg1]

JI_df = pd.DataFrame(columns=main_df.columns, index=main_df.columns)
for k1 in pos_indices.keys():
    for k2 in pos_indices.keys():
        if k1 != k2:    
            pos_JI = Jaccard(pos_indices[k1], pos_indices[k2])
            neg_JI = Jaccard(neg_indices[k1], neg_indices[k2])
            JI_df.loc[k1][k2] = (pos_JI + neg_JI)/2
JI_df = JI_df.apply(pd.to_numeric)
np.fill_diagonal(JI_df.values, np.nan)

#################################################### Figure
sns.heatmap(JI_df, cmap = 'RdYlBu', linewidths=0.5, annot=False, center = 0.2) #PiYG #center = 0.2 , norm=LogNorm()
plt.axvline(x=4, color='black', linestyle='--')
plt.axvline(x=8, color='black', linestyle='--')
plt.axvline(x=12, color='black', linestyle='--')
plt.axhline(y=4, color='black', linestyle='--')
plt.axhline(y=8, color='black', linestyle='--')
plt.axhline(y=12, color='black', linestyle='--')
plt.title(f"EDAME_Max JI")
plt.savefig(f'{out_dir}/Heatmap_JIdf_p{edge_p}_n{edge_n}.png', dpi = 300, bbox_inches = 'tight')
plt.close()

#################################################### Save similarity
JI_df.to_csv(f'{out_dir}/EDAME_JIdf_p{edge_p}_n{edge_n}.csv')

#####################################################################################################################