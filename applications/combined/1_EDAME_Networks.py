import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import os
from itertools import *
import glob 
import networkx as nx

#####################################################################################################################

dir = os.getcwd()
out_dir = f"{dir}/CumulativeNetworks_AllStudies"
os.makedirs(out_dir, exist_ok=True)

#####################################################################################################################

conds = ['healthy', 'CRC', 'IBD', 'T2D']
pns = [[4,2], [5,2], [6,2], ['noedge','noedge']]

#####################################################################################################################

for cond in conds:
    for pn in pns:
        species_name = list(pd.read_csv('MatchSpecies.csv')['Species'])
        files = glob.glob(f'{dir}/*p{pn[0]}_n{pn[1]}*{cond}*IntMat.csv')
        it_files = glob.glob(f'{dir}/*p{pn[0]}_n{pn[1]}*{cond}*Iterations.csv')
        dfs = [pd.read_csv(f, index_col=0) for f in files]
        iters = [list(pd.read_csv(If)['0'])[-1]for If in it_files]

        dfs = [t*df for t,df in zip(iters, dfs)]

        nodelabels = dict(zip(range(len(species_name)), species_name))
        network_df = dfs[0] + dfs[1] + dfs[2] + dfs[3]

        network_df.columns = species_name
        network_df.index = network_df.columns
        network_df.to_csv(f'{out_dir}/EDAME_Network_{cond}_p{pn[0]}n{pn[1]}.csv')

        net_df = network_df.where(np.triu(np.ones(network_df.shape), +1).astype(bool))
        net_df = net_df.stack().reset_index()
        net_df.columns = ['source','target','weight']
        net_df = net_df[net_df['weight'] != 0]
        weighted_edges = [(i,j,k) for i,j,k in zip(net_df['source'], net_df['target'], net_df['weight'])]

        g = nx.Graph()
        g.add_nodes_from(nodelabels.values())
        pos = nx.circular_layout(g)
        g.add_weighted_edges_from(weighted_edges)

        widths = nx.get_edge_attributes(g, 'weight')
        widths_dict = {k:abs(v)**2.2 for k,v in widths.items()}
        cols =  [(k,'crimson') if v < 0 else (k,'royalblue') for k,v in nx.get_edge_attributes(g, 'weight').items()]
        cols_dict = dict(zip([i[0] for i in cols], [i[1] for i in cols]))

        plt.figure(figsize=(7,5.6))
        #nodes
        nx.draw_networkx_nodes(g
                               ,pos
                               ,node_size = 0
                               )
        #edges
        nx.draw_networkx_edges(g,pos,
                                edgelist = widths.keys(),
                                width = list(widths_dict.values()),
                                alpha = 0.8,
                                edge_color = list(cols_dict.values()),
                                style = '-'
                                )
        #labels
        #for labeling outside the node
        offset = .1
        pos_labels = {}
        keys = pos.keys()
        for key in keys:
            x, y = pos[key]
            pos_labels[key] = (x+offset, y+offset)
                                
        #Please note, the code below uses the original idea of re-calculating a dictionary of adjusted label positions per node.
        label_ratio = 1.0/7.0
        pos_labels = {} 
        #For each node in the Graph
        for aNode in g.nodes():
            #Get the node's position from the layout
            x,y = pos[aNode]
            #Get the node's neighbourhood
            N = g[aNode]
            #Find the centroid of the neighbourhood. The centroid is the average of the Neighbourhood's node's x and y coordinates respectively.
            #Please note: This could be optimised further
            cx = sum(map(lambda x:pos[x][0], N)) / len(pos)
            cy = sum(map(lambda x:pos[x][1], N)) / len(pos)
            #Get the centroid's 'direction' or 'slope'. That is, the direction TOWARDS the centroid FROM aNode.
            slopeY = (y-cy)
            slopeX = (x-cx)
            #Position the label at some distance along this line. Here, the label is positioned at about 1/7th of the distance.
            pos_labels[aNode] = (x+slopeX*label_ratio, y+slopeY*label_ratio)

        #Finally, redraw the labels at their new position.
        nx.draw_networkx_labels(g, pos=pos_labels, font_size=14, font_family='Times New Roman')

        # to get nodes inside boundaries
        axis = plt.gca()
        # maybe smaller factors work as well, but 1.1 works fine for this minimal example
        axis.set_xlim([1.3*x for x in axis.get_xlim()])
        axis.set_ylim([1.2*y for y in axis.get_ylim()])
        plt.title(f'EDAME Cummulative network for {cond} samples (p{pn[0]}-n{pn[1]})')

        plt.savefig(f'{out_dir}/EDAME_Network_{cond}_p{pn[0]}n{pn[1]}.png', dpi = 100)
        plt.close()

#####################################################################################################################