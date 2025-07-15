# coding: cp1252

import pandas as pd
import networkx as nx
import numpy as np
import pickle as pk
from collections import Counter

def degree_preserving_randomization(G):
    import random
    '''
    Randomizes a network provided by an edge list 
    producing neither self links nor duplicate links.
    The degree sequence will stay the same.
    INPUT:
    --- network
         
    OUTPUT:
    --- randomized network
    '''
    
    # make new set copy from edgelist
    edges = set( [tuple(e) for e in list(G.edges()) ]) 

    # get list of stubs
    stubs = [ ]
    [ stubs.extend(e) for e in edges ]

    # get a Counter object that counts the stubs for every node
    stub_counter = dict(Counter(stubs))

    # initialize the new edge list
    new_edges = set()

    # get available nodes (nodes that have nonzero stub count)
    nodes = np.array([ stub for stub,count in stub_counter.items() if count!=0 ])

    # loop till the number of available nodes is zero
    while len(nodes)>0:

        # initialize dummy values for new edge
        first,second = -1,-1

        # choose edges that are not self-links (only possible if len(nodes)>1)
        while first == second and len(nodes)>1:
            first,second = np.random.choice(nodes,size=(2,),replace=False)

        # if the chosen (source,target) is are not the same
        # and not yet connected 
        # and there is more than one node with available stubs
        if first!=second and \
           (first,second) not in new_edges and \
           (second,first) not in new_edges and \
           len(nodes)>1:
            new_edges.add((first,second))
            stub_counter[first] -= 1
            stub_counter[second] -= 1
        else:
            # if not, pop a random edge and put its nodes 
            # back in the stub pool
            edge = random.sample(new_edges,1)[0]
            new_edges.remove(edge)
            stub_counter[edge[0]] += 1
            stub_counter[edge[1]] += 1

        # get available nodes (nodes that have nonzero stub count)
        nodes = np.array([ stub for stub,count in stub_counter.items() if count!=0 ])
        
    new_G = nx.Graph()
    for edge in new_edges:
        new_G.add_edge(edge[0],edge[1])
        
    for node in G.nodes():
        if node not in nodes:
            new_G.add_node(node)
    return new_G

    

#Let's import the network
backbone_ss_exposure_network = nx.read_weighted_edgelist("output/backbone_exp_graph_significant_weighted.edgelist")
backbone_ss_gene_network = nx.read_weighted_edgelist("output/backbone_gene_graph_significant_weighted.edgelist")

random_een = degree_preserving_randomization(backbone_ss_exposure_network)
random_ggn = degree_preserving_randomization(backbone_ss_gene_network)

random_een_df = nx.to_pandas_edgelist(random_een,'Exp A', 'Exp B')
random_een_df.to_csv('output/random_een.tsv', sep = '\t')

random_ggn_df = nx.to_pandas_edgelist(random_ggn,'Gene A', 'Gene B')
random_ggn_df.to_csv('output/random_ggn.tsv', sep = '\t')


#To check whether the EEN has a hierarchical structure, we will compute the clusering coefficient for each node 
#and correlated with the degree in the EEN and comparing with a random rewiring network


EEN_clustering_coefficient_node=[]
EEN_degree_node=[]
random_rewired_EEN_clustering_coefficient_node=[]
EEN_predicted_clustering_coefficient=[]

EEN_degree_dict=dict(nx.degree(backbone_ss_exposure_network))
EEN_degree_dict_sorted=dict(sorted(EEN_degree_dict.items(), key=lambda item: item[1]))

for node,degree in EEN_degree_dict_sorted.items():
    EEN_degree_node.append(degree)
    EEN_predicted_clustering_coefficient.append(1/degree)
    EEN_clustering_coefficient_node.append(nx.clustering(backbone_ss_exposure_network, node))
    random_rewired_EEN_clustering_coefficient_node.append(nx.clustering(random_een, node))

EEN_characteristics_dict = {}
EEN_characteristics_dict["EEN_degree"] = EEN_degree_node
EEN_characteristics_dict["EEN_clustering_coefficient"] = EEN_clustering_coefficient_node
EEN_characteristics_dict["random_EEN_clustering_coefficient"] = random_rewired_EEN_clustering_coefficient_node
EEN_characteristics_dict["predicted_EEN_clustering_coefficient"] = EEN_predicted_clustering_coefficient

with open('output/EEN_characteristics_dict.pickle', 'wb') as handle:
    pk.dump(EEN_characteristics_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

GGN_clustering_coefficient_node=[]
GGN_degree_node=[]
random_rewired_GGN_clustering_coefficient_node=[]
GGN_predicted_clustering_coefficient=[]

GGN_degree_dict = dict(nx.degree(backbone_ss_gene_network))
GGN_degree_dict_sorted = dict(sorted(GGN_degree_dict.items(), key=lambda item: item[1]))

for node,degree in GGN_degree_dict_sorted.items():
    GGN_degree_node.append(degree)
    GGN_predicted_clustering_coefficient.append(1/degree)
    GGN_clustering_coefficient_node.append(nx.clustering(backbone_ss_gene_network, node))
    random_rewired_GGN_clustering_coefficient_node.append(nx.clustering(random_ggn, node))


GGN_characteristics_dict = {}
GGN_characteristics_dict["GGN_degree"] = GGN_degree_node
GGN_characteristics_dict["GGN_clustering_coefficient"] = GGN_clustering_coefficient_node
GGN_characteristics_dict["random_GGN_clustering_coefficient"] = random_rewired_GGN_clustering_coefficient_node
GGN_characteristics_dict["predicted_GGN_clustering_coefficient"] = GGN_predicted_clustering_coefficient


with open('output/GGN_characteristics_dict.pickle', 'wb') as handle:
    pk.dump(GGN_characteristics_dict, handle, protocol=pk.HIGHEST_PROTOCOL)


