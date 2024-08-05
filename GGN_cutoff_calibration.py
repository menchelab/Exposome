import os
import sys
import numpy as np
import networkx as nx
import itertools as it
import random as rd
import statsmodels
import pickle as pk
import os.path
import pandas as pd
from collections import (defaultdict,Counter)
import time
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy import integrate

#Here we define the disparity filter for backboning our future networks.
def disparity_filter(G, weight='weight'):
    ''' Compute significance scores (alpha) for weighted edges in G as defined in Serrano et al. 2009
        Args
            G: Weighted NetworkX graph
        Returns
            Weighted graph with a significance score (alpha) assigned to each edge
        References
            M. A. Serrano et al. (2009) Extracting the Multiscale backbone of complex weighted networks. PNAS, 106:16, pp. 6483-6488.
    '''
    if nx.is_directed(G): #directed case
        N = nx.DiGraph()
        for u in G:
            k_out = G.out_degree(u)
            k_in = G.in_degree(u)
            if k_out > 1:
                sum_w_out = sum(np.absolute(G[u][v][weight]) for v in G.successors(u))
                for v in G.successors(u):
                    w = G[u][v][weight]
                    p_ij_out = float(np.absolute(w))/sum_w_out
                    alpha_ij_out = 1 - (k_out-1) * integrate.quad(lambda x: (1-x)**(k_out-2), 0, p_ij_out)[0]
                    N.add_edge(u, v, weight = w, alpha_out=float('%.4f' % alpha_ij_out))
            elif k_out == 1 and G.in_degree(G.successors(u)[0]) == 1:
                #we need to keep the connection as it is the only way to maintain the connectivity of the network
                v = G.successors(u)[0]
                w = G[u][v][weight]
                N.add_edge(u, v, weight = w, alpha_out=0., alpha_in=0.)
                #there is no need to do the same for the k_in, since the link is built already from the tail
            if k_in > 1:
                sum_w_in = sum(np.absolute(G[v][u][weight]) for v in G.predecessors(u))
                for v in G.predecessors(u):
                    w = G[v][u][weight]
                    p_ij_in = float(np.absolute(w))/sum_w_in
                    alpha_ij_in = 1 - (k_in-1) * integrate.quad(lambda x: (1-x)**(k_in-2), 0, p_ij_in)[0]
                    N.add_edge(v, u, weight = w, alpha_in=float('%.4f' % alpha_ij_in))
        return N
    else: #undirected case
        B = nx.Graph()
        for u in G:
            k = len(G[u])
            if k > 1:
                sum_w = sum(np.absolute(G[u][v][weight]) for v in G[u])
                for v in G[u]:
                    w = G[u][v][weight]
                    p_ij = float(np.absolute(w))/sum_w
                    alpha_ij = 1 - (k-1) * integrate.quad(lambda x: (1-x)**(k-2), 0, p_ij)[0]
                    B.add_edge(u, v, weight = w, alpha=float('%.4f' % alpha_ij))
        return B


def overlap_jaccard(list1,list2):
    intersction_term= len(set(list1) & set(list2))
    denominator = len(set(list1).union(set(list2)))
    overlap_jaccard_coeff = intersction_term/denominator
    return overlap_jaccard_coeff


#Let's import the chemical-gene interactions from CTD (downloaded on 5th April 2021)
chem_gene_df = pd.read_csv("input/CTD/CTD_chem_gene_ixns.tsv",delimiter= '\t', skipinitialspace=True)
#Here, we filter for only the interactions that regards the H. Sapiens
chem_homo = chem_gene_df[(chem_gene_df['Organism'] == 'Homo sapiens')]

#Let's build also the vice-versa: the gene-chem dictionary
gene_chem = {}
for i,v in chem_homo.iterrows():
    try:
        gene_chem[v["GeneSymbol"]] |= {v["ChemicalID"]}
    except KeyError as e:
        gene_chem[v["GeneSymbol"]] = set([v["ChemicalID"]])

#Let's remove the elements which have 0 exposures associated with
gene_chem_cleaned = {}
tot_exp_list=[]
for k,v in gene_chem.items():
    if len(v)>0:
        gene_chem_cleaned[k]=v
        for exp in v:
            tot_exp_list.append(exp)
    else:
        pass
gene_chem_cleaned = {k: gene_chem_cleaned[k] for k in gene_chem_cleaned if type(k)==str}

#Let's import the precomputed analysis
with open('intermediate/gene_graph_fisher.pickle', 'rb') as handle:
    gene_graph_fisher = pk.load(handle)

unfiltered_weighted_gene_graph_significant = nx.Graph()
unfiltered_weighted_gene_graph_fisher_dict = {}
for i,pval in gene_graph_fisher.items():
    ji = overlap_jaccard(gene_chem_cleaned[i[0]],gene_chem_cleaned[i[1]])
    if ji>0:
        unfiltered_weighted_gene_graph_significant.add_edge(*i)
        unfiltered_weighted_gene_graph_significant[i[0]][i[1]]['weight']=ji
        unfiltered_weighted_gene_graph_fisher_dict[i]=pval

nx.write_weighted_edgelist(unfiltered_weighted_gene_graph_significant, "output/unfiltered_weighted_gene_graph_significant.edgelist")


#I want to check here how it will changes the number of nodes and the number of edges, changing the threshold
#This function allows us to correct for multiple hypothesis
from statsmodels.sandbox.stats.multicomp import multipletests
import math
def fdr_adjustment(list_of_pvals):
    return multipletests(list_of_pvals,method='fdr_bh')[1] #the benjamin hochberg method is used

fdr_threshold=[0.05,0.01,0.001,0.0001,0.00001]   #different alpha threshold
alpha_backbone_threshold=[0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,0.01]   #different alpha threshold
gene_keys_list=list(unfiltered_weighted_gene_graph_fisher_dict.keys())
gene_pvalues_list=list(unfiltered_weighted_gene_graph_fisher_dict.values())
adj_pvals=fdr_adjustment(gene_pvalues_list)
gene_graph_fisher_adjusted={} #here we obtain a dictionary of pairs with adjusted p-values
for el in range(len(gene_keys_list)):
    gene_graph_fisher_adjusted[gene_keys_list[el]]=adj_pvals[el]
comb_threshold_dict_nodes={}
comb_threshold_dict_edges={}
for alpha_t in fdr_threshold:
    weighted_gene_graph_significant=nx.Graph()
    for gene_p,fdr in gene_graph_fisher_adjusted.items():
        if float(fdr)<alpha_t:
            weighted_gene_graph_significant.add_edge(*gene_p)
            weighted_gene_graph_significant[gene_p[0]][gene_p[1]]['weight']=unfiltered_weighted_gene_graph_significant[gene_p[0]][gene_p[1]]['weight']
    weighted_gene_graph_significant_dif = disparity_filter(weighted_gene_graph_significant)
    for alpha_b in alpha_backbone_threshold:
        backbone_gene_graph_significant = nx.Graph([(u, v, d) for u, v, d in weighted_gene_graph_significant_dif.edges(data=True) if d['alpha'] < alpha_b])   #let's apply the backboning threshold approach
        comb_threshold_dict_nodes["fdr_"+str(alpha_t),"backbone_"+str(alpha_b)]=backbone_gene_graph_significant.number_of_nodes()
        comb_threshold_dict_edges["fdr_"+str(alpha_t),"backbone_"+str(alpha_b)]=backbone_gene_graph_significant.number_of_edges()

#Let's write the results in the intermidiate folder
with open('intermediate/post_rev_comb_threshold_dict_edges_GGN_ji.pickle', 'wb') as handle:
    pk.dump(comb_threshold_dict_edges, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('intermediate/post_rev_comb_threshold_dict_nodes_GGN_ji.pickle', 'wb') as handle:
    pk.dump(comb_threshold_dict_nodes, handle, protocol=pk.HIGHEST_PROTOCOL)
