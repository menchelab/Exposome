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

#This function performs a fisher test for the exposures
def overlap_set_exposures(setA, setB,tot_gene):    
    """
    Accepts to lists
    M is the population size (previously N)
    n is the number of successes in the population
    N is the sample size (previously n)
    x is still the number of drawn “successes”
    """
    M=  tot_gene #total number of genes
    n= len(setA)
    N= len(setB)
    x= len(setA.intersection(setB))
    return str(stats.hypergeom.sf(x-1, M, n, N))

def overlap_jaccard(list1,list2):
    intersction_term= len(set(list1) & set(list2))
    denominator = len(set(list1).union(set(list2)))
    overlap_jaccard_coeff = intersction_term/denominator
    return overlap_jaccard_coeff

#This function allows us to correct for multiple hypothesis
from statsmodels.sandbox.stats.multicomp import multipletests
import math
def fdr_adjustment(list_of_pvals):
    return multipletests(list_of_pvals,method='fdr_bh')[1] #the benjamin hochberg method is used

#Let's import the random bipartite network
with open('output/chem_to_genes_null.pickle', 'rb') as handle:
    chem_to_genes_null = pk.load(handle)

chem_gene_cleaned_random = {}
tot_chem_random_list=[]
for chem,genelist in chem_to_genes_null.items():
    chem_gene_cleaned_random[chem] = set(genelist)
    tot_chem_random_list.append(chem)

gene_chem_random_cleaned = {}
for chem,geneset in chem_gene_cleaned_random.items():
    for g in geneset:
        try:
            gene_chem_random_cleaned[g] |= {chem}
        except KeyError as e:
            gene_chem_random_cleaned[g] = set([chem])


#Let's define the background for future hypergeometric test (the total number of genes in our dataset)
tot_exp_hq = len(set(tot_chem_random_list))

#Here we calculate all the statistical significance of the overlap between all possible exposure pairs with a 
#Fisher test
random_combination_list=list(it.combinations(gene_chem_random_cleaned.keys(),2))
gene_random_graph_fisher={}
for i in random_combination_list:
    fisher_test=overlap_set_exposures(gene_chem_random_cleaned[i[0]],gene_chem_random_cleaned[i[1]],tot_exp_hq)
    gene_random_graph_fisher[i]=float(fisher_test)

with open('intermediate/gene_random_graph_fisher.pickle', 'wb') as handle:
    pk.dump(gene_random_graph_fisher, handle, protocol=pk.HIGHEST_PROTOCOL)

#Here, we define the unweighted projected graph
unweighted_random_gene_graph_significant=nx.Graph()
for i in random_combination_list:
    ji = overlap_jaccard(gene_chem_random_cleaned[i[0]],gene_chem_random_cleaned[i[1]])
    if ji>0:     #this ensures that there is a shared gene in the original bipartite graph
        unweighted_random_gene_graph_significant.add_edge(*i)
        unweighted_random_gene_graph_significant[i[0]][i[1]]['weight']=ji
        
nx.write_weighted_edgelist(unweighted_random_gene_graph_significant, "intermediate/unfiltered_weighted_random_gene_graph_significant.edgelist")

fdr_threshold=[0.05,0.01,0.001,0.0001,0.00001]   #different alpha threshold
alpha_backbone_threshold=[0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.09,0.08,0.07,0.06,0.05,0.01]   #different alpha threshold

unfiltered_weighted_random_gene_graph_fisher_dict = {}                     #Let's extract the Fisher's exact test pvals from all pairs to keep only those
for e in unweighted_random_gene_graph_significant.edges():    #that are present in the EEN
    try:
        unfiltered_weighted_random_gene_graph_fisher_dict[e]=gene_random_graph_fisher[(e[0],e[1])]
    except:
        unfiltered_weighted_random_gene_graph_fisher_dict[e]=gene_random_graph_fisher[(e[1],e[0])]

gene_keys_list=list(unfiltered_weighted_random_gene_graph_fisher_dict.keys())
gene_pvalues_list=list(unfiltered_weighted_random_gene_graph_fisher_dict.values())
adj_pvals=fdr_adjustment(gene_pvalues_list)
exp_graph_fisher_adjusted={} #here we obtain a dictionary of pairs with adjusted p-values
for el in range(len(gene_keys_list)):
    exp_graph_fisher_adjusted[gene_keys_list[el]]=adj_pvals[el]
comb_threshold_dict_nodes={}
comb_threshold_dict_edges={}
for alpha_t in fdr_threshold:
    exp_graph_ji_significant={}  #here, we select only the statistically significant
    for chem,fdr in exp_graph_fisher_adjusted.items():
        if float(fdr)<alpha_t:
            ji= overlap_jaccard(gene_chem_random_cleaned[chem[0]],gene_chem_random_cleaned[chem[1]])      #we are computing the odds ratio
            exp_graph_ji_significant[chem]=ji                                                     #and used as an edge weight
        else:
            pass
    weighted_exp_graph_significant=nx.Graph()
    for gene,ji_score in exp_graph_ji_significant.items():
        weighted_exp_graph_significant.add_edge(*gene)
        weighted_exp_graph_significant[gene[0]][gene[1]]['weight']=ji_score
    weighted_exp_graph_significant_dif = disparity_filter(weighted_exp_graph_significant)
    for alpha_b in alpha_backbone_threshold:
        backbone_exp_graph_significant = nx.Graph([(u, v, d) for u, v, d in weighted_exp_graph_significant_dif.edges(data=True) if d['alpha'] < alpha_b])   #let's apply the backboning threshold approach
        comb_threshold_dict_nodes["fdr_"+str(alpha_t),"backbone_"+str(alpha_b)]=backbone_exp_graph_significant.number_of_nodes()
        comb_threshold_dict_edges["fdr_"+str(alpha_t),"backbone_"+str(alpha_b)]=backbone_exp_graph_significant.number_of_edges()


with open('output/post_rev_comb_threshold_dict_edges_random_GGN_ji.pickle', 'wb') as handle:
    pk.dump(comb_threshold_dict_edges, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/post_rev_comb_threshold_dict_nodes_random_GGN_ji.pickle', 'wb') as handle:
    pk.dump(comb_threshold_dict_nodes, handle, protocol=pk.HIGHEST_PROTOCOL)

