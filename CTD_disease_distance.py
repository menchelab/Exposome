import os
import sys
import numpy as np
import networkx as nx
import itertools as it
import random as rd
import scipy.stats as st
import os.path
import pandas as pd
from collections import (defaultdict,Counter)
import time
# import statsmodels.sandbox.stats.multicomp as mc
import matplotlib.pyplot as plt
# %matplotlib inline
import pickle as pk


chem_vocabolary = pd.read_csv("CTD_chemicals_cleaned.txt",delimiter= '\t',
           skipinitialspace=True)

chem_gene_df = pd.read_csv("CTD_chem_gene_ixns.tsv",delimiter= '\t',
            skipinitialspace=True)

chem_homo = chem_gene_df[(chem_gene_df['Organism'] == 'Homo sapiens')]

chem_gene = {}
for i,v in chem_homo.iterrows():
    try:
        chem_gene[v["ChemicalID"]] |= {v["GeneSymbol"]}
    except KeyError as e:
        chem_gene[v["ChemicalID"]] = set([v["GeneSymbol"]])

#Here, we remove the elements that do not perturb any genes:
chem_gene_cleaned = {}
tot_gene_list=[]
for k,v in chem_gene.items():
    if len(v)>0:
        chem_gene_cleaned[k]=v
        for gene in v:
            tot_gene_list.append(gene)
    else:
        pass


chem_voc_dict = {}
for i,v in chem_vocabolary.iterrows():
    try:
        chem_voc_dict[v["ChemicalID"]] = v["ParentIDs"].split("|")
    except:
        chem_voc_dict[v["ChemicalID"]] = v["ParentIDs"]

from collections.abc import Mapping
# Empty directed graph
G = nx.DiGraph()

for k in chem_voc_dict.keys():
    try:
        for v in chem_voc_dict[k]:
            G.add_edge(k, v)
    except:
        G.add_edge(k,"D")

classes_list=[]
for i in G.nodes():
    if len(nx.predecessor(G, i))==2:
        classes_list.append(i)

#Let's clean the list from the nan elements
classes_list_cleaned = [x for x in classes_list if str(x) != 'nan']


ppi = pd.read_csv("julia_symbol_lcc.csv",delimiter= ',',
           skipinitialspace=True)

G_ppi = nx.from_pandas_edgelist(ppi, 'symbol1', 'symbol2')

G_ppi_lcc = G_ppi.subgraph(max(nx.connected_components(G_ppi), key=len))  # extract lcc graph
print(G_ppi_lcc.number_of_nodes())
print(G_ppi_lcc.number_of_edges())

chem_gene_dictio_cleaned_ppi={}
for k,v in chem_gene_cleaned.items():
    new_list=[]
    for gene in v:
        if gene in G_ppi_lcc.nodes():
            new_list.append(gene)
        else:
            pass
    if len(new_list)>0:
        chem_gene_dictio_cleaned_ppi[k]=new_list
    else:
        pass


#Loading the gene associations
gene_associations = pd.read_csv("all_gene_disease_associations.tsv",
           delimiter= '\t',
           skipinitialspace=True)

gene_associations_filtered=gene_associations[gene_associations['score']>0.3]

diseases_genes_associated = {}
for i,v in gene_associations_filtered.iterrows():
    try:
        diseases_genes_associated[v["diseaseName"]].append(v["geneSymbol"])
    except KeyError:
        diseases_genes_associated[v["diseaseName"]] = [v["geneSymbol"]]

diseases_genes_associated_ppi={}
for k,v in diseases_genes_associated.items():
    gene_list=[]
    for gene in v:
        if gene in G_ppi_lcc.nodes():
            gene_list.append(gene)
        else:
            pass
    if len(gene_list)>0:
        diseases_genes_associated_ppi[k]=gene_list
    else:
        pass

#Let's import the spl dictionary
with open('ppi_spl.pickle', 'rb') as handle:
    spl = pk.load(handle)

def calculate_closest_distance(spl, nodes_from, nodes_to):
    values_outer = []
    for node_from in nodes_from:
        values = []
        for node_to in nodes_to:
            if node_from==node_to:
                val =0
            else:
                try:
                    val = spl[node_from,node_to]
                except:
                    val = spl[node_to,node_from]
            values.append(val)
        d = min(values)
        #print d,
        values_outer.append(d)
    d = np.mean(values_outer)
    #print d
    return d


exp_dis_distance={}

for exp in chem_gene_dictio_cleaned_ppi.keys():
    for dis in diseases_genes_associated_ppi.keys():
        exp_gene_list_cleaned=[x for x in chem_gene_dictio_cleaned_ppi[exp] if str(x) != 'nan']
        dis_gene_list_cleaned=[x for x in diseases_genes_associated_ppi[dis] if str(x) != 'nan']
        exp_dis_distance[exp,dis]=calculate_closest_distance(spl,exp_gene_list_cleaned,dis_gene_list_cleaned)


with open('exp_disease_distance.pickle', 'wb') as handle:
    pk.dump(exp_dis_distance, handle, protocol=pk.HIGHEST_PROTOCOL)
