import numpy as np
import networkx as nx
import itertools as it
import random as rd
import pickle as pk
import pandas as pd
import netmedpy

#Here, we define the PPI
ppi = pd.read_csv("input/PPI/autocore_ppi_symbol_lcc.csv",delimiter= ',',
           skipinitialspace=True)
G_ppi = nx.from_pandas_edgelist(ppi, 'symbol1', 'symbol2')

#Let's remove eventual nans and convert each string to index for speeding up the proces
node_to_index = {}
node_to_remove = []
i = 0
for n in G_ppi.nodes():
    if str(n)!='nan':
        node_to_index[n] = i
        i = i+1
    else:
        node_to_remove.append(n)

G_ppi.remove_nodes_from(node_to_remove)
G_ppi_relabel = nx.relabel_nodes(G_ppi, node_to_index)

#Let's import the pre-calculated spl dictionary
with open('intermediate/ppi_spl.pickle', 'rb') as handle:
    spl = pk.load(handle)

#Let's convert the distance metrics into indexes
spl_index = {}
for node_pair,val in spl.items():
    if str(node_pair[0])!='nan' and str(node_pair[1])!='nan':
        spl_index[(node_to_index[node_pair[0]],node_to_index[node_pair[1]])] = val

#Let's import the diseases
gene_associations = pd.read_csv("input/Disease/all_gene_disease_associations.tsv",
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
        if gene in G_ppi.nodes():
            gene_list.append(gene)
        else:
            pass
    if len(gene_list)>0:
        diseases_genes_associated_ppi[k]=gene_list
    else:
        pass

#Let's convert the gene name lists into indexes
diseases_genes_associated_ppi_index = {}
for dis,genelist in diseases_genes_associated_ppi.items():
    genelist_index = [node_to_index[gene] for gene in genelist]
    diseases_genes_associated_ppi_index[dis]=genelist_index


chem_gene_df = pd.read_csv("input/CTD/CTD_chem_gene_ixns.tsv",delimiter= '\t',
            skipinitialspace=True)

chem_homo = chem_gene_df[(chem_gene_df['Organism'] == 'Homo sapiens')]

chem_gene = {}
for i,v in chem_homo.iterrows():
    try:
        chem_gene[v["ChemicalID"]] |= {v["GeneSymbol"]}
    except KeyError as e:
        chem_gene[v["ChemicalID"]] = set([v["GeneSymbol"]])

#Here, we remove the elements that do not perturb any genes:
chem_gene_ppi_index = {}
for chem,geneset in chem_gene.items():
    intersected_geneset = set(geneset)&set(G_ppi.nodes())
    intersected_geneset_index = [node_to_index[gene] for gene in intersected_geneset]
    if len(intersected_geneset_index)>1:
        chem_gene_ppi_index[chem]=intersected_geneset_index


chem_disease_proximity_dict = {}
for chem,geneset1 in chem_gene_ppi_index.items():
    for dis,geneset2 in diseases_genes_associated_ppi_index.items():
        prx_dict = netmedpy.proximity(G_ppi_relabel, geneset1,geneset2, spl_index,
                                          null_model="degree_match",n_iter=1000,
                                          symmetric=False)
        chem_disease_proximity_dict[chem,dis] = prx_dict

with open('output/chem_disease_proximity_dict.pickle', 'wb') as handle:
    pk.dump(chem_disease_proximity_dict, handle, protocol=pk.HIGHEST_PROTOCOL)
