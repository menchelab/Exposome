import random as rd
import pandas as pd
import networkx as nx
import numpy as np
import pickle as pk
import re
import itertools as it
import math
from collections import (defaultdict,Counter)
import obonet

#Let's import the communities
with open('output/Communities/Louvain/weighted_ji_lead_genes_gobp_fine_grained_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    fine_grained_communities_gobp = pk.load(handle)
with open('output/Communities/Louvain/weighted_ji_lead_genes_gobp_middle_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    middle_communities_gobp = pk.load(handle)
with open('output/Communities/Louvain/weighted_ji_lead_genes_gobp_broad_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    broad_communities_gobp = pk.load(handle)

with open('output/Communities/Louvain/weighted_ji_lead_genes_gomf_fine_grained_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    fine_grained_communities_gomf = pk.load(handle)
with open('output/Communities/Louvain/weighted_ji_lead_genes_gomf_middle_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    middle_communities_gomf = pk.load(handle)
with open('output/Communities/Louvain/weighted_ji_lead_genes_gomf_broad_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    broad_communities_gomf = pk.load(handle)

with open('output/Communities/Louvain/weighted_ji_lead_genes_gocc_fine_grained_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    fine_grained_communities_gocc = pk.load(handle)
with open('output/Communities/Louvain/weighted_ji_lead_genes_gocc_middle_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    middle_communities_gocc = pk.load(handle)
with open('output/Communities/Louvain/weighted_ji_lead_genes_gocc_broad_louvain_cluster_significant_fdr.pickle', 'rb') as handle:
    broad_communities_gocc = pk.load(handle)

#Let's import the Go Resnik dictionary
with open('output/GO_resnik.pickle', 'rb') as handle:
    GO_resnik = pk.load(handle)

#Let's import the ontology
url = 'http://purl.obolibrary.org/obo/go.obo'
graph = obonet.read_obo(url)
graph_up = nx.DiGraph.reverse(graph)

#Let's divide it into the three brunches
GOMF = graph_up.subgraph(list(nx.descendants(graph_up,'GO:0003674')))
GOCC = graph_up.subgraph(list(nx.descendants(graph_up,'GO:0005575')))
GOBP = graph_up.subgraph(list(nx.descendants(graph_up,'GO:0008150')))

#Here, we define a function that allows us to identify enrichment terms in the GO tree
def find_enclosed(s):
    import re
    # find all matches
    matches = re.findall(r"\((.*?)\)", s)
    matches = [x for x in matches if str(x).startswith("GO:")]
    # if there are no matches return None
    if len(matches) == 0:
        return None
    # if it is a valid number change its type to a number
    for i in range(len(matches)):
        try:
            matches[i] = int(matches[i])
        except:
            pass
    # if there is only one match return it without a list
    if len(matches) ==  1:
        return matches[0]
    return matches


def terms_extraction(enriched_terms):
    id_list1=[]
    for el in enriched_terms:
        id_list1.append(find_enclosed(el))
    return id_list1

#This function is used to calculate the resnik values between a list of terms
def resnik_terms_calculation(terms,resnik_dict):
    resnik_id_list1=[]
    pairwise_list = list(it.combinations(terms, 2))
    for i,j in pairwise_list:
        try:
            resnik_id_list1.append(round(resnik_dict[i,j],3))
        except:
            try:
                resnik_id_list1.append(round(resnik_dict[j,i],3))
            except:
                pass
    return resnik_id_list1

#Let's calculate the resnik for each community

fine_grained_communities_gomf_resnik={}
for com,term_list in fine_grained_communities_gomf.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        fine_grained_communities_gomf_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)
middle_communities_gomf_resnik={}
for com,term_list in middle_communities_gomf.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        middle_communities_gomf_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)
broad_communities_gomf_resnik={}
for com,term_list in broad_communities_gomf.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        broad_communities_gomf_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)

fine_grained_communities_gocc_resnik={}
for com,term_list in fine_grained_communities_gocc.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        fine_grained_communities_gocc_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)
middle_communities_gocc_resnik={}
for com,term_list in middle_communities_gocc.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        middle_communities_gocc_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)
broad_communities_gocc_resnik={}
for com,term_list in broad_communities_gocc.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        broad_communities_gocc_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)

fine_grained_communities_gobp_resnik={}
for com,term_list in fine_grained_communities_gobp.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        fine_grained_communities_gobp_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)
middle_communities_gobp_resnik={}
for com,term_list in middle_communities_gobp.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        middle_communities_gobp_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)
broad_communities_gobp_resnik={}
for com,term_list in broad_communities_gobp.items():
    if type(term_list[1])!=str:
        enriched_terms_id=terms_extraction(term_list[1])
        broad_communities_gobp_resnik[com]=resnik_terms_calculation(enriched_terms_id,GO_resnik)

#Let's write the results

with open('output/Communities/Louvain/weighted_fine_grained_communities_gomf_resnik.pickle', 'wb') as handle:
    pk.dump(fine_grained_communities_gomf_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_middle_communities_gomf_resnik.pickle', 'wb') as handle:
    pk.dump(middle_communities_gomf_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_broad_communities_gomf_resnik.pickle', 'wb') as handle:
    pk.dump(broad_communities_gomf_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)

with open('output/Communities/Louvain/weighted_fine_grained_communities_gocc_resnik.pickle', 'wb') as handle:
    pk.dump(fine_grained_communities_gocc_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_middle_communities_gocc_resnik.pickle', 'wb') as handle:
    pk.dump(middle_communities_gocc_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_broad_communities_gocc_resnik.pickle', 'wb') as handle:
    pk.dump(broad_communities_gocc_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)

with open('output/Communities/Louvain/weighted_fine_grained_communities_gobp_resnik.pickle', 'wb') as handle:
    pk.dump(fine_grained_communities_gobp_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_middle_communities_gobp_resnik.pickle', 'wb') as handle:
    pk.dump(middle_communities_gobp_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_broad_communities_gobp_resnik.pickle', 'wb') as handle:
    pk.dump(broad_communities_gobp_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)


#Let's go for calculating the random expectation
#This function will create a random expectation for a list of enriched terms
#it returns the avg random expected resnik values and a z-score for significance

def random_resnik_expectation(terms, resnik_id_list1,resnik_dict, ontology_tree):   #it needs the list of enriched terms, calculated resnik and HPO tree
    term_sample_size = len(set(terms))
    avg_random_resnik_list = []
    if term_sample_size>2:
        while len(avg_random_resnik_list) < 1000:
            term_sample = rd.sample(list(ontology_tree.nodes()), term_sample_size)
            resnik_randomid_list1 = []
            pairwise_list = list(it.combinations(term_sample, 2))
            for i, j in pairwise_list:
                try:
                    try:
                        resnik_randomid_list1.append(round(resnik_dict[i,j],3))
                    except:
                        resnik_randomid_list1.append(round(resnik_dict[j,i],3))
                except:
                    pass
            resnik_randomid_list1_cleaned = [x for x in resnik_randomid_list1 if str(x) != 'nan']
            if str(np.mean(resnik_randomid_list1_cleaned)) != 'nan':
                avg_random_resnik_list.append(np.mean(resnik_randomid_list1_cleaned))
            else:
                pass
        mu = np.mean(avg_random_resnik_list)
        std = np.std(avg_random_resnik_list)
        z = (np.mean(resnik_id_list1) - mu) / std
    else:
        z = 'nan'
    return avg_random_resnik_list,z

rd_fine_grained_communities_gomf_resnik={}
for com,term_list in fine_grained_communities_gomf.items():
    if type(term_list[1])!=str:
        rd_fine_grained_communities_gomf_resnik[com]=random_resnik_expectation(term_list[1],fine_grained_communities_gomf_resnik[com],GO_resnik,GOMF)

rd_middle_communities_gomf_resnik={}
for com,term_list in middle_communities_gomf.items():
    if type(term_list[1])!=str:
        rd_middle_communities_gomf_resnik[com]=random_resnik_expectation(term_list[1],middle_communities_gomf_resnik[com],GO_resnik,GOMF)

rd_broad_communities_gomf_resnik={}
for com,term_list in broad_communities_gomf.items():
    if type(term_list[1])!=str:
        rd_broad_communities_gomf_resnik[com]=random_resnik_expectation(term_list[1],broad_communities_gomf_resnik[com],GO_resnik,GOMF)

rd_fine_grained_communities_gocc_resnik={}
for com,term_list in fine_grained_communities_gocc.items():
    if type(term_list[1])!=str:
        rd_fine_grained_communities_gocc_resnik[com]=random_resnik_expectation(term_list[1],fine_grained_communities_gocc_resnik[com],GO_resnik,GOCC)

rd_middle_communities_gocc_resnik={}
for com,term_list in middle_communities_gocc.items():
    if type(term_list[1])!=str:
        rd_middle_communities_gocc_resnik[com]=random_resnik_expectation(term_list[1],middle_communities_gocc_resnik[com],GO_resnik,GOCC)

rd_broad_communities_gocc_resnik={}
for com,term_list in broad_communities_gocc.items():
    if type(term_list[1])!=str:
        rd_broad_communities_gocc_resnik[com]=random_resnik_expectation(term_list[1],broad_communities_gocc_resnik[com],GO_resnik,GOCC)

rd_fine_grained_communities_gobp_resnik={}
for com,term_list in fine_grained_communities_gobp.items():
    if type(term_list[1])!=str:
        rd_fine_grained_communities_gobp_resnik[com]=random_resnik_expectation(term_list[1],fine_grained_communities_gobp_resnik[com],GO_resnik,GOBP)

rd_middle_communities_gobp_resnik={}
for com,term_list in middle_communities_gobp.items():
    if type(term_list[1])!=str:
        rd_middle_communities_gobp_resnik[com]=random_resnik_expectation(term_list[1],middle_communities_gobp_resnik[com],GO_resnik,GOBP)

rd_broad_communities_gobp_resnik={}
for com,term_list in broad_communities_gobp.items():
    if type(term_list[1])!=str:
        rd_broad_communities_gobp_resnik[com]=random_resnik_expectation(term_list[1],broad_communities_gobp_resnik[com],GO_resnik,GOBP)

with open('output/Communities/Louvain/weighted_rd_fine_grained_communities_gomf_resnik.pickle', 'wb') as handle:
    pk.dump(rd_fine_grained_communities_gomf_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_middle_communities_gomf_resnik.pickle', 'wb') as handle:
    pk.dump(rd_middle_communities_gomf_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_broad_communities_gomf_resnik.pickle', 'wb') as handle:
    pk.dump(rd_broad_communities_gomf_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_fine_grained_communities_gocc_resnik.pickle', 'wb') as handle:
    pk.dump(rd_fine_grained_communities_gocc_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_middle_communities_gocc_resnik.pickle', 'wb') as handle:
    pk.dump(rd_middle_communities_gocc_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_broad_communities_gocc_resnik.pickle', 'wb') as handle:
    pk.dump(rd_broad_communities_gocc_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_fine_grained_communities_gobp_resnik.pickle', 'wb') as handle:
    pk.dump(rd_fine_grained_communities_gobp_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_middle_communities_gobp_resnik.pickle', 'wb') as handle:
    pk.dump(rd_middle_communities_gobp_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/weighted_rd_broad_communities_gobp_resnik.pickle', 'wb') as handle:
    pk.dump(rd_broad_communities_gobp_resnik, handle, protocol=pk.HIGHEST_PROTOCOL)
