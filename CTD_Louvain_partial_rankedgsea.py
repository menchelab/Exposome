import pandas as pd
import networkx as nx
import numpy as np
import pickle as pk
import gseapy as gp


#Let's import the chemical-gene interactions from CTD (downloaded on 5th April 2021)
chem_gene_df = pd.read_csv("input/CTD/CTD_chem_gene_ixns.tsv",delimiter= '\t', skipinitialspace=True)
#Here, we filter for only the interactions that regards the H. Sapiens
chem_homo = chem_gene_df[(chem_gene_df['Organism'] == 'Homo sapiens')]

#Let's import the network
final_backbone_exp_graph_significant_combo_df = pd.read_csv("output/final_backbone_exp_graph_significant_combo_df.tsv", sep="\t",index_col=0)
backbone_ss_exposure_network = nx.from_pandas_edgelist(final_backbone_exp_graph_significant_combo_df, 'Exp A', 'Exp B')

#This cells create a dictionary where each key is a chemical compound and the correspondent value is a genelist
chem_gene = {}
for i,v in chem_homo.iterrows():
    try:
        chem_gene[v["ChemicalID"]] |= {v["GeneSymbol"]}
    except KeyError as e:
        chem_gene[v["ChemicalID"]] = set([v["GeneSymbol"]])

#Here, we keep only the exposures which perturb at least one gene
chem_gene_cleaned = {}
tot_gene_list=[]
for k,v in chem_gene.items():
    if len(v)>0:
        chem_gene_cleaned[k]=v
        for gene in v:
            tot_gene_list.append(gene)
    else:
        pass


with open('output/Communities/Louvain/ee_first_louvain_iteration_exposures.pickle', 'rb') as handle:
    ee_first_louvain_iteration_exposures = pk.load(handle)

with open('output/Communities/Louvain/ee_second_louvain_iteration_exposures.pickle', 'rb') as handle:
    ee_second_louvain_iteration_exposures = pk.load(handle)

with open('output/Communities/Louvain/ee_third_louvain_iteration_exposures.pickle', 'rb') as handle:
    ee_third_louvain_iteration_exposures = pk.load(handle)


def check_existence(x):
    try:
        x
        a=1
    except:
        a=0
    return a


def enr_df_dict(gene_list):   #This function returns the enrichment df for GOBP,GOCC,GOMF,KeGG

    import gseapy as gp

    libraries=['GO_Biological_Process_2021']
    c = 0
    while c<10:
        try:

            geneset_to_enrich_GOBP_df = gp.enrichr(gene_list=gene_list,
                                 gene_sets=libraries,
                                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                                 outdir=None, # don't write to disk
                                )

        except:
            c=c+1
        if check_existence(geneset_to_enrich_GOBP_df)==1:
            break

    libraries=['GO_Molecular_Function_2021']
    c = 0
    while c<10:
        try:
            geneset_to_enrich_GOMF_df = gp.enrichr(gene_list=gene_list,
                                 gene_sets=libraries,
                                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                                 outdir=None, # don't write to disk
                                )

        except:
            c=c+1
        if check_existence(geneset_to_enrich_GOMF_df)==1:
            break

    libraries=['GO_Cellular_Component_2021']
    c = 0
    while c<10:
        try:
            geneset_to_enrich_GOCC_df = gp.enrichr(gene_list=gene_list,
                             gene_sets=libraries,
                             organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                             outdir=None, # don't write to disk
                            )

        except:
            c=c+1
        if check_existence(geneset_to_enrich_GOCC_df)==1:
            break


    libraries=['KEGG_2021_Human']
    c = 0
    while c<10:
        try:
            geneset_to_enrich_KEGG_df = gp.enrichr(gene_list=gene_list,
                             gene_sets=libraries,
                             organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                             outdir=None, # don't write to disk
                            )

        except:
            c=c+1
        if check_existence(geneset_to_enrich_KEGG_df)==1:
            break

    enrdf_dict={}
    enrdf_dict['KEGG']=geneset_to_enrich_KEGG_df
    enrdf_dict['GOBP']=geneset_to_enrich_GOBP_df
    enrdf_dict['GOMF']=geneset_to_enrich_GOMF_df
    enrdf_dict['GOCC']=geneset_to_enrich_GOCC_df
    return enrdf_dict

def enr_fdr_list(enrdf_dict):
    sig_dict={}
    for lib, enr_df in enrdf_dict.items():
        lib_significant_terms=[]
        for i,v in enr_df.res2d.iterrows():
            if v['Adjusted P-value']<0.05:
                lib_significant_terms.append(v['Term'])
        sig_dict[lib]=lib_significant_terms
    return sig_dict



def enr_ranking(pre_rank_df):    #This function calculates the enrichment terms using GSEA with ranking

    pre_rank_df_significant_fdr={}
    c=0
    while c<10:
        try:
            kegg_pre_res = gp.prerank(rnk=pre_rank_df, gene_sets='KEGG_2021_Human',
                             processes=4,min_size=5,max_size=5000,outdir=None,   #we impose at least 5 enriched common genes
                             permutation_num=1000) # reduce number to speed up testing
            pre_rank_df_significant_fdr['KEGG']=list(kegg_pre_res.res2d[kegg_pre_res.res2d['fdr']<0.05].index)

        except:
            c=c+1
        if check_existence(kegg_pre_res)==1:
            if 'fdr' in kegg_pre_res.res2d==True:
                break
    c=0
    while c<10:
        try:

            gobp_pre_res = gp.prerank(rnk=pre_rank_first_louvain_cluster_genelist_df[com], gene_sets='GO_Biological_Process_2021',
                             processes=4,min_size=5,max_size=5000,outdir=None,
                             permutation_num=1000) # reduce number to speed up testing
            pre_rank_df_significant_fdr["GOBP"]=list(gobp_pre_res.res2d[gobp_pre_res.res2d['fdr']<0.05].index)

        except:
            c=c+1
        if check_existence(gobp_pre_res)==1:
            if 'fdr' in gobp_pre_res.res2d==True:
                break
    c=0
    while c<10:
        try:
            gomf_pre_res = gp.prerank(rnk=pre_rank_first_louvain_cluster_genelist_df[com], gene_sets='GO_Molecular_Function_2021',
                             processes=4,min_size=5,max_size=5000,outdir=None,
                             permutation_num=1000) # reduce number to speed up testing
            pre_rank_df_significant_fdr["GOMF"]=list(gomf_pre_res.res2d[gomf_pre_res.res2d['fdr']<0.05].index)

        except:
            c=c+1
        if check_existence(gomf_pre_res)==1:
            if 'fdr' in gomf_pre_res.res2d==True:
                break
    c=0
    while c<10:
        try:
            gocc_pre_res = gp.prerank(rnk=pre_rank_first_louvain_cluster_genelist_df[com], gene_sets='GO_Cellular_Component_2021',
                             processes=4,min_size=5,max_size=5000,outdir=None,
                             permutation_num=1000) # reduce number to speed up testing
            pre_rank_df_significant_fdr["GOCC"]=list(gocc_pre_res.res2d[gocc_pre_res.res2d['fdr']<0.05].index)

        except:
            c=c+1
        if check_existence(gocc_pre_res)==1:
            if 'fdr' in gocc_pre_res.res2d==True:
                break

    return pre_rank_df_significant_fdr



def enr_ranking_leading_genes(pre_rank_df):   #This function returns the leading genes

    leading_genes_tot=[]
    c=0
    while c<10:
        try:
            kegg_pre_res = gp.prerank(rnk=pre_rank_df, gene_sets='KEGG_2021_Human',
                             processes=4,min_size=1,max_size=5000,outdir=None,   #we impose at least 5 enriched common genes
                             permutation_num=1000) # reduce number to speed up testing

            leading_gene_kegg_list=[]
            for gene_str in kegg_pre_res.res2d['ledge_genes'].tolist():
                leading_gene_kegg_list.extend(gene_str.split(";"))
            leading_genes_tot.extend(leading_gene_kegg_list)

        except:
            c=c+1
        if check_existence(kegg_pre_res)==1:
            if 'fdr' in kegg_pre_res.res2d:
                break
    c=0
    while c<10:
        try:

            gobp_pre_res = gp.prerank(rnk=pre_rank_df, gene_sets='GO_Biological_Process_2021',
                             processes=4,min_size=1,max_size=5000,outdir=None,
                             permutation_num=1000) # reduce number to speed up testing
            leading_gene_gobp_list=[]
            for gene_str in gobp_pre_res.res2d['ledge_genes'].tolist():
                leading_gene_gobp_list.extend(gene_str.split(";"))
            leading_genes_tot.extend(leading_gene_gobp_list)
        except:
            c=c+1
        if check_existence(gobp_pre_res)==1:
            if 'fdr' in gobp_pre_res.res2d:
                break
    c=0
    while c<10:
        try:
            gomf_pre_res = gp.prerank(rnk=pre_rank_df, gene_sets='GO_Molecular_Function_2021',
                             processes=4,min_size=1,max_size=5000,outdir=None,
                             permutation_num=1000) # reduce number to speed up testing
            leading_gene_gomf_list=[]
            for gene_str in gomf_pre_res.res2d['ledge_genes'].tolist():
                leading_gene_gomf_list.extend(gene_str.split(";"))
            leading_genes_tot.extend(leading_gene_gomf_list)
        except:
            c=c+1
        if check_existence(gomf_pre_res)==1:
            if 'fdr' in gomf_pre_res.res2d:
                break
    c=0
    while c<10:
        try:
            gocc_pre_res = gp.prerank(rnk=pre_rank_df, gene_sets='GO_Cellular_Component_2021',
                             processes=4,min_size=1,max_size=5000,outdir=None,
                             permutation_num=1000) # reduce number to speed up testing
            leading_gene_gocc_list=[]
            for gene_str in gocc_pre_res.res2d['ledge_genes'].tolist():
                leading_gene_gocc_list.extend(gene_str.split(";"))
            leading_genes_tot.extend(leading_gene_gocc_list)
        except:
            c=c+1
        if check_existence(gocc_pre_res)==1:
            if 'fdr' in gocc_pre_res.res2d:
                break

    return set(leading_genes_tot)

#I can also consider a pre-rank for making a choice on the gene cutoff
from collections import Counter
pre_rank_first_louvain_cluster_genelist={}
for com,expset in ee_first_louvain_iteration_exposures.items():
    gene_list=[]
    rank_genelist=[]
    G_sub = nx.subgraph(backbone_ss_exposure_network,expset)   #let's create the subgraph corresponding to the cluster
    for edge in G_sub.edges():
        overlap_geneset=chem_gene_cleaned[edge[0]] & chem_gene_cleaned[edge[1]]
        for gene in overlap_geneset:
            gene_list.append(gene)
    count_dict=Counter(gene_list)
    for gene in set(gene_list):
        rank_genelist.append([gene,count_dict[gene]])

    pre_rank_first_louvain_cluster_genelist[com]=rank_genelist

pre_rank_first_louvain_cluster_genelist_df={}
for com,genelist in pre_rank_first_louvain_cluster_genelist.items():
    rank_df=pd.DataFrame(columns=('gene', 'rank'))
    for n in range(len(genelist)):
        rank_df.loc[n] = [genelist[n][0], float(genelist[n][1])]
    rank_df=rank_df.sort_values(by=['rank'],ascending=False)
    pre_rank_first_louvain_cluster_genelist_df[com]=rank_df

pre_rank_second_louvain_cluster_genelist={}
for com,expset in ee_second_louvain_iteration_exposures.items():
    gene_list=[]
    rank_genelist=[]
    G_sub = nx.subgraph(backbone_ss_exposure_network,expset)   #let's create the subgraph corresponding to the cluster
    for edge in G_sub.edges():
        overlap_geneset=chem_gene_cleaned[edge[0]] & chem_gene_cleaned[edge[1]]
        for gene in overlap_geneset:
            gene_list.append(gene)
    count_dict=Counter(gene_list)
    for gene in set(gene_list):
        rank_genelist.append([gene,count_dict[gene]])

    pre_rank_second_louvain_cluster_genelist[com]=rank_genelist

pre_rank_second_louvain_cluster_genelist_df={}
for com,genelist in pre_rank_second_louvain_cluster_genelist.items():
    rank_df=pd.DataFrame(columns=('gene', 'rank'))
    for n in range(len(genelist)):
        rank_df.loc[n] = [genelist[n][0], float(genelist[n][1])]
    rank_df=rank_df.sort_values(by=['rank'],ascending=False)
    pre_rank_second_louvain_cluster_genelist_df[com]=rank_df

pre_rank_third_louvain_cluster_genelist={}
for com,expset in ee_third_louvain_iteration_exposures.items():
    gene_list=[]
    rank_genelist=[]
    G_sub = nx.subgraph(backbone_ss_exposure_network,expset)   #let's create the subgraph corresponding to the cluster
    for edge in G_sub.edges():
        overlap_geneset=chem_gene_cleaned[edge[0]] & chem_gene_cleaned[edge[1]]
        for gene in overlap_geneset:
            gene_list.append(gene)
    count_dict=Counter(gene_list)
    for gene in set(gene_list):
        rank_genelist.append([gene,count_dict[gene]])

    pre_rank_third_louvain_cluster_genelist[com]=rank_genelist

pre_rank_third_louvain_cluster_genelist_df={}
for com,genelist in pre_rank_third_louvain_cluster_genelist.items():
    rank_df=pd.DataFrame(columns=('gene', 'rank'))
    for n in range(len(genelist)):
        rank_df.loc[n] = [genelist[n][0], float(genelist[n][1])]
    rank_df=rank_df.sort_values(by=['rank'],ascending=False)
    pre_rank_third_louvain_cluster_genelist_df[com]=rank_df


print(len(pre_rank_first_louvain_cluster_genelist_df[391]))


#First louvain interation
#Now we will run the traditional enrichment analysis for (Fisher's exact test + BH correction for those communties
#that perturb less than 20 genes, otherwise we will collect the leading genes from the GSEA with ranking
#and then performing again a traditional enrichment analysis

lead_genes_kegg_first_louvain_cluster_significant_fdr={}
lead_genes_gobp_first_louvain_cluster_significant_fdr={}
lead_genes_gocc_first_louvain_cluster_significant_fdr={}
lead_genes_gomf_first_louvain_cluster_significant_fdr={}
for com in list(pre_rank_first_louvain_cluster_genelist_df.keys()):
    num_of_genes=len(pre_rank_first_louvain_cluster_genelist_df[com])
    if num_of_genes>20:
        try:
            lead_genelist =list(enr_ranking_leading_genes(pre_rank_first_louvain_cluster_genelist_df[com]))
            come_enr_df_dict=enr_df_dict(lead_genelist)
            come_enr_fdr_df_dict=enr_fdr_list(come_enr_df_dict)
            lead_genes_kegg_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['KEGG']]
            lead_genes_gobp_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOBP']]
            lead_genes_gomf_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOMF']]
            lead_genes_gocc_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOCC']]

        except:
            lead_genes_kegg_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gobp_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gomf_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gocc_first_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']

    else:
        com_genelist=pre_rank_first_louvain_cluster_genelist_df[com]['gene'].tolist()
        try:
            come_enr_df_dict=enr_df_dict(com_genelist)
            come_enr_fdr_df_dict=enr_fdr_list(come_enr_df_dict)
            lead_genes_kegg_first_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['KEGG']]
            lead_genes_gobp_first_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOBP']]
            lead_genes_gomf_first_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOMF']]
            lead_genes_gocc_first_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOCC']]
        except:
            lead_genes_kegg_first_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gobp_first_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gomf_first_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gocc_first_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
        print(com,len(pre_rank_first_louvain_cluster_genelist_df[com]),'not ranked')


#second louvain interation
#Now we will run the traditional enrichment analysis for (Fisher's exact test + BH correction for those communties
#that perturb less than 20 genes, otherwise we will collect the leading genes from the GSEA with ranking
#and then performing again a traditional enrichment analysis

lead_genes_kegg_second_louvain_cluster_significant_fdr={}
lead_genes_gobp_second_louvain_cluster_significant_fdr={}
lead_genes_gocc_second_louvain_cluster_significant_fdr={}
lead_genes_gomf_second_louvain_cluster_significant_fdr={}
for com in list(pre_rank_second_louvain_cluster_genelist_df.keys()):
    num_of_genes=len(pre_rank_second_louvain_cluster_genelist_df[com])
    if num_of_genes>20:
        try:
            lead_genelist =list(enr_ranking_leading_genes(pre_rank_second_louvain_cluster_genelist_df[com]))
            come_enr_df_dict=enr_df_dict(lead_genelist)
            come_enr_fdr_df_dict=enr_fdr_list(come_enr_df_dict)
            lead_genes_kegg_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['KEGG']]
            lead_genes_gobp_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOBP']]
            lead_genes_gomf_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOMF']]
            lead_genes_gocc_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOCC']]

        except:
            lead_genes_kegg_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gobp_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gomf_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gocc_second_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']

    else:
        com_genelist=pre_rank_second_louvain_cluster_genelist_df[com]['gene'].tolist()
        try:
            come_enr_df_dict=enr_df_dict(com_genelist)
            come_enr_fdr_df_dict=enr_fdr_list(come_enr_df_dict)
            lead_genes_kegg_second_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['KEGG']]
            lead_genes_gobp_second_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOBP']]
            lead_genes_gomf_second_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOMF']]
            lead_genes_gocc_second_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOCC']]
        except:
            lead_genes_kegg_second_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gobp_second_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gomf_second_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gocc_second_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
        print(com,len(pre_rank_second_louvain_cluster_genelist_df[com]),'not ranked')


#third louvain interation
#Now we will run the traditional enrichment analysis for (Fisher's exact test + BH correction for those communties
#that perturb less than 20 genes, otherwise we will collect the leading genes from the GSEA with ranking
#and then performing again a traditional enrichment analysis

lead_genes_kegg_third_louvain_cluster_significant_fdr={}
lead_genes_gobp_third_louvain_cluster_significant_fdr={}
lead_genes_gocc_third_louvain_cluster_significant_fdr={}
lead_genes_gomf_third_louvain_cluster_significant_fdr={}
for com in list(pre_rank_third_louvain_cluster_genelist_df.keys()):
    num_of_genes=len(pre_rank_third_louvain_cluster_genelist_df[com])
    if num_of_genes>20:
        try:
            lead_genelist =list(enr_ranking_leading_genes(pre_rank_third_louvain_cluster_genelist_df[com]))
            come_enr_df_dict=enr_df_dict(lead_genelist)
            come_enr_fdr_df_dict=enr_fdr_list(come_enr_df_dict)
            lead_genes_kegg_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['KEGG']]
            lead_genes_gobp_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOBP']]
            lead_genes_gomf_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOMF']]
            lead_genes_gocc_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),come_enr_fdr_df_dict['GOCC']]

        except:
            lead_genes_kegg_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gobp_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gomf_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']
            lead_genes_gocc_third_louvain_cluster_significant_fdr[com]=[len(lead_genelist),'no enrichment with leading genes']

    else:
        com_genelist=pre_rank_third_louvain_cluster_genelist_df[com]['gene'].tolist()
        try:
            come_enr_df_dict=enr_df_dict(com_genelist)
            come_enr_fdr_df_dict=enr_fdr_list(come_enr_df_dict)
            lead_genes_kegg_third_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['KEGG']]
            lead_genes_gobp_third_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOBP']]
            lead_genes_gomf_third_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOMF']]
            lead_genes_gocc_third_louvain_cluster_significant_fdr[com]=[num_of_genes,come_enr_fdr_df_dict['GOCC']]
        except:
            lead_genes_kegg_third_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gobp_third_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gomf_third_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
            lead_genes_gocc_third_louvain_cluster_significant_fdr[com]=[num_of_genes,'not found']
        print(com,len(pre_rank_third_louvain_cluster_genelist_df[com]),'not ranked')


with open('output/Communities/Louvain/lead_genes_kegg_third_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_kegg_third_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gobp_third_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gobp_third_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gomf_third_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gomf_third_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gocc_third_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gocc_third_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)



with open('output/Communities/Louvain/lead_genes_kegg_second_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_kegg_second_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gobp_second_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gobp_second_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gomf_second_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gomf_second_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gocc_second_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gocc_second_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)


with open('output/Communities/Louvain/lead_genes_kegg_first_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_kegg_first_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gobp_first_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gobp_first_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gomf_first_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gomf_first_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
with open('output/Communities/Louvain/lead_genes_gocc_first_louvain_cluster_significant_fdr.pickle', 'wb') as handle:
    pk.dump(lead_genes_gocc_first_louvain_cluster_significant_fdr, handle, protocol=pk.HIGHEST_PROTOCOL)
