import pandas as pd
import networkx as nx
import numpy as np
import pickle as pk
import itertools as it
import rdkit
from pubchempy import *
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit import Chem
from rdkit.Chem import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

def calculate_tanimoto_coefficient(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is not None and mol2 is not None:
        fp1 = FingerprintMols.FingerprintMol(mol1)  # Get the fingerprint of the first molecule
        fp2 = FingerprintMols.FingerprintMol(mol2)  # Get the fingerprint of the second molecule
        tanimoto_coefficient = DataStructs.TanimotoSimilarity(fp1, fp2)
        return tanimoto_coefficient
    return None

#Let's import the chemical-gene interactions from CTD (downloaded on 5th April 2021)
chem_gene_df = pd.read_csv("input/CTD/CTD_chem_gene_ixns.tsv",delimiter= '\t', skipinitialspace=True)
#Here, we filter for only the interactions that regards the H. Sapiens
chem_homo = chem_gene_df[(chem_gene_df['Organism'] == 'Homo sapiens')]

#Let's import the Chemical Vocabolary
chem_vocabolary = pd.read_csv("input/CTD/CTD_chemicals_cleaned.txt",delimiter= '\t',
           skipinitialspace=True)

#Let's convert the chemical names from CTD to smiles
chem_name_set=set(chem_homo['# ChemicalName'])
chem_smile_list=[]
chem_name_smile_dict={}
for chem in chem_name_set:
    try:
        for compound in get_compounds(chem, 'name'):
            smile=compound.isomeric_smiles
            chem_smile_list.append(smile)
            chem_name_smile_dict[smile]=chem
    except:
        try:
            syn_list = chem_vocabolary[chem_vocabolary['# ChemicalName']==chem]['Synonyms'].values[0].split("|")
            for syn in syn_list:
                try:
                    for compound in get_compounds(syn, 'name'):
                        smile=compound.isomeric_smiles
                        chem_smile_list.append(smile)
                        chem_name_smile_dict[smile]=chem
                except:
                    pass
        except:
            pass
with open('output/chem_name_smile_dict.pickle', 'wb') as handle:
    pk.dump(chem_name_smile_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

#Let's define a dictionary that will convert the chemical name in chemical ID
chem_name_id_conversion = {}
for i,v in chem_gene_df.iterrows():
        chem_name_id_conversion[v["# ChemicalName"]] = v["ChemicalID"]

#Let's invert the chem_name_smile_dict
def invert_dict(d):
    inverted_dict = {}
    for key, value in d.items():
        if value not in inverted_dict:
            inverted_dict[value] = [key]
        else:
            inverted_dict[value].append(key)
    return inverted_dict

chem_name_smile_dict_inv=invert_dict(chem_name_smile_dict)

#Let's convert the chem name into chem ids
chem_id_smile_dict={}
for chem_name,smile_list in chem_name_smile_dict_inv.items():
    chem_id_smile_dict[chem_name_id_conversion[chem_name]]=smile_list

with open('output/chem_id_smile_dict.pickle', 'wb') as handle:
    pk.dump(chem_id_smile_dict, handle, protocol=pk.HIGHEST_PROTOCOL)

chem_id_list=list(chem_id_smile_dict.keys())
chem_id_pairwise=list(it.combinations(chem_id_list,2))

chem_tanimoto_dict={}
for chem_id_pair in chem_id_pairwise:
    try:
        taninomoto_coeff=calculate_tanimoto_coefficient(chem_id_smile_dict[chem_id_pair[0]][0],chem_id_smile_dict[chem_id_pair[1]][0])
        chem_tanimoto_dict[chem_id_pair]=taninomoto_coeff
    except:
        pass

with open('output/chem_tanimoto_dict.pickle', 'wb') as handle:
    pk.dump(chem_tanimoto_dict, handle, protocol=pk.HIGHEST_PROTOCOL)
