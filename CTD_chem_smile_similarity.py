import pandas as pd
import networkx as nx
import numpy as np
import pickle as pk
import itertools as it
from rdkit import Chem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

def calculate_tanimoto_coefficient(smiles1, smiles2):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    if mol1 is not None and mol2 is not None:
        fp1 = Chem.RDKFingerprint(mol1)
        fp2 = Chem.RDKFingerprint(mol2)
        tanimoto_coefficient = TanimotoSimilarity(fp1, fp2)
        return tanimoto_coefficient
    return None

#Let's import the chem_smile_dictionary
with open('output/chem_id_smile_dict.pickle', 'rb') as handle:
    chem_id_smile_dict = pk.load(handle)

unique_chem_id_smile_dict={}
for chem_id,smile_list in chem_id_smile_dict.items():
    if len(smile_list)==1:
        unique_chem_id_smile_dict[chem_id]=smile_list[0]
chem_id_list=list(unique_chem_id_smile_dict.keys())
chem_id_pairwise=list(it.combinations(chem_id_list,2))

chem_tanimoto_dict={}
for chem_id_pair in chem_id_pairwise:
    try:
        taninomoto_coeff=calculate_tanimoto_coefficient(unique_chem_id_smile_dict[chem_id_pair[0]],unique_chem_id_smile_dict[chem_id_pair[1]])
        chem_tanimoto_dict[chem_id_pair]=taninomoto_coeff
    except:
        pass

with open('output/chem_tanimoto_dict.pickle', 'wb') as handle:
    pk.dump(chem_tanimoto_dict, handle, protocol=pk.HIGHEST_PROTOCOL)
