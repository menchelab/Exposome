import os
import sys
import numpy as np
import networkx as nx
import itertools as it
import random as rd
import colorsys
import umap
import pickle as pk
import pymysql as mysql
from scipy.interpolate import interpn
from scipy.stats import gaussian_kde
from prettytable import PrettyTable
#from fisher import pvalue
import scipy.stats as st
import scipy.cluster.hierarchy as sch
import scipy.spatial.distance as dist
from sklearn.preprocessing import normalize
from scipy.cluster.hierarchy import fcluster
from sklearn import (manifold, datasets, decomposition, ensemble,
                     discriminant_analysis, random_projection,cluster)
import os.path
import pandas as pd
from collections import (defaultdict,Counter)
import time
# import statsmodels.sandbox.stats.multicomp as mc
import matplotlib.pyplot as plt
# %matplotlib inline
import pickle


ppi = pd.read_csv("input/PPI/autocore_symbol_lcc.csv"delimiter= ',',
           skipinitialspace=True)

G_ppi = nx.from_pandas_edgelist(ppi, 'symbol1', 'symbol2')

G_ppi_lcc = G_ppi.subgraph(max(nx.connected_components(G_ppi), key=len))  # extract lcc graph
print(G_ppi_lcc.number_of_nodes())
print(G_ppi_lcc.number_of_edges())
spl={}

pairwise=it.combinations(G_ppi_lcc.nodes(), 2)
for pair in pairwise:
    spl[pair]= nx.shortest_path_length(G_ppi_lcc, pair[0], pair[1])

with open('ppi_spl.pickle', 'wb') as handle:
    pickle.dump(spl, handle, protocol=pickle.HIGHEST_PROTOCOL)
