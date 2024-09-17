
import pickle as pk
import networkx as nx
from pecanpy.graph import AdjlstGraph
from pecanpy import pecanpy as node2vec

filtered_weighted_exp_graph_significant = nx.read_weighted_edgelist("backbone_exp_graph_significant_weighted.edgelist")

filtered_unweighted_exp_graph_significant=nx.Graph()
for e in filtered_weighted_exp_graph_significant.edges():
    filtered_unweighted_exp_graph_significant.add_edge(*e)

g = AdjlstGraph()
for e in list(filtered_weighted_exp_graph_significant.edges()):
    weight = filtered_weighted_exp_graph_significant[e[0]][e[1]]["weight"]
    g.add_edge(e[0],e[1],weight)
g.save("filtered_weighted_exp_graph_significant_new.edg")
g = node2vec.PreComp(p=0.5, q=1, workers=20, verbose=True,extend=True)
g.read_edg("filtered_weighted_exp_graph_significant_new.edg", weighted=True, directed=False)
g.preprocess_transition_probs()

filtered_weighted_emd = g.embed(dim=64,num_walks=100,walk_length=50)
#Let's write the results
with open('filtered_weighted_emd.pickle', 'wb') as handle:
    pk.dump(filtered_weighted_emd, handle, protocol=pk.HIGHEST_PROTOCOL)

g = AdjlstGraph()
for e in list(filtered_unweighted_exp_graph_significant.edges()):
    g.add_edge(e[0],e[1])
g.save("filtered_unweighted_exp_graph_significant_new.edg")

g = node2vec.PreComp(p=0.5, q=1, workers=20, verbose=True,extend=True)
g.read_edg("filtered_unweighted_exp_graph_significant_new.edg", weighted=False, directed=False)
g.preprocess_transition_probs()
filtered_unweighted_emd = g.embed(dim=64,num_walks=100,walk_length=50)
#Let's write the results
with open('filtered_unweighted_emd.pickle', 'wb') as handle:
    pk.dump(filtered_unweighted_emd, handle, protocol=pk.HIGHEST_PROTOCOL)


unfiltered_weighted_exp_graph_significant = nx.read_weighted_edgelist("unfiltered_weighted_exp_graph_significant.edgelist")

unfiltered_unweighted_exp_graph_significant=nx.Graph()
for e in unfiltered_weighted_exp_graph_significant.edges():
    unfiltered_unweighted_exp_graph_significant.add_edge(*e)

g = AdjlstGraph()
for e in list(unfiltered_unweighted_exp_graph_significant.edges()):
    g.add_edge(e[0],e[1])
g.save("unfiltered_unweighted_exp_graph_significant_new.edg")
g = node2vec.PreComp(p=0.5, q=1, workers=20, verbose=False,extend=True)
g.read_edg("unfiltered_unweighted_exp_graph_significant_new.edg", weighted=False, directed=False)
g.preprocess_transition_probs()
unfiltered_unweighted_emd = g.embed(dim=64,num_walks=100,walk_length=50)
#Let's write the results
with open('unfiltered_unweighted_emd.pickle', 'wb') as handle:
    pk.dump(unfiltered_unweighted_emd, handle, protocol=pk.HIGHEST_PROTOCOL)

g = AdjlstGraph()
for e in list(unfiltered_weighted_exp_graph_significant.edges()):
    weight = unfiltered_weighted_exp_graph_significant[e[0]][e[1]]["weight"]
    g.add_edge(e[0],e[1],weight)
g.save("unfiltered_weighted_exp_graph_significant_new.edg")
g = node2vec.PreComp(p=0.5, q=1, workers=20, verbose=False,extend=True)
g.read_edg("unfiltered_weighted_exp_graph_significant_new.edg", weighted=True, directed=False)
g.preprocess_transition_probs()
unfiltered_weighted_emd = g.embed(dim=64,num_walks=100,walk_length=50)
#Let's write the results
with open('unfiltered_weighted_emd.pickle', 'wb') as handle:
    pk.dump(unfiltered_weighted_emd, handle, protocol=pk.HIGHEST_PROTOCOL)
