import os
import random
from time import time
import pandas as pd
import numpy as np
import networkx as nx
import argparse
import pickle
import math

from utils.utils import *
from utils.graph_creation import *
from utils.motif_counts import *
from utils.visualization import *
import graph_tool.all as gt

sns.set_context("talk")


from utils.motif_strength_parallel_cython import count_three_neuron_motifs_weighted


if __name__ == "__main__":

# 	parser = argparse.ArgumentParser(description='Counting motifs in randomized graph')
# 	parser.add_argument('--neuropil_id', type=int, required=True, help='id')
# 	args = parser.parse_args()
	neuropil_names = ["whole"]
	# construct graph

	ver = 630
	# fdir = "~/seungmount/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_{}.feather".format(ver, ver)
	fdir = "/mnt/cup/labs/seung/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_ntavg_{}.feather".format(ver, ver)
	syn_table = pd.read_feather(fdir)

	neuropil_name = "whole"
# 	syn_table = syn_table[syn_table["neuropil"] == neuropil_name] 
    
	# build index dictionary
	cellids =  np.unique(syn_table[["pre_pt_root_id", "post_pt_root_id"]])
	nid2cid = {i: cid for i, cid in enumerate(cellids)}
	cid2nid = {cid: i for i, cid in enumerate(cellids)}

	syn_table["pre_nid"] = pd.Series([cid2nid[cid] for cid in syn_table["pre_pt_root_id"]], 
		 index=syn_table.index)
	syn_table["post_nid"] = pd.Series([cid2nid[cid] for cid in syn_table["post_pt_root_id"]], 
		 index=syn_table.index)

	edge_list = [(e[0], e[1], e[2], e[3]) for e in 
		 syn_table[["pre_nid", "post_nid", "neuropil", "syn_count"]].values]

	merged_syn_table = syn_table[["pre_nid", "post_nid", "syn_count"]
		 ].groupby(by=["pre_nid", "post_nid"]).sum().reset_index()

	merged_edge_list = [(e[0], e[1], e[2]) for e in merged_syn_table.values]

	g = gt.Graph()
	g.add_vertex(len(cellids))
	e_syn_count = g.new_ep("int32_t")
	g.add_edge_list(merged_edge_list, eprops=[e_syn_count])
	g.vp["cellid"] = g.new_vp("int64_t")
	g.vp["cellid"].a = cellids
	g.ep["syn_count"] = e_syn_count

	strong_filter = g.new_ep('bool')
	for e in g.edges():
		strong_filter[e] = True if g.ep["syn_count"][e] >=5 else False

	g.set_edge_filter(strong_filter)
	u = gt.Graph(g, prune=True)

	u_nx = pyintergraph.gt2nx(u)
    
	g = gt.Graph()
	u = gt.Graph()

	edge_weight_dict = {} 
	for n1, n2, data in u_nx.edges(data=True):
		edge_weight_dict[(n1, n2)] = data['syn_count']
        
	motif_strength_mean, motif_strength_std = count_three_neuron_motifs_weighted(set(u_nx.nodes()), edge_weight_dict, u_nx, motifs, n_processes=20)

	with open('saved/neuropil-motifs-strength-{}/{}-v{}-motif-strength-mean.pickle'.format(ver, neuropil_name, ver), 'wb') as handle:
		pickle.dump(motif_strength_mean, handle, protocol=pickle.HIGHEST_PROTOCOL)

	with open('saved/neuropil-motifs-strength-{}/{}-v{}-motif-strength-std.pickle'.format(ver, neuropil_name, ver), 'wb') as handle:
		pickle.dump(motif_strength_std, handle, protocol=pickle.HIGHEST_PROTOCOL)







