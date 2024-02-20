import os
import random
from time import time
import pandas as pd
import numpy as np
import networkx as nx
import argparse
import pickle

from utils.utils import *
from utils.graph_creation import *
from utils.motif_counts import *
from utils.visualization import *
import graph_tool.all as gt

sns.set_context("talk")

if __name__ == "__main__":

	parser = argparse.ArgumentParser(description='Counting motifs in randomized graph')
	parser.add_argument('--exp_id', type=int, required=True, help='id')
	args = parser.parse_args()
	# construct graph

	ver = 630
	# fdir = "~/seungmount/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_{}.feather".format(ver, ver)
	fdir = "/mnt/cup/labs/seung/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_ntavg_{}.feather".format(ver, ver)
	syn_table = pd.read_feather(fdir)

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

	# get neuropil subgraph
	np_motif_profiles = {}
	np_cfg_motif_profiles = {}
	
	neuropil_name = 'whole'
    
	strong_filter = g.new_ep('bool')
	for e in g.edges():
		strong_filter[e] = True if g.ep["syn_count"][e] >=5 else False
	g.set_edge_filter(strong_filter)
	u = gt.Graph(g, prune=True)
	g.clear_filters()
    
	if args.exp_id < 1:
		# get observed motif counts
		np_motif_profiles[neuropil_name] = gt.motifs(u, k=3, p=1.0, 
			motif_list=None, return_maps=False)

		with open('saved/neuropil-motifs-630/{}-v630-motif.pickle'.format(neuropil_name), 'wb') as handle:
			pickle.dump(np_motif_profiles, handle, protocol=pickle.HIGHEST_PROTOCOL)

	np_cfg_motif_profiles[neuropil_name] = []
	u_cfg = gt.Graph(u, prune=True)
    
	for _ in range(5):
		gt.random_rewire(u_cfg, model='configuration')
		np_cfg_motif_profiles[neuropil_name].append(gt.motifs(u_cfg, 
			k=3, p=1.0, motif_list=None, return_maps=False)[1])

	with open('saved/neuropil-motifs-630/{}-v630-cfg-motif-{}.pickle'.format(neuropil_name, args.exp_id), 'wb') as handle:
		pickle.dump(np_cfg_motif_profiles, handle, protocol=pickle.HIGHEST_PROTOCOL)

