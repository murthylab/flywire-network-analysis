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
	parser.add_argument('--neuropil_id', type=int, required=True, help='id')
	args = parser.parse_args()
# 	neuropil_names = ["FB","EB","PB","NO","AMMC_L","AMMC_R","FLA_L","FLA_R","CAN_L","CAN_R","PRW","SAD","GNG","AL_L","AL_R","LH_L","LH_R","MB_CA_L","MB_CA_R","MB_PED_L","MB_PED_R","MB_VL_L","MB_VL_R","MB_ML_L","MB_ML_R","BU_L","BU_R","GA_L","GA_R","LAL_L","LAL_R","SLP_L","SLP_R","SIP_L","SIP_R","SMP_L","SMP_R","CRE_L","CRE_R","SCL_L","SCL_R","ICL_L","ICL_R","IB_L","IB_R","ATL_L","ATL_R","VES_L","VES_R","EPA_L","EPA_R","GOR_L","GOR_R","SPS_L","SPS_R","IPS_L","IPS_R","AOTU_L","AOTU_R","AVLP_L","AVLP_R","PVLP_L","PVLP_R","PLP_L","PLP_R","WED_L","WED_R","ME_L","ME_R","AME_L","AME_R","LO_L","LO_R","LOP_L","LOP_R","LA_L","LA_R","OCG"]
	neuropil_names = ["PB","NO","AMMC_L","FLA_R"]
# 	neuropil_names = [ 'AME_R']
	# construct graph

	ver = 630
	# fdir = "~/seungmount/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_{}.feather".format(ver, ver)
	fdir = "/mnt/cup/labs/seung/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_ntavg_{}.feather".format(ver, ver)
	syn_table = pd.read_feather(fdir)

	neuropil_name = neuropil_names[args.neuropil_id]
	syn_table = syn_table[syn_table["neuropil"] == neuropil_name] 
    
	# build index dictionary
	cellids =  np.unique(syn_table[["pre_pt_root_id", "post_pt_root_id"]])
	nid2cid = {i: cid for i, cid in enumerate(cellids)}
	cid2nid = {cid: i for i, cid in enumerate(cellids)}
    
	AL_motif = None
	with open('/mnt/cup/labs/seung/research/runzhey/FlyWireAnalysis/saved/neuropil-motifs-314/AL_L-v314-motif.pickle', 'rb') as handle:
		AL_motif = pickle.load(handle)

	def_motif_list = AL_motif['AL_L'][0]

	syn_table["pre_nid"] = pd.Series([cid2nid[cid] for cid in syn_table["pre_pt_root_id"]], 
		 index=syn_table.index)
	syn_table["post_nid"] = pd.Series([cid2nid[cid] for cid in syn_table["post_pt_root_id"]], 
		 index=syn_table.index)
    
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

	np_motif_profiles = {}
	np_cfg_motif_profiles = {}

	strong_filter = g.new_ep('bool')
	for e in g.edges():
		strong_filter[e] = True if g.ep["syn_count"][e] >=5 else False

	g.set_edge_filter(strong_filter)
	u = gt.Graph(g, prune=True)

	# get observed motif counts
	np_motif_profiles[neuropil_name] = gt.motifs(u, k=3, p=1.0, motif_list=def_motif_list, return_maps=False)

	np_cfg_motif_profiles[neuropil_name] = []
	u_cfg = gt.Graph(u, prune=True)
	for _ in range(100):
		gt.random_rewire(u_cfg, model='configuration')
		np_cfg_motif_profiles[neuropil_name].append(gt.motifs(u_cfg, 
			k=3, p=1.0, motif_list=def_motif_list, return_maps=False)[1])

	with open('saved/neuropil-motifs-{}/{}-v{}-motif.pickle'.format(ver, neuropil_name, ver), 'wb') as handle:
		pickle.dump(np_motif_profiles, handle, protocol=pickle.HIGHEST_PROTOCOL)

	with open('saved/neuropil-motifs-{}/{}-v{}-cfg-motif.pickle'.format(ver, neuropil_name, ver), 'wb') as handle:
		pickle.dump(np_cfg_motif_profiles, handle, protocol=pickle.HIGHEST_PROTOCOL)






