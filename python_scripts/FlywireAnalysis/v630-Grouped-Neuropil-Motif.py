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
    neuropil_names = [ 'AL_L', 'AL_R', 'AME_L', 'AMMC_L', 'AMMC_R', 'AOTU_L', 
                       'AOTU_R', 'ATL_L', 'ATL_R', 'AVLP_L', 'AVLP_R', 'BU_L', 
                       'BU_R', 'CAN_L', 'CAN_R', 'CRE_L', 'CRE_R', 'EB', 'EPA_L', 
                       'EPA_R', 'FB', 'FLA_L', 'FLA_R', 'GA_L', 'GA_R', 'GNG', 'GOR_L', 
                       'GOR_R', 'IB_L', 'IB_R', 'ICL_L', 'ICL_R', 'IPS_L', 'IPS_R', 
                       'LAL_L', 'LAL_R', 'LH_L', 'LH_R', 'LOP_L', 'LOP_R', 'LO_L', 
                       'LO_R', 'MB_CA_L', 'MB_CA_R', 'MB_ML_L', 'MB_ML_R', 'MB_PED_L', 
                       'MB_PED_R', 'MB_VL_L', 'MB_VL_R', 'ME_L', 'ME_R', 'NO', 'PB', 
                       'PLP_L', 'PLP_R', 'PRW', 'PVLP_L', 'PVLP_R', 'SAD', 'SCL_L', 
                       'SCL_R', 'SIP_L', 'SIP_R', 'SLP_L', 'SLP_R', 'SMP_L', 'SMP_R', 
                       'SPS_L', 'SPS_R', 'VES_L', 'VES_R', 'WED_L', 'WED_R']
    
    grouped_neuropil_names = {
        'INP_L': ['ATL_L','CRE_L','IB_L','SCL_L','ICL_L'],
        'INP_R': ['ATL_R','CRE_R','IB_R','SCL_R','ICL_R'],
        'LX_L': ['BU_L','LAL_L', 'GA_L'],
        'LX_R': ['BU_R','LAL_R', 'GA_R'],
        'VMNP_L': ['VES_L','EPA_L','GOR_L','SPS_L','IPS_L'],
        'VMNP_R': ['VES_R','EPA_R','GOR_R','SPS_R','IPS_R'],
        'PENP': ['PRW', 'CAN_L', 'CAN_R', 'FLA_L', 'FLA_R'],
        'MB_L': ['MB_CA_L','MB_ML_L', 'MB_PED_L', 'MB_VL_L'],
        'MB_R': ['MB_CA_R','MB_ML_R', 'MB_PED_R', 'MB_VL_R']
    }
    
    flatten = lambda x: [item for sublist in x for item in sublist]
    covered_ = flatten([grouped_neuropil_names[k] for k in grouped_neuropil_names.keys()])
    # after grouping there should be 44 grouped neuropils
    for np_name in neuropil_names:
        if not np_name in covered_:
            grouped_neuropil_names[np_name] = [np_name]
    
    # construct graph

    ver = 447
    # fdir = "~/seungmount/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_{}.feather".format(ver, ver)
    fdir = "/mnt/cup/labs/seung/research/runzhey/datasets/FlyWire/{}/syn_proof_analysis_{}.feather".format(ver, ver)
    syn_table = pd.read_feather(fdir)

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
    survival = pd.read_csv("/mnt/cup/labs/seung/research/runzhey/FlyWireAnalysis/saved/survival_447.csv")
    np_motif_profiles = {}
    np_cfg_motif_profiles = {}
    
    neuropil_name = neuropil_names[args.neuropil_id]
    neuropil_cids = list(survival[(survival["max_in_region"]==neuropil_name) & (survival["max_out_region"]==neuropil_name)]['root_id'])

    neuropil_filter = g.new_vp('bool')
    for v in g.vertices():
        neuropil_filter[v] = True if g.vp["cellid"][v] in neuropil_cids else False

    g.set_vertex_filter(neuropil_filter)
    u = gt.Graph(g, prune=True)
    g.clear_filters()

    strong_filter = u.new_ep('bool')
    for e in u.edges():
        strong_filter[e] = True if u.ep["syn_count"][e] >=5 else False

    u.set_edge_filter(strong_filter)
    u = gt.Graph(u, prune=True)

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
