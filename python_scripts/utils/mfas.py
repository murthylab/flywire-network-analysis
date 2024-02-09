import numpy as np
import graph_tool.all as gt
from scipy.sparse.linalg import spsolve

def get_sorted_label(u, msfa):
    
    fu = u.copy()
    mfas_filter = fu.new_edge_property('bool')
    mfas_filter.a = msfa
    fu.set_edge_filter(mfas_filter)
    
    sort = gt.topological_sort(fu)

    return sort

def get_sorted_matrix(u, msfa):
    
    fu = u.copy()
    mfas_filter = fu.new_edge_property('bool')
    mfas_filter.a = msfa
    fu.set_edge_filter(mfas_filter)
    
    sort = gt.topological_sort(fu)
    A = gt.adjacency(u, weight=u.ep['#synapses'])[sort][:,sort].todense()

    return A


def mfas_naive_mc(u, n=100, mode="add", weighted=True):
    
    best_mfas = None
    best_len = 0
    best_mfas_size = u.ep["#synapses"].a.sum()
    info = {}
    
    info["min_size_curve"] = []
    
    for _ in range(n):
        fu = u.copy()
        mfas = fu.new_edge_property('bool')
        mfas.a = False if mode == "add" else True
        fu.set_edge_filter(mfas)

        if mode == "add":
            idxs = np.random.permutation(np.nonzero(1-mfas.a)[0])
            for idx in idxs:
                mfas.a[idx] = 1
                if not gt.is_DAG(fu):
                    mfas.a[idx] = 0
                
        elif mode == "sub":
            idxs = np.random.permutation(np.nonzero(mfas.a)[0])
            for idx in idxs:
                mfas.a[idx] = 0
                if gt.is_DAG(fu):
                    break
                
        mfas_size = fu.ep["#synapses"].a[mfas.a == 0].sum()
        
        if weighted:
            if mfas_size < best_mfas_size:
                best_mfas = mfas.a
                best_mfas_size = mfas_size
        else:
            if best_len < mfas.a.sum():
                best_mfas = mfas.a
                best_mfas_size = mfas_size
                best_len = mfas.a.sum()
        
        info["min_size_curve"].append(best_mfas_size)

    return best_mfas, best_mfas_size, info


def mfas_upperbound(u, mode="diff", weighted=True):
    
    best_mfas = None
    best_mfas_size = u.ep["#synapses"].a.sum()
    info = {}
    
    fu = u.copy()
    sub = fu.new_vertex_property('bool')
    sub.a = True
    fu.set_vertex_filter(sub)

    degree_rank = fu.new_vertex_property('int32_t')
    rank_l, rank_r = 0, u.num_vertices()
    
    if weighted:
        weight = fu.ep["#synapses"]
    else:
        weight = None
    
    while sub.a.sum() > 0:
         
        while True:
            allchecked = True
            for v in fu.vertices():
                if v.out_degree(weight=weight) == 0:
                    degree_rank[v] = rank_r
                    rank_r -= 1
                    sub[v] = False
            if allchecked:
                break
                
        while True:
            allchecked = True
            for v in fu.vertices():
                if v.in_degree(weight=weight) == 0:
                    allchecked = False
                    degree_rank[v] = rank_l
                    rank_l += 1
                    sub[v] = False
            if allchecked:
                break
                
        if mode == "diff":
            delta = [v.out_degree(weight=weight) 
                      - v.in_degree(weight=weight) for v in fu.vertices()]
        elif mode == "ratio":
            delta = [v.out_degree(weight=weight) 
                      / v.in_degree(weight=weight) for v in fu.vertices()]
            
        if len(delta) > 0:
            v = list(fu.vertices())[np.argmax(delta)]
            degree_rank[v] = rank_l
            rank_l += 1
            sub[v] = False
    
    fu.clear_filters()
    mfas = fu.new_edge_property('bool')
    mfas.a = True
    
    for e in fu.edges():
        if degree_rank[e.source()] >= degree_rank[e.target()]:
            mfas[e] = False
        
    fu.set_edge_filter(mfas)

    mfas_size = fu.ep["#synapses"].a[mfas.a == 0].sum()

    if mfas_size < best_mfas_size:
        best_mfas = mfas.a
        best_mfas_size = mfas_size

    return best_mfas, best_mfas_size, info


def mfas_quadratic(u, weighted=True):
    
    ## min Σ Aij * (ri - rj - 1)^2
    ## -> Lr = din - dout
    
    best_mfas = None
    best_mfas_size = u.ep["#synapses"].a.sum()
    info = {}
    
    fu = u.copy()
    
    wcc, _ =  gt.label_components(fu, directed=False)
    comp_labels = np.unique(wcc.a)
    wcc_filters = []
    for c in comp_labels:
        c_filter = fu.new_vertex_property("bool")
        for v in fu.vertices():
            c_filter[v] = True if wcc[v] == c else False
        wcc_filters.append(c_filter)
        
        
    if weighted:
        weight = fu.ep["#synapses"]
    else:
        weight = None

    node_rank = fu.new_vertex_property('float')
    
    for wcc_f in wcc_filters:
        fu.set_vertex_filter(wcc_f)
        mask = (wcc_f.a == 1)
        degree_diff = fu.new_vertex_property('int32_t')
        for v in fu.vertices():
            degree_diff[v] = v.out_degree(weight=weight) - v.in_degree(weight=weight)
        A = gt.adjacency(fu, weight=weight)
        L = gt.laplacian(fu, weight=weight) - A.T
        if fu.num_vertices() > 1:
            wcc_rank = spsolve(L, degree_diff.a[mask])
            node_rank.a[mask] = wcc_rank + node_rank.a.max() - wcc_rank.min() + 1
        else:
            wcc_rank = 0
            node_rank.a[mask] = node_rank.a.max() + 1
        fu.clear_filters()

    mfas = fu.new_edge_property('bool')
    mfas.a = True
    
    info["r"] = node_rank.a
    
    for e in fu.edges():
        if node_rank[e.source()] <= node_rank[e.target()]:
            mfas[e] = False
        
    fu.set_edge_filter(mfas)

    mfas_size = fu.ep["#synapses"].a[mfas.a == 0].sum()

    if mfas_size < best_mfas_size:
        best_mfas = mfas.a
        best_mfas_size = mfas_size

    return best_mfas, best_mfas_size, info