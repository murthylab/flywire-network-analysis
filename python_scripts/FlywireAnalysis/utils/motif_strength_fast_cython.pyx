import math
import networkx as nx
from collections import defaultdict
from cython cimport boundscheck, wraparound

# Cython-specific annotations
cimport cython
from libc.math cimport sqrt

# Define the Triplet class
cdef class Triplet:
    cdef public dict vertices
    cdef public set edges

    def __init__(self, vertices, edges=None):
        vertices.sort()
        self.vertices = {k: v for k, v in enumerate(vertices)}
        self.edges = set()
        if edges is not None:
            for e in edges:
                self.add_edge(*e)

    @property
    def ids(self):
        return {v: k for k, v in self.vertices.items()}

    def get_id(self, v):
        return self.ids[v]

    def get_edges(self):
        return [(self.vertices[e[0]], self.vertices[e[1]]) for e in self.edges]

    def add_edge(self, u, v):
        e = (self.get_id(u), self.get_id(v))
        self.edges.add(e)

    def __len__(self):
        return len(self.edges)

    def __repr__(self):
        return str(self.vertices) + ' ' + str(sorted(self.edges))

    def __hash__(self):
        return hash(self.__repr__())

    def __eq__(self, X):
        return self.vertices == X.vertices and self.edges == X.edges

# Define the OnlineStatistics class
cdef class OnlineStatistics:
    cdef public int count
    cdef public float mean
    cdef public float M2

    def __init__(self):
        self.count = 0
        self.mean = 0.0
        self.M2 = 0.0

    cpdef void update(self, float value):
        self.count += 1
        delta = value - self.mean
        self.mean += delta / self.count
        self.M2 += delta * (value - self.mean)

    cpdef float get_mean(self):
        if self.count < 1:
            return float('nan')
        return self.mean

    cpdef float get_variance(self):
        if self.count < 2:
            return float('nan')
        else:
            return self.M2 / self.count

    cpdef float get_standard_deviation(self):
        return sqrt(self.get_variance())

@cython.boundscheck(False)
@cython.wraparound(False)
def collect_three_neuron_motifs_weighted(V, E_w, u, motifs):
    cdef dict matches = {k: OnlineStatistics() for k in motifs.keys()}
    cdef set tri = set()
    cdef set E = set(E_w.keys())
    cdef list E_list = list(E)
    cdef int a, b, c
    cdef float t_sum_weight, t_avg_weight
    cdef int t_count_weight
    cdef set predecessors, successors
    cdef list candi_c
    cdef Triplet t

    for a, b in E_list:
        t_sum_weight = E_w[(a, b)]
        t_count_weight = 1
        predecessors = set(u.predecessors(b))
        successors = set(u.successors(b))
        candi_c = list(predecessors.union(successors))

        for c in candi_c:
            if c != a and c != b:
                t = Triplet([a, b, c])
                t.add_edge(a, b)

                if (a, c) in E:
                    t_sum_weight += E_w[(a, c)]
                    t_count_weight += 1
                    t.add_edge(a, c)

                if (c, a) in E:
                    t_sum_weight += E_w[(c, a)]
                    t_count_weight += 1
                    t.add_edge(c, a)

                if (b, c) in E:
                    t_sum_weight += E_w[(b, c)]
                    t_count_weight += 1
                    t.add_edge(b, c)

                if (c, b) in E_w.keys():
                    t_sum_weight += E_w[(c, b)]
                    t_count_weight += 1
                    t.add_edge(c, b)

                if (b, a) in E_w.keys():
                    t_sum_weight += E_w[(b, a)]
                    t_count_weight += 1
                    t.add_edge(b, a)

                t_avg_weight = t_sum_weight / t_count_weight

                if t_count_weight >= 2:
                    for k, m in motifs.items():
                        if match(t, m):
                            matches[k].update(t_avg_weight)

    return matches

@cython.boundscheck(False)
@cython.wraparound(False)
def count_three_neuron_motifs_weighted(V, E, u, motifs):
    cdef dict tri_motifs_stats = collect_three_neuron_motifs_weighted(V, E, u, motifs)
    cdef dict tri_mean_weight = {k: v.get_mean() for k, v in tri_motifs_stats.items() if v.get_mean() > 0}
    cdef dict tri_std_weight = {k: v.get_standard_deviation() for k, v in tri_motifs_stats.items() if v.get_standard_deviation() > 0}

    return tri_mean_weight, tri_std_weight

def match_edges(A, B, reflections=range(3)):
    if A == B:
        return True
    if len(reflections) > 1:
        for k in reflections:
            if match_edges(reflect(A, k), B, reflections=[r for r in reflections if r != k]):
                return True
    return False

def match(A, B):
    return match_edges(A.edges, B.edges)

def reflect(edges, identity):
    reflect = [x for x in range(3) if x != identity]
    remap = {identity: identity, reflect[0]: reflect[1], reflect[1]: reflect[0]}
    return set([(remap[u], remap[v]) for u, v in edges])