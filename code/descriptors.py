"""Comparison graph entropies and topological indices (Table 4 of the manuscript)."""
import math
from collections import Counter
import numpy as np
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from metric_basis import distance_matrix


def _H(ps):
    return -sum(p * math.log2(p) for p in ps if p > 0)


def automorphism_orbits(G):
    nodes = list(G.nodes())
    orbit = {v: {v} for v in nodes}
    for iso in GraphMatcher(G, G).isomorphisms_iter():
        for v, w in iso.items():
            orbit[v].add(w)
    return {frozenset(o) for o in orbit.values()}


def descriptors(G):
    G = nx.convert_node_labels_to_integers(G)
    n, m = G.number_of_nodes(), G.number_of_edges()
    D = np.array(distance_matrix(G))
    deg = dict(G.degree())
    s = D.sum(axis=1)                       # distance sums s_v
    W = D.sum() / 2                         # Wiener index
    out = {"diameter": int(D.max()), "radius": int(D.max(axis=1).min()), "wiener": W}
    out["M1"] = sum(d * d for d in deg.values())
    out["M2"] = sum(deg[u] * deg[v] for u, v in G.edges())
    out["randic"] = sum((deg[u] * deg[v]) ** -0.5 for u, v in G.edges())
    mu = m - n + 1                          # cyclomatic number
    out["balaban"] = m / (mu + 1) * sum((s[u] * s[v]) ** -0.5 for u, v in G.edges())
    out["platt"] = sum(deg[u] + deg[v] - 2 for u, v in G.edges())
    out["H_degree"] = _H([d / (2 * m) for d in deg.values()])
    out["H_distance"] = _H([sv / (2 * W) for sv in s])
    nk = Counter(D[np.triu_indices(n, 1)].tolist())
    out["H_bonchev"] = -sum(c * (k / W) * math.log2(k / W) for k, c in nk.items())
    out["H_orbit"] = _H([len(o) / n for o in automorphism_orbits(G)])
    L = nx.laplacian_matrix(G).toarray().astype(float) / (2 * m)
    lam = np.linalg.eigvalsh(L)
    out["H_vonneumann"] = _H([x for x in lam if x > 1e-12])
    return {k: v + 0.0 for k, v in out.items()}   # + 0.0 turns -0.0 into 0.0
