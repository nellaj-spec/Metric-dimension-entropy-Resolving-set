"""Metric-dimension entropy (Algorithm 1 of the manuscript).

Membership of a vertex v with respect to a resolving set W:
    mu_W(v) = 1 / (1 + min_{u in W} d(v, u))
Probabilities p_W(v) = mu_W(v) / sum_x mu_W(x); entropy I_W(G) = -sum p log2 p.
The basis-averaged entropy I(G) (Definition 5.1) is the mean of I_W(G) over all
metric bases W of G; I_min and I_max are the minimum and maximum.
"""
import math
import networkx as nx
from metric_basis import all_metric_bases, distance_matrix


def entropy_for_set(D, W):
    """I_W(G) for a resolving set W given the distance matrix D (list of lists)."""
    n = len(D)
    mu = [1.0 / (1.0 + min(D[v][u] for u in W)) for v in range(n)]
    s = sum(mu)
    return -sum((m / s) * math.log2(m / s) for m in mu)


def metric_dimension_entropy(G, bases=None):
    """Return a dict with beta, number of bases, I, I_min, I_max and I_norm.

    If `bases` is given (e.g. a single basis from integer programming), the
    values are computed over that list only; the field 'n_bases' then reports
    the length of the list and 'averaged' is False.
    """
    G = nx.convert_node_labels_to_integers(G)
    n = G.number_of_nodes()
    D = distance_matrix(G)
    averaged = bases is None
    if bases is None:
        bases = all_metric_bases(D)
    vals = [entropy_for_set(D, W) for W in bases]
    I = sum(vals) / len(vals)
    return {
        "n": n,
        "m": G.number_of_edges(),
        "beta": len(bases[0]),
        "n_bases": len(bases),
        "I": I,
        "I_min": min(vals),
        "I_max": max(vals),
        "I_norm": I / math.log2(n) if n > 1 else float("nan"),
        "averaged": averaged,
    }


def eccentricity_entropy_for_set(D, W):
    """Entropy under the eccentricity-based efficiency of Section 2.3.

    A basis element u has efficiency ecc(u)/n; any other vertex takes the
    efficiency of the basis element nearest to it (the average if several
    basis elements are equally near).
    """
    n = len(D)
    ecc = [max(row) for row in D]
    mu = []
    for v in range(n):
        if v in W:
            mu.append(ecc[v] / n)
        else:
            dmin = min(D[v][u] for u in W)
            near = [u for u in W if D[v][u] == dmin]
            mu.append(sum(ecc[u] / n for u in near) / len(near))
    s = sum(mu)
    return -sum((m / s) * math.log2(m / s) for m in mu)
