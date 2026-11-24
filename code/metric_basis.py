"""Finding and enumerating metric bases.

* distance_matrix  - all-pairs shortest paths by breadth-first search (Step 1)
* all_metric_bases - exhaustive enumeration of every minimum resolving set
                     (Steps 2-3, used for all graphs whose bases can be enumerated)
* ilp_metric_basis - one optimal metric basis by 0-1 integer linear programming
                     (hitting-set formulation, solved with SciPy/HiGHS)
"""
from itertools import combinations
import numpy as np
import networkx as nx


def distance_matrix(G):
    nodes = list(G.nodes())
    idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes)
    D = [[0] * n for _ in range(n)]
    for s, dist in nx.all_pairs_shortest_path_length(G):
        for t, d in dist.items():
            D[idx[s]][idx[t]] = d
    return D


def is_resolving(D, W):
    seen = set()
    for v in range(len(D)):
        key = tuple(D[v][u] for u in W)
        if key in seen:
            return False
        seen.add(key)
    return True


def all_metric_bases(D):
    """Every resolving set of minimum size (list of tuples of vertex indices)."""
    n = len(D)
    if n == 1:
        return [(0,)]
    for k in range(1, n):
        found = [W for W in combinations(range(n), k) if is_resolving(D, W)]
        if found:
            return found
    raise ValueError("no resolving set found")


def ilp_metric_basis(D, time_limit=None):
    """One optimal metric basis via min sum x_v s.t. sum_{v in R(x,y)} x_v >= 1."""
    from scipy.optimize import milp, LinearConstraint, Bounds
    n = len(D)
    rows = []
    for x, y in combinations(range(n), 2):
        rows.append([1.0 if D[v][x] != D[v][y] else 0.0 for v in range(n)])
    A = np.array(rows)
    opts = {} if time_limit is None else {"time_limit": time_limit}
    res = milp(c=np.ones(n), constraints=LinearConstraint(A, lb=1, ub=np.inf),
               integrality=np.ones(n), bounds=Bounds(0, 1), options=opts)
    if not res.success:
        raise RuntimeError(res.message)
    return tuple(int(v) for v in np.flatnonzero(res.x > 0.5))
