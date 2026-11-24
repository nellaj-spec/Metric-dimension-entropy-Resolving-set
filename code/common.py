"""Shared helpers: paths, graph (de)serialisation and the alkane generator."""
import json
from pathlib import Path
import networkx as nx

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
for p in (RESULTS, FIGURES):
    p.mkdir(exist_ok=True)


def all_alkane_trees(nmax=10):
    """All constitutional isomers of the alkanes C2..C_nmax as trees (max degree 4)."""
    out = []
    for n in range(2, nmax + 1):
        trees = [nx.path_graph(2)] if n == 2 else list(nx.nonisomorphic_trees(n))
        out += [T for T in trees if max(d for _, d in T.degree()) <= 4]
    return out


def tree_to_smiles(T):
    """SMILES of the carbon skeleton of an alkane tree."""
    start = min((v for v in T if T.degree(v) == 1), default=next(iter(T)))

    def rec(v, parent):
        kids = [w for w in T[v] if w != parent]
        s = "C"
        for w in kids[:-1]:
            s += "(" + rec(w, v) + ")"
        if kids:
            s += rec(kids[-1], v)
        return s
    return rec(start, None)


def edges_to_json(G):
    G = nx.convert_node_labels_to_integers(G)
    return json.dumps(sorted([min(u, v), max(u, v)] for u, v in G.edges()))


def graph_from_json(s, n=None):
    G = nx.Graph()
    if n:
        G.add_nodes_from(range(n))
    G.add_edges_from(json.loads(s))
    return G
