"""Numerical checks of Section 6: closed forms (Proposition 6.2) and bounds (Theorem 6.1)."""
import math
import networkx as nx
import pandas as pd
from common import DATA, RESULTS, graph_from_json
from metric_basis import distance_matrix, all_metric_bases
from entropy import entropy_for_set, metric_dimension_entropy


def closed_forms(n):
    H = sum(1 / i for i in range(1, n + 1))
    path = math.log2(H) + sum(math.log2(i) / i for i in range(1, n + 1)) / H
    S = n - 7 / 6
    star = math.log2(S) + (0.5 + math.log2(3) / 3) / S
    comp = math.log2(n - 0.5) + 1 / (2 * n - 1)
    return path, star, comp


def bounds_hold(G):
    G = nx.convert_node_labels_to_integers(G)
    D = distance_matrix(G)
    n, diam = len(D), max(max(r) for r in D)
    ok = True
    for W in all_metric_bases(D):
        b = len(W)
        I = entropy_for_set(D, W)
        lo = math.log2(b + (n - b) / (1 + diam))
        hi = math.log2(n + b) - 2 * b / (n + b)
        ok &= lo - 1e-9 <= I <= hi + 1e-9 and hi < math.log2(n)
    return ok


def main():
    lines = []
    for n in range(3, 13):
        p, s, c = closed_forms(n)
        P = metric_dimension_entropy(nx.path_graph(n))["I"]
        K = metric_dimension_entropy(nx.complete_graph(n))["I"]
        line = f"n={n:2d}  path {P:.6f}/{p:.6f}  complete {K:.6f}/{c:.6f}"
        assert abs(P - p) < 1e-9 and abs(K - c) < 1e-9
        if n >= 4:
            St = metric_dimension_entropy(nx.star_graph(n - 1))["I"]
            assert abs(St - s) < 1e-9
            line += f"  star {St:.6f}/{s:.6f}"
        lines.append(line)
    mol = pd.read_csv(DATA / "molecules_205.csv")
    bad = [r.name for r in mol.itertuples() if r.n_vertices > 1 and not bounds_hold(graph_from_json(r.edge_list, r.n_vertices))]
    lines.append(f"Theorem 6.1 bounds hold for every metric basis of all 205 molecular graphs: {not bad}")
    txt = "Proposition 6.2 (computed / closed form)\n" + "\n".join(lines) + "\n"
    open(RESULTS / "verify_theory.txt", "w").write(txt)
    print(txt)


if __name__ == "__main__":
    main()
