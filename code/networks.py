"""Engineered and benchmark networks (Table 11).

All metric bases are enumerated for every network except the Zachary karate club,
for which C(34, 11) subsets would have to be examined; there one optimal metric
basis is found by integer linear programming and the single-basis entropy I_W(G)
is reported (it is NOT the basis average I(G) of Definition 5.1).
"""
import math
import networkx as nx
import pandas as pd
from common import DATA, RESULTS
from metric_basis import distance_matrix, all_metric_bases, ilp_metric_basis
from entropy import entropy_for_set, eccentricity_entropy_for_set


def ieee14():
    e = pd.read_csv(DATA / "ieee14_bus_edges.csv")
    return nx.Graph(list(zip(e.from_bus, e.to_bus)))


NETWORKS = [("IEEE 14-bus power system", ieee14, False),
            ("4 x 4 grid (sensor field)", lambda: nx.grid_2d_graph(4, 4), False),
            ("5 x 5 grid (sensor field)", lambda: nx.grid_2d_graph(5, 5), False),
            ("Hypercube Q4 (interconnect)", lambda: nx.hypercube_graph(4), False),
            ("Petersen graph", nx.petersen_graph, False),
            ("Florentine families (social)", nx.florentine_families_graph, False),
            ("Zachary karate club (social)", lambda: nx.Graph(nx.karate_club_graph().edges()), True)]


def main():
    rows = []
    for name, make, use_ilp in NETWORKS:
        G = nx.convert_node_labels_to_integers(make())
        D = distance_matrix(G)
        n, m = G.number_of_nodes(), G.number_of_edges()
        bases = [ilp_metric_basis(D)] if use_ilp else all_metric_bases(D)
        b = len(bases[0])
        I = sum(entropy_for_set(D, W) for W in bases) / len(bases)
        Iecc = sum(eccentricity_entropy_for_set(D, W) for W in bases) / len(bases)
        rows.append({"network": name, "n": n, "m": m, "D": nx.diameter(G), "beta": b,
                     "bases": "1 (ILP)" if use_ilp else len(bases),
                     "value": "I_W(G), single basis" if use_ilp else "I(G), basis average",
                     "I_bits": round(I, 3), "I_norm": round(I / math.log2(n), 3),
                     "I_ecc/log2n": round(Iecc / math.log2(n), 3),
                     "upper_bound": round(math.log2(n + b) - 2 * b / (n + b), 3),
                     "ilp_basis_vertices": " ".join(map(str, bases[0])) if use_ilp else ""})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "table11_networks.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
