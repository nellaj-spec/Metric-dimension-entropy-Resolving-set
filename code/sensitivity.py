"""Sensitivity of the entropy to single structural perturbations (Table 10, Figure 16).

For each graph every edge addition, every edge deletion and every vertex deletion
that leaves the graph connected is applied, and the basis-averaged entropy is
recomputed with all metric bases. For vertex deletion the change of I_norm is
reported, otherwise the change of I.
"""
import itertools
import networkx as nx
from entropy import metric_dimension_entropy


def perturbations(G):
    G = nx.convert_node_labels_to_integers(G)
    for u, v in itertools.combinations(G.nodes(), 2):
        if not G.has_edge(u, v):
            H = G.copy(); H.add_edge(u, v)
            yield "edge addition", H
    for u, v in list(G.edges()):
        H = G.copy(); H.remove_edge(u, v)
        if nx.is_connected(H):
            yield "edge deletion", H
    for v in list(G.nodes()):
        H = G.copy(); H.remove_node(v)
        if H.number_of_nodes() > 1 and nx.is_connected(H):
            yield "vertex deletion", nx.convert_node_labels_to_integers(H)


def sensitivity(G):
    base = metric_dimension_entropy(G)
    rows = []
    for op, H in perturbations(G):
        e = metric_dimension_entropy(H)
        key = "I_norm" if op == "vertex deletion" else "I"
        rows.append({"operation": op, "delta": e[key] - base[key],
                     "beta_changed": e["beta"] != base["beta"]})
    return base, rows


# ---------------------------------------------------------------------------
# Table 10 / Figure 16
GRAPHS = [("Path P12", lambda: nx.path_graph(12)),
          ("Cycle C12", lambda: nx.cycle_graph(12)),
          ("Grid 3 x 4", lambda: nx.grid_2d_graph(3, 4)),
          ("Random tree, n=12", lambda: nx.random_labeled_tree(12, seed=3)),
          ("Erdos-Renyi, n=12, p=0.3", lambda: nx.erdos_renyi_graph(12, 0.3, seed=11)),
          ("Barabasi-Albert, n=12, m=2", lambda: nx.barabasi_albert_graph(12, 2, seed=5)),
          ("2,2,4-Trimethylpentane", None)]


def main():
    import pandas as pd
    from common import RESULTS
    from smiles_graph import smiles_to_graph
    out, raw = [], []
    for name, make in GRAPHS:
        G = smiles_to_graph("CC(C)CC(C)(C)C") if make is None else make()
        base, rows = sensitivity(G)
        df = pd.DataFrame(rows)
        for r in rows:
            raw.append({"graph": name, **r})
        for op, g in df.groupby("operation", sort=False):
            out.append({"graph": name, "beta0": base["beta"], "I0": round(base["I"], 3), "operation": op,
                        "cases": len(g), "mean_delta": round(g.delta.mean(), 3),
                        "min_delta": round(g.delta.min(), 3), "max_delta": round(g.delta.max(), 3),
                        "pct_delta_positive": round(100 * (g.delta > 0).mean()),
                        "pct_beta_changed": round(100 * g.beta_changed.mean())})
    pd.DataFrame(out).to_csv(RESULTS / "table10_sensitivity.csv", index=False)
    pd.DataFrame(raw).to_csv(RESULTS / "sensitivity_all_perturbations.csv", index=False)
    print(pd.DataFrame(out).to_string(index=False))
    print("total perturbed graphs:", len(raw))


if __name__ == "__main__":
    main()
