"""Build the 205-molecule data set and compute all descriptors.

Outputs
  data/molecules_205.csv    id, name, class, SMILES, n, m and the edge list of the
                            hydrogen-suppressed graph (JSON, 0-based vertex indices)
  data/descriptors_205.csv  metric dimension, number of metric bases, I, I_min, I_max,
                            I_norm, the five comparison entropies and six indices
"""
import csv
import networkx as nx
import pandas as pd
from common import DATA, all_alkane_trees, tree_to_smiles, edges_to_json
from alkane_names import name_alkane, parse_alkane
from smiles_graph import smiles_to_graph
from entropy import metric_dimension_entropy
from descriptors import descriptors


def main():
    mols = []
    for T in all_alkane_trees(10):
        name = name_alkane(T)
        assert nx.is_isomorphic(parse_alkane(name), T), name
        smi = tree_to_smiles(T)
        assert nx.is_isomorphic(smiles_to_graph(smi), T), smi
        mols.append({"name": name, "class": "Alkanes C2-C10 (all isomers)", "smiles": smi, "G": T})
    for r in csv.DictReader(open(DATA / "ring_molecules.csv")):
        mols.append({"name": r["name"], "class": r["class"], "smiles": r["smiles"],
                     "G": smiles_to_graph(r["smiles"])})
    assert len(mols) == 205
    assert len({m["name"] for m in mols}) == 205

    mol_rows, desc_rows = [], []
    for i, m in enumerate(mols, start=1):
        G = nx.convert_node_labels_to_integers(m["G"])
        assert nx.is_connected(G)
        mid = f"M{i:03d}"
        mol_rows.append({"id": mid, "name": m["name"], "class": m["class"], "smiles": m["smiles"],
                         "n_vertices": G.number_of_nodes(), "n_edges": G.number_of_edges(),
                         "edge_list": edges_to_json(G)})
        e = metric_dimension_entropy(G)
        d = {"id": mid, "name": m["name"], "class": m["class"], "n": e["n"], "m": e["m"],
             "metric_dimension": e["beta"], "n_metric_bases": e["n_bases"],
             "I": e["I"], "I_min": e["I_min"], "I_max": e["I_max"], "I_norm": e["I_norm"]}
        d.update(descriptors(G))
        desc_rows.append(d)
        print(f"{mid} {m['name']:40s} n={e['n']:2d} beta={e['beta']} bases={e['n_bases']:4d} I={e['I']:.4f}")
    pd.DataFrame(mol_rows).to_csv(DATA / "molecules_205.csv", index=False)
    pd.DataFrame(desc_rows).to_csv(DATA / "descriptors_205.csv", index=False, float_format="%.10g")
    print("wrote data/molecules_205.csv and data/descriptors_205.csv")


if __name__ == "__main__":
    main()
