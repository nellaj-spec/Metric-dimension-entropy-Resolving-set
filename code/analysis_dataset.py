"""Whole data set: Table 3, basis dependence (Section 7.3) and ring molecules (Section 7.6)."""
import pandas as pd
from scipy import stats
from common import DATA, RESULTS

ORDER = ["Alkanes C2-C10 (all isomers)", "Cycloalkanes and cages", "Benzene derivatives",
         "Polycyclic aromatic hydrocarbons", "Drug-like molecules"]


def row(label, g):
    sp = g.I_max - g.I_min
    return {"class": label, "molecules": len(g), "n": f"{g.n.min()}-{g.n.max()}",
            "beta": f"{g.metric_dimension.min()}-{g.metric_dimension.max()}",
            "I_range_bits": f"{g.I.min():.3f}-{g.I.max():.3f}",
            "I_norm_mean": round(g.I_norm.mean(), 3), "I_norm_sd": round(g.I_norm.std(ddof=0), 3),
            "mean_Imax_minus_Imin": round(sp.mean(), 3)}


def main():
    d = pd.read_csv(DATA / "descriptors_205.csv")
    rows = [row(c, d[d["class"] == c]) for c in ORDER] + [row("All", d)]
    pd.DataFrame(rows).to_csv(RESULTS / "table3_dataset.csv", index=False)

    sp = d.I_max - d.I_min
    k = sp.idxmax()
    rings = d[d["class"] != ORDER[0]]
    m = pd.read_csv(DATA / "molecules_205.csv")
    from common import graph_from_json
    import networkx as nx
    reps = []
    for _, r in m[m["class"] != ORDER[0]].iterrows():
        G = graph_from_json(r.edge_list, r.n_vertices)
        if not any(nx.is_isomorphic(G, H) for H in reps):
            reps.append(G)
    with open(RESULTS / "stats_dataset.txt", "w") as f:
        f.write(f"I_max - I_min over 205 molecules: mean {sp.mean():.3f} bits "
                f"({100 * (sp / d.I).mean():.1f}% of I), max {sp.max():.3f} bits "
                f"({100 * sp[k] / d.I[k]:.1f}%, {d.name[k]})\n")
        f.write(f"Ring-containing molecules: {len(rings)}; Spearman rho(I, n) = {stats.spearmanr(rings.I, rings.n)[0]:.3f}\n")
        f.write(f"Distinct (non-isomorphic) graphs among the 56 ring-containing molecules: {len(reps)}\n")
        for name in ["cubane", "adamantane", "coronene", "mesitylene", "aspirin", "tetracene",
                     "anthracene", "phenanthrene", "chrysene", "triphenylene"]:
            r = d[d.name == name].iloc[0]
            f.write(f"  {name:13s} I = {r.I:.3f} bits, I_norm = {r.I_norm:.3f}\n")
    print(open(RESULTS / "stats_dataset.txt").read())


if __name__ == "__main__":
    main()
