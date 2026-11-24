"""Structure-property study: leave-one-out prediction of alkane boiling points (Table 9).

The 67 experimental normal boiling points (data/boiling_points_67.csv) are taken
from the supplementary material of reference [31]; each name is matched to its
structure by parsing it and testing graph isomorphism.
"""
import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from common import DATA, RESULTS, graph_from_json
from alkane_names import parse_alkane

FIX = {"2-methyl propane": "2-methylpropane", "2-methyl butane": "2-methylbutane",
       "2,2 dimethylpropane": "2,2-dimethylpropane", "2,2-methylbutane": "2,2-dimethylbutane",
       "heptanes": "heptane", "3-methyalpentane": "3-methylpentane"}
SETS = {"n": ["n"], "n, W": ["n", "wiener"], "n, I": ["n", "I"],
        "n, I, I_min, I_max, beta (proposed)": ["n", "I", "I_min", "I_max", "metric_dimension"],
        "n, W, M1, M2, R, J (classical indices)": ["n", "wiener", "M1", "M2", "randic", "balaban"],
        "n and the five entropies of Table 4": ["n", "H_degree", "H_distance", "H_bonchev", "H_orbit", "H_vonneumann"],
        "Classical indices and proposed descriptors": ["n", "wiener", "M1", "M2", "randic", "balaban",
                                                       "I", "I_min", "I_max", "metric_dimension"]}
SEED = 0


def load():
    bp = pd.read_csv(DATA / "boiling_points_67.csv")
    mol = pd.read_csv(DATA / "molecules_205.csv")
    desc = pd.read_csv(DATA / "descriptors_205.csv")
    alk = mol[mol["class"].str.startswith("Alkanes")]
    graphs = {r.id: graph_from_json(r.edge_list, r.n_vertices) for r in alk.itertuples()}
    ids = []
    for nm in bp.name_in_source:
        key = nm.lower()
        G = parse_alkane(FIX.get(key, key))
        hit = [i for i, H in graphs.items() if H.number_of_nodes() == G.number_of_nodes() and nx.is_isomorphic(G, H)]
        assert len(hit) == 1, nm
        ids.append(hit[0])
    assert len(set(ids)) == 67
    d = desc.set_index("id").loc[ids].reset_index()
    d["bp"] = bp.bp_celsius.values
    return d


def scores(y, p):
    return 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum(), np.sqrt(((y - p) ** 2).mean())


def main():
    d = load()
    y = d.bp.values
    rows, preds = [], {}
    for lab, cols in SETS.items():
        X = d[cols].values
        pl = cross_val_predict(LinearRegression(), X, y, cv=LeaveOneOut())
        pr = cross_val_predict(RandomForestRegressor(n_estimators=300, random_state=SEED), X, y, cv=LeaveOneOut())
        preds[lab] = pl
        (r2l, el), (r2r, er) = scores(y, pl), scores(y, pr)
        rows.append({"feature set": lab, "linear R2": round(r2l, 3), "linear RMSE (C)": round(el, 2),
                     "random forest R2": round(r2r, 3), "random forest RMSE (C)": round(er, 2)})
    pd.DataFrame(rows).to_csv(RESULTS / "table9_boiling_points.csv", index=False)
    out = d[["id", "name", "n", "bp"]].copy()
    for k in ("n", "n, I, I_min, I_max, beta (proposed)", "n, W, M1, M2, R, J (classical indices)"):
        out["loo_pred: " + k] = preds[k]
    out.to_csv(RESULTS / "boiling_point_predictions.csv", index=False)

    a = np.abs(y - preds["n, W, M1, M2, R, J (classical indices)"])
    b = np.abs(y - preds["Classical indices and proposed descriptors"])
    with open(RESULTS / "stats_boiling_points.txt", "w") as f:
        f.write(f"Alkanes per carbon number: {d.groupby('n').size().to_dict()}\n")
        f.write(f"Pearson r(I, bp), all 67: {stats.pearsonr(d.I, d.bp)[0]:.2f}\n")
        for n in (8, 9):
            s = d[d.n == n]
            f.write(f"C{n} ({len(s)} alkanes): r(I, bp) = {stats.pearsonr(s.I, s.bp)[0]:.2f}, "
                    f"r(Randic, bp) = {stats.pearsonr(s.randic, s.bp)[0]:.2f}\n")
        f.write(f"Wilcoxon signed-rank, |errors| classical vs classical+proposed: p = {stats.wilcoxon(a, b).pvalue:.2f}\n")
        f.write(f"Random forests: 300 trees, random_state = {SEED}\n")
    print(open(RESULTS / "stats_boiling_points.txt").read())
    print(pd.DataFrame(rows).to_string(index=False))


if __name__ == "__main__":
    main()
