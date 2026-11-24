"""Alkane analyses: Tables 5-8 and the statistics of Sections 7.3-7.5."""
import pandas as pd
from scipy import stats
from common import DATA, RESULTS

ALK = "Alkanes C2-C10 (all isomers)"
MEASURES = [("H_distance", "Distance entropy"), ("H_vonneumann", "von Neumann entropy"),
            ("balaban", "Balaban index J"), ("H_bonchev", "Bonchev-Trinajstic entropy"),
            ("I", "Proposed I(G)"), ("randic", "Randic index R"), ("wiener", "Wiener index W"),
            ("I_max", "Proposed, I_max(G) only"), ("M2", "Zagreb index M2"),
            ("H_orbit", "Orbit entropy"), ("H_degree", "Degree entropy"),
            ("M1", "Zagreb index M1"), ("platt", "Platt index F")]
CORR = [("n", "Order n"), ("metric_dimension", "Metric dimension"), ("diameter", "Diameter"),
        ("radius", "Radius"), ("wiener", "Wiener index"), ("M1", "Zagreb index M1"),
        ("M2", "Zagreb index M2"), ("randic", "Randic index"), ("balaban", "Balaban index"),
        ("platt", "Platt index"), ("H_degree", "Degree entropy"), ("H_distance", "Distance entropy"),
        ("H_bonchev", "Bonchev-Trinajstic entropy"), ("H_orbit", "Orbit entropy"),
        ("H_vonneumann", "von Neumann entropy")]


def main():
    d = pd.read_csv(DATA / "descriptors_205.csv")
    a = d[d["class"] == ALK].copy()
    assert len(a) == 149

    # Table 5
    rows = []
    for n, g in a.groupby("n"):
        lin = g[g.metric_dimension == 1].iloc[0]
        br = g[g.metric_dimension > 1].I_norm
        rows.append({"n": n, "isomers": len(g),
                     "beta": f"{g.metric_dimension.min()}-{g.metric_dimension.max()}",
                     "I_n_alkane": round(lin.I, 3), "I_norm_n_alkane": round(lin.I_norm, 3),
                     "I_norm_branched_mean": round(br.mean(), 3) if len(br) else None,
                     "I_norm_branched_min": round(br.min(), 3) if len(br) else None,
                     "I_norm_branched_max": round(br.max(), 3) if len(br) else None})
    pd.DataFrame(rows).to_csv(RESULTS / "table5_alkanes_by_carbon_number.csv", index=False)

    # Table 6
    oc = a[a.n == 8].sort_values("I")
    oc[["name", "metric_dimension", "n_metric_bases", "I", "I_norm", "wiener"]].round(3).to_csv(
        RESULTS / "table6_octanes.csv", index=False)

    # Table 7: number of distinct values (8 decimals) per carbon number / number of isomers
    rows = []
    for col, lab in MEASURES:
        r = {"measure": lab}
        tot = 0
        for n, g in a.groupby("n"):
            k = g[col].round(8).nunique()
            tot += k
            if n in (8, 9, 10):
                r[f"C{n} ({len(g)})"] = round(k / len(g), 3)
        r["C2-C10 (149)"] = round(tot / 149, 3)
        r["distinct values"] = tot
        rows.append(r)
    pd.DataFrame(rows).to_csv(RESULTS / "table7_discrimination.csv", index=False)

    # Table 8
    rows = []
    subsets = [("All alkanes (149)", a), ("Octanes (18)", a[a.n == 8]),
               ("Nonanes (35)", a[a.n == 9]), ("Decanes (75)", a[a.n == 10])]
    for col, lab in CORR:
        r = {"descriptor": lab}
        for sl, s in subsets:
            if s[col].nunique() < 2:
                r[sl] = "-"
            else:
                r[sl] = f"{stats.pearsonr(s.I, s[col])[0]:.2f} / {stats.spearmanr(s.I, s[col])[0]:.2f}"
        rows.append(r)
    pd.DataFrame(rows).to_csv(RESULTS / "table8_correlations.csv", index=False)

    # statistics quoted in the text
    unb, br = a[a.metric_dimension == 1].I_norm, a[a.metric_dimension > 1].I_norm
    mw = stats.mannwhitneyu(unb, br)
    sp = stats.spearmanr(a.I_norm, a.metric_dimension)
    dec = a[a.n == 10]
    with open(RESULTS / "stats_alkanes.txt", "w") as f:
        f.write(f"Mann-Whitney unbranched (n={len(unb)}) vs branched (n={len(br)}): U={mw.statistic:.0f}, p={mw.pvalue:.2e}\n")
        f.write(f"Spearman I_norm vs metric dimension (149 alkanes): rho={sp.statistic:.3f}, p={sp.pvalue:.2e}\n")
        f.write("Mean I_norm by metric dimension: " + str(a.groupby('metric_dimension').I_norm.mean().round(3).to_dict()) + "\n")
        for col in ("H_orbit", "H_distance"):
            s = stats.spearmanr(dec.I, dec[col])
            f.write(f"Decanes, Spearman I vs {col}: rho={s.statistic:.3f}, p={s.pvalue:.3f}\n")
    print(open(RESULTS / "stats_alkanes.txt").read())
    print("wrote results/table5-8")


if __name__ == "__main__":
    main()
