"""Figures 10-16 of the manuscript (data figures). Figures 1-9 are explanatory drawings."""
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from common import DATA, RESULTS, FIGURES

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})
ALK = "Alkanes C2-C10 (all isomers)"
LAB = {"I": "Proposed I(G)", "n": "Order n", "metric_dimension": "Metric dimension", "diameter": "Diameter",
       "wiener": "Wiener", "M1": "Zagreb M1", "M2": "Zagreb M2", "randic": "Randic", "balaban": "Balaban J",
       "platt": "Platt", "H_degree": "Degree entropy", "H_distance": "Distance entropy",
       "H_bonchev": "Bonchev-Trinajstic", "H_orbit": "Orbit entropy", "H_vonneumann": "von Neumann"}


def main():
    d = pd.read_csv(DATA / "descriptors_205.csv")
    a = d[d["class"] == ALK]
    rng = np.random.default_rng(0)

    # Figure 10
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))
    br = a.metric_dimension > 1
    ax[0].scatter(a.n[br] + rng.uniform(-.15, .15, br.sum()), a.I[br], s=10, alpha=.7, label="branched")
    ax[0].scatter(a.n[~br], a.I[~br], s=28, marker="s", color="C3", label="unbranched")
    ns = np.arange(2, 11)
    ax[0].plot(ns, np.log2(ns), "k--", lw=1, label="log2 n")
    ax[0].set(xlabel="Number of carbon atoms", ylabel="I(G) (bits)", title="(a)")
    ax[0].legend(frameon=False)
    groups = [a.I_norm[a.metric_dimension == b] for b in (1, 2, 3, 4)]
    ax[1].boxplot(groups, tick_labels=["1", "2", "3", "4"])
    ax[1].set(xlabel="Metric dimension", ylabel="Normalised entropy", title="(b)")
    fig.savefig(FIGURES / "Figure10.png"); plt.close(fig)

    # Figure 11
    t7 = pd.read_csv(RESULTS / "table7_discrimination.csv").sort_values("C2-C10 (149)")
    fig, ax = plt.subplots(figsize=(5, 4))
    cols = ["C3" if m.startswith("Proposed I(G)") else "C0" for m in t7.measure]
    ax.barh(t7.measure, t7["C2-C10 (149)"], color=cols)
    ax.set(xlabel="Discrimination ratio (149 alkanes)", xlim=(0, 1.05))
    fig.savefig(FIGURES / "Figure11.png"); plt.close(fig)

    # Figure 12
    keys = list(LAB)
    C = a[keys].corr(method="spearman").values
    fig, ax = plt.subplots(figsize=(6.5, 5.5))
    im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(keys)), [LAB[k] for k in keys], rotation=90)
    ax.set_yticks(range(len(keys)), [LAB[k] for k in keys])
    for i in range(len(keys)):
        for j in range(len(keys)):
            ax.text(j, i, f"{C[i, j]:.2f}", ha="center", va="center", fontsize=5.5)
    fig.colorbar(im, ax=ax, shrink=.8, label="Spearman rho")
    fig.savefig(FIGURES / "Figure12.png"); plt.close(fig)

    # Figure 13
    dec = a[a.n == 10]
    fig, ax = plt.subplots(1, 3, figsize=(9, 2.8), sharey=True)
    for k, (col, lab) in enumerate([("wiener", "Wiener index"), ("diameter", "Diameter"),
                                    ("metric_dimension", "Metric dimension")]):
        x = dec[col] + (rng.uniform(-.12, .12, len(dec)) if col != "wiener" else 0)
        ax[k].scatter(x, dec.I, s=10)
        rho = stats.spearmanr(dec.I, dec[col])[0]
        ax[k].set(xlabel=lab, title=f"({'abc'[k]}) rho = {rho:.2f}")
    ax[0].set_ylabel("I(G) (bits)")
    fig.savefig(FIGURES / "Figure13.png"); plt.close(fig)

    # Figure 14
    order = [ALK, "Cycloalkanes and cages", "Benzene derivatives", "Polycyclic aromatic hydrocarbons", "Drug-like molecules"]
    fig, ax = plt.subplots(figsize=(6, 3.4))
    ax.boxplot([d.I_norm[d["class"] == c] for c in order],
               tick_labels=["Alkanes", "Cycloalkanes\n& cages", "Benzene\nderivatives", "PAHs", "Drug-like"])
    ax.set_ylabel("Normalised entropy")
    fig.savefig(FIGURES / "Figure14.png"); plt.close(fig)

    # Figure 15
    p = pd.read_csv(RESULTS / "boiling_point_predictions.csv")
    sets = [c for c in p.columns if c.startswith("loo_pred")]
    fig, ax = plt.subplots(1, 3, figsize=(9, 3), sharey=True)
    for k, c in enumerate(sets):
        rmse = math.sqrt(((p.bp - p[c]) ** 2).mean())
        ax[k].scatter(p.bp, p[c], s=10)
        lim = [p.bp.min() - 10, p.bp.max() + 10]
        ax[k].plot(lim, lim, "k--", lw=1)
        short = ["n only", "n + proposed descriptors", "n + classical indices"][k]
        ax[k].set(xlabel="Observed b.p. (C)", title=f"({'abc'[k]}) {short}\nRMSE {rmse:.2f} C")
    ax[0].set_ylabel("LOO predicted b.p. (C)")
    fig.savefig(FIGURES / "Figure15.png"); plt.close(fig)

    # Figure 16
    s = pd.read_csv(RESULTS / "sensitivity_all_perturbations.csv")
    s["absd"] = s.delta.abs()
    piv = s.pivot_table(index="graph", columns="operation", values="absd", aggfunc="mean")
    piv = piv.reindex(s.graph.drop_duplicates())
    fig, ax = plt.subplots(figsize=(8, 3.4))
    piv.plot.bar(ax=ax, width=.8)
    ax.set(ylabel="Mean |change of entropy|", xlabel="")
    ax.set_xticklabels(piv.index, rotation=25, ha="right")
    ax.legend(frameon=False)
    fig.savefig(FIGURES / "Figure16.png"); plt.close(fig)
    print("wrote figures/Figure10-16.png")


if __name__ == "__main__":
    main()
