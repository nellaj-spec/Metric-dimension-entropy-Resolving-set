# Metric-dimension entropy of graphs and molecules

Data and code for the article

> N. J. Jacob and Manjusha R, *A Graph-Theoretic Entropy Based on Metric Dimension and Its Application to Molecular Structures*, submitted to *Journal of Cheminformatics*.

The repository contains the complete molecular data set (205 molecules), all computed descriptors, and the Python code that reproduces every computed table, statistic and data figure of the article.

## The measure in brief

For a connected graph *G* and a metric basis *W*, each vertex *v* receives the membership
μ_W(v) = 1 / (1 + min_{u∈W} d(v, u)). After normalisation, p_W(v) = μ_W(v) / Σ_x μ_W(x), and the entropy is
I_W(G) = −Σ_v p_W(v) log₂ p_W(v). The basis-averaged entropy *I(G)* (Definition 5.1) is the mean of I_W(G) over **all** metric bases of *G*; I_min and I_max are the minimum and maximum, and I_norm = I(G) / log₂ n.

## Repository layout

```
data/
  molecules_205.csv        (1) the 205 molecules: id, name, class, SMILES, n, m, edge list   [generated]
  descriptors_205.csv      (2) metric dimension, number of metric bases, I, I_min, I_max, I_norm,
                               five comparison entropies and six topological indices           [generated]
  ring_molecules.csv       input: names, classes and SMILES of the 56 ring-containing molecules
  boiling_points_67.csv    (6) input: experimental normal boiling points of 67 alkanes, ref. [31]
  ieee14_bus_edges.csv     input: branch list of the IEEE 14-bus power test system
code/
  metric_basis.py          (3) distance matrix, enumeration of all metric bases, ILP metric basis
  entropy.py               (4) the proposed entropy (Algorithm 1) and the eccentricity-based variant
  descriptors.py           (5) comparison graph entropies and topological indices
  boiling_point.py         (6) leave-one-out boiling-point prediction                 -> Table 9, Fig. 15
  sensitivity.py           (6) single-perturbation sensitivity analysis                -> Table 10, Fig. 16
  build_dataset.py         builds data/molecules_205.csv and data/descriptors_205.csv
  analysis_dataset.py      data-set summary, basis dependence, ring molecules         -> Table 3, Sec. 7.3, 7.6
  analysis_alkanes.py      alkane analyses                                            -> Tables 5-8, Sec. 7.3-7.5
  networks.py              engineered and benchmark networks                          -> Table 11
  verify_theory.py         numerical check of Proposition 6.2 and Theorem 6.1         -> Section 6
  timing.py                run time of the exact metric-basis computation             -> Section 5.3
  make_figures.py          data figures                                               -> Figures 10-16
  alkane_names.py          IUPAC-style names for the alkane isomers, and a name parser
  smiles_graph.py          minimal SMILES reader (hydrogen-suppressed graph)
  common.py                paths and shared helpers
  run_all.py               runs everything in the right order
results/                   all tables (CSV) and statistics (TXT)                      [generated]
figures/                   Figure10.png ... Figure16.png                              [generated]
```

The numbers (1)–(7) refer to the items listed in the article's *Availability of data and materials* section; item (7) is this README.

## Requirements and installation

Python 3.10 or later and the packages in `requirements.txt` (NetworkX, NumPy, SciPy, pandas, scikit-learn, matplotlib). RDKit is **not** required.

```bash
git clone [REPOSITORY URL]
cd metric-dimension-entropy
python -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Reproducing the results

```bash
python code/run_all.py              # everything except the timing experiment, about 2 minutes
python code/run_all.py --timing     # additionally the run-time experiment of Section 5.3
```

Each script can also be run on its own from inside `code/` (run `build_dataset.py` first, because the other scripts read `data/descriptors_205.csv`).

| Article item | Script | Output |
|---|---|---|
| Data set, Table 3 | `build_dataset.py`, `analysis_dataset.py` | `data/*.csv`, `results/table3_dataset.csv` |
| Basis dependence (Sec. 7.3), ring molecules (Sec. 7.6) | `analysis_dataset.py` | `results/stats_dataset.txt` |
| Tables 5, 6, 7, 8; Mann–Whitney and Spearman tests | `analysis_alkanes.py` | `results/table5_*.csv` … `table8_*.csv`, `results/stats_alkanes.txt` |
| Table 9, Wilcoxon test, within-isomer correlations | `boiling_point.py` | `results/table9_boiling_points.csv`, `results/stats_boiling_points.txt`, `results/boiling_point_predictions.csv` |
| Table 10 | `sensitivity.py` | `results/table10_sensitivity.csv`, `results/sensitivity_all_perturbations.csv` |
| Table 11 | `networks.py` | `results/table11_networks.csv` |
| Proposition 6.2, Theorem 6.1 | `verify_theory.py` | `results/verify_theory.txt` |
| Section 5.3 run times | `timing.py` | `results/timing.txt` |
| Figures 10–16 | `make_figures.py` | `figures/Figure10.png` … `Figure16.png` |

Figures 1–9 are explanatory drawings and are not generated by code. Tables 1 and 2 are hand calculations for complete graphs and cycles; `verify_theory.py` checks the corresponding closed forms.

## Notes on the data and methods

* **Alkanes.** All 149 constitutional isomers of C₂–C₁₀ are generated as the non-isomorphic trees of maximum degree 4 (`networkx.nonisomorphic_trees`). Their names are produced by `alkane_names.py` and checked by parsing each name back into a structure.
* **Ring-containing molecules.** The 56 molecules are read from SMILES in `data/ring_molecules.csv`. Only connectivity is used, so bond orders and atom types are ignored; for this reason some molecules share a graph (for example benzene and cyclohexane), and the 56 molecules give 50 distinct graphs.
* **Edge lists.** `edge_list` in `molecules_205.csv` is a JSON list of vertex pairs with vertices numbered 0 … n−1.
* **Metric bases.** For every molecule all metric bases are enumerated exhaustively. For the Zachary karate club (Table 11), enumeration is infeasible, so one optimal basis is found by integer linear programming and the value reported is the single-basis entropy I_W(G), not the basis average.
* **Discrimination ratio (Table 7).** For each carbon number: the number of distinct values (rounded to 8 decimals) divided by the number of isomers. The overall ratio is the sum of the distinct values over all carbon numbers divided by 149.
* **Boiling points.** The 67 values in `data/boiling_points_67.csv` are copied from the supplementary material of ref. [31], with the names as printed there. The value printed for 2,2-dimethylhexane (10.84 °C) is a typographical error and 106.84 °C is used. Of the 72 alkane isomers C₄–C₉, the source gives no value for 2,2,3-trimethylbutane, 2,2,3,3-tetramethylbutane, 2,2,3-trimethylhexane, 2,4,4-trimethylhexane and 3-ethyl-4-methylhexane.
* **Random elements.** Random forests use `random_state = 0`. The random graphs of Table 10 are `nx.random_labeled_tree(12, seed=3)`, `nx.erdos_renyi_graph(12, 0.3, seed=11)` and `nx.barabasi_albert_graph(12, 2, seed=5)`. Timings in `timing.py` depend on the machine.

## Licence

The code is released under the MIT licence (see `LICENSE`). The data in `data/` may be reused with attribution to the article.

## How to cite

If you use this code or data, please cite the article above (full reference and DOI to be added on publication) and this repository ([Zenodo DOI to be inserted]).
