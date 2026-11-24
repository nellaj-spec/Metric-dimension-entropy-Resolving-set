"""Run-time of Step 3 (exact metric basis by integer programming) on random graphs (Section 5.3).

Connected Erdos-Renyi graphs with expected degree 2.5 ln n. Timings depend on the
machine and will differ from the values in the manuscript; the trend is what matters.
Not run by run_all.py (pass --timing to include it).
"""
import math
import time
import networkx as nx
from common import RESULTS
from metric_basis import distance_matrix, ilp_metric_basis
from entropy import entropy_for_set


def main(sizes=(10, 20, 30, 40, 50), seed=1, time_limit=120):
    lines = []
    for n in sizes:
        p = 2.5 * math.log(n) / (n - 1)
        s = seed
        while True:
            G = nx.erdos_renyi_graph(n, p, seed=s)
            if nx.is_connected(G):
                break
            s += 1
        D = distance_matrix(G)
        t0 = time.perf_counter()
        W = ilp_metric_basis(D, time_limit=time_limit)
        t1 = time.perf_counter()
        entropy_for_set(D, W)
        t2 = time.perf_counter()
        lines.append(f"n={n:3d}  beta={len(W):2d}  ILP {t1 - t0:8.2f} s   entropy {1000 * (t2 - t1):.3f} ms")
        print(lines[-1])
    open(RESULTS / "timing.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
