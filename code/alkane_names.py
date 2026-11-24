"""IUPAC-style names for acyclic alkanes (C2-C10) and a parser for such names.

name_alkane(T): name of the tree T (maximum degree 4).
parse_alkane(name): tree for a substitutive alkane name (used to match the
boiling-point table of reference [31] to the generated structures).
Round-trip consistency is checked in build_dataset.py.
"""
import re
import networkx as nx

PARENT = {1: "methane", 2: "ethane", 3: "propane", 4: "butane", 5: "pentane", 6: "hexane",
          7: "heptane", 8: "octane", 9: "nonane", 10: "decane"}
MULT = {2: "di", 3: "tri", 4: "tetra"}


def _longest_paths(T):
    ends = [v for v in T if T.degree(v) <= 1]
    best, paths = -1, []
    for i, a in enumerate(ends):
        for b in ends[i + 1:]:
            p = nx.shortest_path(T, a, b)
            if len(p) > best:
                best, paths = len(p), [p]
            elif len(p) == best:
                paths.append(p)
    if not paths:                      # single vertex
        paths = [[next(iter(T))]]
    return paths


def _branch_name(T, root, parent):
    """Name of the substituent rooted at `root` (attached to `parent`)."""
    sub = nx.bfs_tree(T.subgraph([v for v in T if v != parent]), root).to_undirected()
    k = sub.number_of_nodes()
    if k == 1:
        return "methyl"
    if k == 2:
        return "ethyl"
    if k == 3:
        return "propyl" if sub.degree(root) == 1 else "isopropyl"
    raise ValueError("substituent with more than 3 carbons")


def _candidates(T):
    for chain in _longest_paths(T):
        for c in (chain, chain[::-1]):
            subs = []
            cs = set(c)
            try:
                for i, v in enumerate(c, start=1):
                    for w in T[v]:
                        if w not in cs:
                            subs.append((i, _branch_name(T, w, v)))
            except ValueError:         # chain leaves a complex substituent: never preferred
                continue
            yield c, sorted(subs)


def name_alkane(T):
    if T.number_of_nodes() == 1:
        return "methane"
    best = None
    for c, subs in _candidates(T):
        key = (-len(c), -len(subs), [i for i, _ in subs],
               [s for _, s in sorted(subs, key=lambda x: (x[1], x[0]))])
        if best is None or key < best[0]:
            best = (key, c, subs)
    _, c, subs = best
    groups = {}
    for i, s in subs:
        groups.setdefault(s, []).append(i)
    parts = []
    for s in sorted(groups):           # alphabetical (multiplying prefixes ignored)
        loc = groups[s]
        parts.append(",".join(map(str, loc)) + "-" + MULT.get(len(loc), "") + s)
    name = "-".join(parts) + PARENT[len(c)]
    return name if not parts else name


PARENT_INV = {v: k for k, v in PARENT.items()}
SUB = {"methyl": [(0, None)], "ethyl": [(0, None), (1, 0)], "propyl": [(0, None), (1, 0), (2, 1)],
       "isopropyl": [(0, None), (1, 0), (2, 0)]}


def parse_alkane(name):
    s = name.lower().replace(" ", "").replace("methyal", "methyl")
    m = re.search(r"(methane|ethane|propane|butane|pentane|hexane|heptane|octane|nonane|decane)s?$", s)
    L = PARENT_INV[m.group(1)]
    T = nx.path_graph(L)
    nxt = L
    for loc, mult, sub in re.findall(r"([\d,]+)-(di|tri|tetra)?(methyl|ethyl|propyl|isopropyl)", s[:m.start()]):
        locs = [int(x) for x in loc.split(",")]
        for x in locs:
            ids = []
            for j, (k, par) in enumerate(SUB[sub]):
                ids.append(nxt); T.add_node(nxt)
                T.add_edge(nxt, x - 1 if par is None else ids[par]); nxt += 1
    return T
