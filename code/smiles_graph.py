"""Minimal SMILES reader returning the hydrogen-suppressed molecular graph.

Only connectivity is kept (bond orders, charges, stereo and aromaticity flags are
ignored), because the entropy is defined on simple graphs.
"""
import re
import networkx as nx

_TOKEN = re.compile(r"\[[^\]]+\]|Br|Cl|[BCNOPSFI]|[bcnops]|\(|\)|%\d\d|\d|[=#:\-/\\.$]")


def smiles_to_graph(smiles):
    G = nx.Graph()
    stack, rings = [], {}
    prev, idx = None, 0
    for tok in _TOKEN.findall(smiles):
        if tok == "(":
            stack.append(prev)
        elif tok == ")":
            prev = stack.pop()
        elif tok[0] == "%" or tok.isdigit():
            r = tok
            if r in rings:
                G.add_edge(prev, rings.pop(r))
            else:
                rings[r] = prev
        elif tok in "=#:-/\\$":
            continue
        elif tok == ".":
            prev = None
        else:
            sym = tok.strip("[]")
            sym = re.match(r"\d*([A-Za-z][a-z]?)", sym).group(1)
            G.add_node(idx, element=sym)
            if prev is not None:
                G.add_edge(prev, idx)
            prev = idx
            idx += 1
    if rings:
        raise ValueError(f"unclosed ring in {smiles}")
    if "".join(_TOKEN.findall(smiles)) != smiles:
        raise ValueError(f"could not tokenise {smiles}")
    return G
