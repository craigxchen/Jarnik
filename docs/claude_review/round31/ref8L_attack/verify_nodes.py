"""Exact integer verification of the node (weighted Prouhet-Tarry-Escott) families.

Construction (proved in the report, Prop. 2): nodes n_1..n_k (distinct integers), weights v_j != 0,
r >= 1.  For every 0/1 vector e with  sum_j e_j v_j n_j^p = 0  (p = 0..r-1)  put
        Y_e(T) = sum_j e_j v_j prod_{l != j} (T - n_l)          (a polynomial of degree <= k-1-r).
Modulo m_j = T - n_j one has Y_e(T) = e_j v_j prod_{l != j}(n_j - n_l): two classes.  Hence
q2 >= lcm_j |T - n_j| ~ T^k and diam <= C T^(k-1-r): ratio -> k/(k-1-r).
"""
import sys, itertools
from math import gcd, log, lcm
from fractions import Fraction as Fr

def patterns(nodes, v, r):
    k = len(nodes); out = []
    for e in itertools.product((0, 1), repeat=k):
        if all(sum(e[j] * v[j] * nodes[j] ** p for j in range(k)) == 0 for p in range(r)):
            out.append(e)
    return out

def points(nodes, v, pats, T):
    k = len(nodes)
    prods = []
    for j in range(k):
        P = 1
        for l in range(k):
            if l != j: P *= (T - nodes[l])
        prods.append(P)
    return [sum(e[j] * v[j] * prods[j] for j in range(k)) for e in pats]

def check(name, nodes, v, r, Ts):
    pats = patterns(nodes, v, r); k = len(nodes)
    for T in Ts:
        Y = points(nodes, v, pats, T)
        assert len(set(Y)) == len(Y), "points not distinct"
        mods = [abs(T - n) for n in nodes]
        two = all(len({y % m for y in Y}) <= 2 for m in mods)
        Q = 1
        for m in mods: Q = lcm(Q, m)
        L = max(Y) - min(Y)
        print(f"{name}: k={k} r={r} d={k-1-r} B={len(Y)} T=10^{len(str(T))-1}: two classes at all nodes={two}; "
              f"log lcm/log diam = {log(Q)/log(L):.4f} (limit {k/(k-1-r):.4f})")

def from_vectors(Z):
    """r=2: integer vectors (x,y), pairwise non-parallel -> nodes/weights after a shear making x != 0."""
    for c in range(0, 50):
        xs = [x + c * y for x, y in Z]
        if all(xs):
            break
    # node n = y/x, weight v = x ; clear denominators by scaling nodes: use T -> rational nodes
    nodes = [Fr(y, x) for (x, y), x in zip(Z, xs)]
    v = [Fr(x) for x in xs]
    D = 1
    for n in nodes: D = lcm(D, n.denominator)
    # scaling nodes by D (n -> D n) multiplies moment p by D^p: zero-moment patterns unchanged
    return [int(n * D) for n in nodes], [int(x) for x in v]

if __name__ == "__main__":
    Ts = [10**6 + 7, 10**30 + 7]
    check("r1k4", [0, 1, 2, 3], [1, 1, -1, -1], 1, Ts)
    check("r1k5", [0, 1, 2, 3, 4], [1, 1, 1, -1, -1], 1, Ts)
    check("r1k6", [0, 1, 2, 3, 4, 5], [1, 1, 1, -1, -1, -1], 1, Ts)
    for name, Z in [("r2k6", [(-4,-4),(-4,-3),(-4,-1),(2,0),(2,1),(2,3)]),
                    ("r2k7", [(-3,-3),(-3,-2),(0,-1),(1,2),(1,3),(2,0),(2,1)]),
                    ("r2k8", [(-3,-3),(-3,-2),(-3,-1),(0,-1),(1,2),(1,3),(2,-1),(2,0)])]:
        nodes, v = from_vectors(Z)
        check(name, nodes, v, 2, Ts)
    if len(sys.argv) > 1:
        pass
