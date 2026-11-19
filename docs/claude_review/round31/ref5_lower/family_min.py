"""Exact minimum point degree of c=0 clusters INSIDE Astra's family (fixed V=F_2[t]_{<d},
Y=T+theta, sides of forms with l(1)=1, division by the full common gcd), small d.

For a selection X of sides: D0 = deg F (common), g = deg gcd, req = max_pairs 2deg(F_i-F_j)-D0,
valid iff g >= req; point degree = D0 - g.  DFS with pruning (req increases, g decreases).
deg(F_i-F_j) = deg P_J(1) + sum_{v in S_i cap S_j} deg(Y+v)  (checked exactly in check_formula.py).
Usage: family_min.py d maxdegT
"""
import sys, itertools
import numpy as np
from f4poly import *
from enum_gcd import prime_data


def run(d, T, best):
    Q = 1 << d
    fs = forms_and_sides(d)
    primes, ordm, degs, fdeg = prime_data(d, T)
    D0 = sum(fdeg) // 2
    sides = []
    for fi, (w, k, c) in enumerate(fs):
        for S in (k, c):
            sides.append((fi, S, tuple(int(x) for x in ordm[:, sorted(S)].sum(axis=1))))
    nS = len(sides)
    PJ1 = {}
    for a, b in itertools.combinations(range(len(fs)), 2):
        J = fs[a][1] & fs[b][1]
        PJ1[(a, b)] = sum(f2deg(1 ^ v) for v in J)
    req = [[None] * nS for _ in range(nS)]
    for x in range(nS):
        for y in range(nS):
            fx, Sx, _ = sides[x]
            fy, Sy, _ = sides[y]
            if fx < fy:
                req[x][y] = 2 * (PJ1[(fx, fy)] + sum(fdeg[v] for v in Sx & Sy)) - D0
    degl = list(int(x) for x in degs)
    BIG = 10**9

    def gdeg(mv):
        return sum(a * b for a, b in zip(mv, degl))

    # DFS over forms in order; each form: skip / side0 / side1
    def visit(fi, chosen, mv, r):
        if fi == len(fs):
            s = len(chosen)
            if s >= 2:
                pd = D0 - gdeg(mv)
                if s not in best or pd < best[s][0]:
                    best[s] = (pd, T, tuple(chosen), D0, r)
            return
        visit(fi + 1, chosen, mv, r)
        for side in (2 * fi, 2 * fi + 1):
            nmv = tuple(min(a, b) for a, b in zip(mv, sides[side][2])) if chosen else sides[side][2]
            nr = r
            for x in chosen:
                nr = max(nr, req[x][side])
            if len(chosen) >= 1 and nr > gdeg(nmv):
                continue  # prune: no extension can be valid
            visit(fi + 1, chosen + [side], nmv, nr)

    visit(0, [], (), -BIG)


if __name__ == "__main__":
    d = int(sys.argv[1])
    maxdeg = int(sys.argv[2])
    Q = 1 << d
    best = {}
    Ts = [0] + [Q * k for k in range(1, 1 << (maxdeg + 1 - d))]
    for T in Ts:
        run(d, T, best)
    print(f"d={d}, T over all classes mod V with deg T <= {maxdeg} ({len(Ts)} classes)")
    print("s  min_point_degree  attained_at_T(deg)  D0  req  ;  s^2*(d-3) ; Astra (V.1) D for m=s=2^(d-1): ")
    for s in sorted(best):
        pd, T, ch, D0, r = best[s]
        print(f"{s:2d}  {pd:6d}   T={T} (deg {f2deg(T) if T else '-inf'})  D0={D0} req={r}  sides={ch}")
    m = Q // 2
    print("V.1 point degree for m =", m, ":", (m - 1) * (d - 2) * 2 ** (d - 1))
