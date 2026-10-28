"""K-degrees of F-curves on Mbar_{0,8}, exact (Fraction).
K = sum_S (s(n-s)/(n-1) - 2) D_S  (Keel-McKernan / Pandharipande; ptolemy.md §4).
F-curve of partition N1|N2|N3|N4: D_S.F = 1 if {S,S^c} = {Ni u Nj, Nk u Nl}; -1 if S = Ni or S^c = Ni (|Ni|>=2); else 0.
Also prints B.F (B = total boundary) and (K+B).F > 0 (ampleness sanity)."""
from fractions import Fraction
from itertools import combinations
n = 8
pts = list(range(n))
def norm(S):
    S = frozenset(S); C = frozenset(pts) - S
    return min(S, C, key=lambda X: (len(X), sorted(X)))
def coef(S):
    s = len(S); return Fraction(s*(n-s), n-1) - 2
def Fdot(parts, S):
    S = norm(S); tot = 0
    for i, j in [(0,1),(0,2),(0,3)]:
        U = parts[i] | parts[j]
        if norm(U) == S: tot += 1
    for P in parts:
        if len(P) >= 2 and norm(P) == S: tot -= 1
    return tot
splits = {norm(S) for r in range(2, n-1) for S in combinations(pts, r)}
assert len(splits) == 2**(n-1) - n - 1
types = [(5,1,1,1),(4,2,1,1),(3,3,1,1),(3,2,2,1),(2,2,2,2)]
for t in types:
    parts = []; c = 0
    for size in t:
        parts.append(frozenset(range(c, c+size))); c += size
    K = sum(coef(S) * Fdot(parts, S) for S in splits)
    beta = sum(Fdot(parts, S) for S in splits)
    KB = K + beta
    print(f"F-curve type {t}: K.F = {K}, B.F = {beta}, (K+B).F = {KB}")
