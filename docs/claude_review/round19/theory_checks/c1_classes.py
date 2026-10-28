"""c1: exact checks of the divisor-class facts used in theory.md (round 19).
All arithmetic exact (Fraction / int).
 (a) beta (D_T.beta = 1 for all 119 boundary divisors) respects every Keel relation on Mbar_{0,8}.
 (b) (K + boundary).beta_k on Mbar_{0,k}, pulled back to Mbar_{0,8}: degrees 16,40,68,98,129;
     successive differences Lambda_k.beta = 16,24,28,30,31; also psi_j.beta = 31 directly.
 (c) K.beta = psi.beta - 2*boundary.beta = 10 (K = psi - 2 delta in genus 0).
"""
from fractions import Fraction as Fr
from itertools import combinations
n = 8
pts = list(range(n))
def splits(n):
    S = []
    for s in range(2, n//2+1):
        for T in combinations(range(n), s):
            T = frozenset(T)
            Tc = frozenset(range(n)) - T
            if s == n - s and min(Tc) < min(T):
                continue
            S.append(T)
    return S
DIV = splits(8)
assert len(DIV) == 119
def canon(T, n=8):
    T = frozenset(T); Tc = frozenset(range(n)) - T
    if len(T) < len(Tc) or (len(T) == len(Tc) and min(T) < min(Tc)): return T
    return Tc
# (a) Keel relations: for distinct i,j,k,l: sum_{ij in T, kl notin T} D_T = sum_{ik in T, jl notin T} D_T
bad = 0; nrel = 0
for (i, j, k, l) in __import__('itertools').permutations(range(8), 4):
    lhs = rhs = 0
    for T in DIV:
        for TT in (T, frozenset(range(8)) - T):
            if i in TT and j in TT and k not in TT and l not in TT: lhs += 1
            if i in TT and k in TT and j not in TT and l not in TT: rhs += 1
    nrel += 1
    if lhs != rhs: bad += 1
print("(a) Keel relations checked:", nrel, "violations by beta:", bad, " (each side has beta-degree 16)")
# (b) (K+bdry).beta_k with K = sum_s (s(k-s)/(k-1) - 2) D_s
out = {}
for k in range(4, 9):
    divs = splits(k)
    Kb = sum(Fr(len(T)*(k-len(T)), k-1) - 2 for T in divs)
    db = len(divs)
    out[k] = (Kb, db, Kb + db, 2**(8-k)*(Kb+db))
    print(f"(b) k={k}: #div={db}, K.beta_k={Kb}, (K+d).beta_k={Kb+db}, pulled back to Mbar08: {2**(8-k)*(Kb+db)}")
prev = 0
for k in range(4, 9):
    v = out[k][3]
    print(f"    Lambda_{k}.beta = {v - prev}   (formula 32-2^(8-k) for k>=5: {32-2**(8-k)})")
    prev = v
# psi_j = sum over D_T with j in T, a,b notin T (a,b fixed != j); beta-degree
for j in range(8):
    a, b = [x for x in range(8) if x != j][:2]
    cnt = 0
    for T in DIV:
        for TT in (T, frozenset(range(8)) - T):
            if j in TT and a not in TT and b not in TT: cnt += 1
    assert cnt == 31
print("(b) psi_j.beta = 31 for all j (from psi_j = sum_{j in T; a,b notin T} D_T)")
Kb8 = out[8][0]
print("(c) K.beta =", Kb8, "; psi.beta - 2*119 =", 8*31 - 2*119)
